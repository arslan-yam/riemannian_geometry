"""NeuroGF (Zhang, Hou, Adikusuma, Wang, He, NeurIPS 2023; github.com/keeganhk/NeuroGF) on our terrains.

Per-model ("overfitting") mode, geodesic-distance branch: the authors' network blocks (FC = Linear +
BatchNorm + ReLU, lifting 3 -> D/4 -> D/2 -> D, distance head on |f(p) - f(q)|, D = 256), their loss
(L1 on distances, global + two local target sets per source, weights 1/1/1), optimiser (AdamW,
lr 1e-2 -> 1e-4 cosine) and sampling sizes (256 sources x (256 global + 32 + 64 local) targets per
iteration). Deviations, all forced by the setting: the training data are generated here (their
prepared data is for their 10 shape models) with exact MMP sweeps from `--sources` random vertices;
the SDF and shortest-path branches are dropped (they serve path queries; in their network the
branches share only the lifting MLP); the epoch budget is reduced from 400 x 200 iterations to
`--epochs` x `--iters` so a model trains in minutes on an Apple M1 (MPS); vertices only (V2V queries).

Reported: build = data generation (exact sweeps) + training, size = parameters (float32), query
time single-pair on CPU and batched, and relative errors on our 1,500 seed-7 pairs (results/pairs).
usage: python3 scripts/competitors/neurogf_terrain.py MESH [--sources 512] [--epochs 60] [--iters 200]
"""
import argparse, json, os, sys, time
import numpy as np
import torch
import torch.nn as nn
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, ROOT)
from rgeo.terrain import read_off
import pygeodesic.geodesic as pg


class FC(nn.Module): #cdbs/general.py of the NeuroGF code
    def __init__(self, ic, oc, is_bn, nl):
        super().__init__()
        self.linear = nn.Linear(ic, oc, bias=False)
        self.bn = nn.BatchNorm1d(oc) if is_bn else None
        self.act = nn.ReLU(inplace=True) if nl == "relu" else None
    def forward(self, x):
        y = self.linear(x)
        if self.bn is not None:
            y = self.bn(y)
        return self.act(y) if self.act is not None else y


class NeuroGFDistance(nn.Module): #lifting + gdf_head of NeuroGF_Offline_Trainer / Online_Querier_GDistOnly
    def __init__(self, D=256):
        super().__init__()
        self.lifting = nn.Sequential(FC(3, D // 4, True, "relu"), FC(D // 4, D // 2, True, "relu"), FC(D // 2, D, True, "relu"))
        self.gdf_head = nn.Sequential(FC(D, D // 4, True, "relu"), FC(D // 4, 64, True, "relu"), FC(64, 1, False, "none"))
    def forward(self, qs, qe):
        f = self.lifting(torch.cat((qs, qe), 0))
        n = qs.shape[0]
        return self.gdf_head((f[:n] - f[n:]).abs()).squeeze(-1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("mesh")
    ap.add_argument("--sources", type=int, default=512)
    ap.add_argument("--epochs", type=int, default=60)
    ap.add_argument("--iters", type=int, default=200)
    ap.add_argument("--device", default="mps" if torch.backends.mps.is_available() else "cpu")
    a = ap.parse_args()
    V, F = read_off(os.path.join(ROOT, "data", "terrains", f"{a.mesh}.off"))
    n = len(V)
    #normalise into [-1, 1]^3 with one scale (their meshes are normalised the same way)
    ctr = 0.5 * (V.max(0) + V.min(0)); scale = 0.5 * float((V.max(0) - V.min(0)).max())
    Vn = ((V - ctr) / scale).astype(np.float32)
    rng = np.random.default_rng(0)
    t0 = time.perf_counter()
    geo = pg.PyGeodesicAlgorithmExact(V, F)
    src = rng.choice(n, size=min(a.sources, n), replace=False)
    D = np.zeros((len(src), n), dtype=np.float32)
    for i, s in enumerate(src):
        D[i] = geo.geodesicDistances(np.array([s]), None)[0] / scale
    t_data = time.perf_counter() - t0
    #local target sets per source: the 32 and the 5% nearest vertices (their local_1 / local_2)
    d2 = ((Vn[src][:, None, :2] - Vn[None, :, :2]) ** 2).sum(-1)
    order = np.argsort(d2, 1)
    loc1 = order[:, 1:33]; loc2 = order[:, 1:max(65, n // 20)]
    dev = torch.device(a.device)
    net = NeuroGFDistance().to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=1e-2, weight_decay=1e-8)
    sch = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=a.epochs, eta_min=1e-4)
    L1 = nn.L1Loss()
    Vt = torch.tensor(Vn, device=dev)
    def batch(table, k): #table: (sources, m) target ids per source, or None for all vertices
        si = rng.integers(len(src), size=256)
        if table is None:
            t = rng.integers(n, size=(256, k))
        else:
            t = table[si[:, None], rng.integers(table.shape[1], size=(256, k))]
        s = np.repeat(src[si], k); t = t.ravel()
        d = D[np.repeat(si, k), t]
        return Vt[torch.as_tensor(s, device=dev)], Vt[torch.as_tensor(t, device=dev)], torch.as_tensor(d, device=dev)
    t1 = time.perf_counter()
    for ep in range(a.epochs):
        net.train()
        for _ in range(a.iters):
            #one forward pass over the global and both local sets together, as in their training
            #script (BatchNorm statistics then see the same mix in every batch), one L1 term per set
            sets = [batch(table, k) for table, k in ((None, 256), (loc1, 32), (loc2, 64))]
            out = net(torch.cat([x[0] for x in sets]), torch.cat([x[1] for x in sets]))
            loss, o = 0, 0
            for qs, _, dg in sets:
                loss = loss + L1(out[o:o + len(dg)], dg); o += len(dg)
            opt.zero_grad(); loss.backward(); opt.step()
        sch.step()
        if (ep + 1) % 10 == 0:
            print(f"  epoch {ep + 1}/{a.epochs}: loss {loss.item():.6f} ({time.perf_counter() - t1:.0f} s)", flush=True)
    t_train = time.perf_counter() - t1
    #evaluation on the same 1,500 pairs as the index (truth from the C++ run's cache)
    pairs = np.loadtxt(os.path.join(ROOT, "results", "pairs", f"{a.mesh}.seed7.txt"), dtype=int)
    truth = np.loadtxt(os.path.join(ROOT, "results", "cpp", f"{a.mesh}.seed7.truth"))
    net.eval().to("cpu")
    Vc = torch.tensor(Vn)
    with torch.no_grad():
        t2 = time.perf_counter(); est = net(Vc[pairs[:, 0]], Vc[pairs[:, 1]]).numpy() * scale
        t_batch = (time.perf_counter() - t2) / len(pairs)
        t3 = time.perf_counter()
        for s, t in pairs[:300]:
            net(Vc[[s]], Vc[[t]])
        t_single = (time.perf_counter() - t3) / 300
    e = est / truth - 1
    out = dict(method="NeuroGF (distance branch)", mesh=a.mesh, n=int(n), sources=int(len(src)), epochs=a.epochs, iters=a.iters,
               build_s=t_data + t_train, data_s=t_data, train_s=t_train, index_MB=sum(p.numel() for p in net.parameters()) * 4 / 1e6,
               query_us=t_single * 1e6, query_us_batch=t_batch * 1e6, err_mean=float(np.abs(e).mean() * 100),
               err_signed_mean=float(e.mean() * 100), err_max=float(np.abs(e).max() * 100), err_p99=float(np.percentile(np.abs(e), 99) * 100))
    print(json.dumps(out))
    os.makedirs(os.path.join(ROOT, "results", "neurogf"), exist_ok=True)
    json.dump(out, open(os.path.join(ROOT, "results", "neurogf", f"{a.mesh}.json"), "w"), indent=1)


if __name__ == "__main__":
    main()
