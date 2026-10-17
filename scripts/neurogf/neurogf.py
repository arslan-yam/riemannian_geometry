"""NeuroGF (Zhang et al., NeurIPS 2023, github.com/keeganhk/NeuroGF) on our terrains and our query sets.

Their per-model ("overfitting") network and training are copied verbatim from external/NeuroGF:
FC / MLP from cdbs/general.py, the models from code/neurogf_ovft/neurogf_ovft_models.py, and the training
loops and hyperparameters from train/pretraining_sdf_querier.py and train/training_neurogf.py. Their
package imports compiled CUDA extensions this code path never uses, so it is not imported. The training
data comes from cpp/tools/ngf_data: exact MMP distances and shortest paths, in their layout.

Changes from their scripts, each forced by the setting (nothing is tuned):
  * Signed distance: a terrain is an open surface, so the SDF head learns the unsigned distance to the
    surface (to a dense surface sample, through a k-d tree). It only enters through the path head's
    on-surface constraint (cstt_dcp), not through the distances.
  * Data loading: the training arrays are held in memory (their scripts read one .npy per source per
    epoch). The per-epoch sampling is theirs: 256 sources x 256 / 32 / 64 targets, and 128 sources x 32
    global / 32 local paths, then 200 iterations on that batch.
  * separate_refining.py is not run. It refines the SDF and path heads with the lifting frozen, so the
    distance queries do not change.
  * The best checkpoint is chosen, as in theirs, by the mean relative error on held-out validation
    sources, never on the test queries. Test-pair vertices are excluded from all training data.
Coordinates are normalised into [-1, 1]^3 (their SDF sampling assumes it): centre of the bounding box,
divided by half its largest side. Distances are scaled the same way, and predictions are scaled back.

usage: python3 neurogf.py --mesh m.off --data DIR --queries PREFIX --out DIR [--epochs 400] [--sdf-epochs 100]
       PREFIX: the exported query files (PREFIX.V2V.pairs / .truth / .detour, PREFIX.A2A.query / .truth)
writes OUT/neurogf.json (rewritten as it goes, "complete" at the end), the checkpoints and the predictions
"""
import argparse, json, os, sys, time
# thread pools sized to the cores this job may use (a Slurm job sees every core of the node, but runs on its
# allocation: 48 BLAS / OpenMP / k-d tree threads on one core thrash). Set before numpy and torch load.
NCPU = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ[_v] = str(NCPU)
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from scipy.spatial import cKDTree
torch.set_num_threads(NCPU)

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
from rgeo.terrain import read_off


# ---------------------------------------------------------------- their layers and models (verbatim)
class MLP(nn.Module):
    def __init__(self, ic, oc, is_bn, nl):
        super(MLP, self).__init__()
        assert isinstance(is_bn, bool)
        assert nl in ['none', 'relu', 'tanh', 'sigmoid']
        self.is_bn = is_bn
        self.nl = nl
        self.conv = nn.Conv2d(in_channels=ic, out_channels=oc, kernel_size=1, bias=False)
        if self.is_bn:
            self.bn = nn.BatchNorm2d(oc)
        if nl == 'relu':
            self.activate = nn.ReLU(inplace=True)
        if nl == 'tanh':
            self.activate = nn.Tanh()
        if nl == 'sigmoid':
            self.activate = nn.Sigmoid()
    def forward(self, x):
        x = x.permute(0, 2, 1).contiguous().unsqueeze(-1)
        y = self.conv(x)
        if self.is_bn:
            y = self.bn(y)
        if self.nl != 'none':
            y = self.activate(y)
        y = y.squeeze(-1).permute(0, 2, 1).contiguous()
        return y


class FC(nn.Module):
    def __init__(self, ic, oc, is_bn, nl):
        super(FC, self).__init__()
        assert isinstance(is_bn, bool)
        assert nl in ['none', 'relu', 'tanh', 'sigmoid']
        self.is_bn = is_bn
        self.nl = nl
        self.linear = nn.Linear(ic, oc, bias=False)
        if self.is_bn:
            self.bn = nn.BatchNorm1d(oc)
        if nl == 'relu':
            self.activate = nn.ReLU(inplace=True)
        if nl == 'tanh':
            self.activate = nn.Tanh()
        if nl == 'sigmoid':
            self.activate = nn.Sigmoid()
    def forward(self, x):
        y = self.linear(x)
        if self.is_bn:
            y = self.bn(y)
        if self.nl != 'none':
            y = self.activate(y)
        return y


class SDist_Querier_Offline_Pretrainer(nn.Module):
    def __init__(self):
        super(SDist_Querier_Offline_Pretrainer, self).__init__()
        fc_1 = FC(3, 32, True, 'relu')
        fc_2 = FC(32, 64, True, 'relu')
        fc_3 = FC(64, 128, True, 'relu')
        fc_4 = FC(128, 256, True, 'relu')
        fc_5 = FC(256, 512, True, 'relu')
        fc_6 = FC(512, 128, True, 'relu')
        fc_7 = FC(128, 1, False, 'none')
        self.fc = nn.Sequential(fc_1, fc_2, fc_3, fc_4, fc_5, fc_6, fc_7)
    def forward(self, sdf_q):
        return self.fc(sdf_q).squeeze(-1)


class NeuroGF_Offline_Trainer(nn.Module):
    def __init__(self, D, num_pp_g, num_pp_l):
        super(NeuroGF_Offline_Trainer, self).__init__()
        self.num_pp_g = num_pp_g
        self.num_pp_l = num_pp_l
        self.lifting = nn.Sequential(FC(3, D//4, True, 'relu'), FC(D//4, D//2, True, 'relu'), FC(D//2, D, True, 'relu'))
        self.sdf_head = nn.Sequential(FC(D, D//4, True, 'relu'), FC(D//4, 64, True, 'relu'), FC(64, 1, False, 'none'))
        self.gdf_head = nn.Sequential(FC(D, D//4, True, 'relu'), FC(D//4, 64, True, 'relu'), FC(64, 1, False, 'none'))
        self.cdw_fc = FC(2*D, D, True, 'relu')
        mlp_1 = MLP(D+3, 128, True, 'relu')
        mlp_2 = MLP(128, 64, True, 'relu')
        mlp_3 = MLP(64, 32, True, 'relu')
        mlp_4 = MLP(32, 3, False, 'none')
        self.spf_head = nn.Sequential(mlp_1, mlp_2, mlp_3, mlp_4)
    def forward(self, sdf_q, gdf_q_g, gdf_q_l1, gdf_q_l2, spf_q_g, spf_q_l, lines_g, lines_l):
        N = sdf_q.size(0)
        Mg, Ml1, Ml2 = gdf_q_g.size(0), gdf_q_l1.size(0), gdf_q_l2.size(0)
        Kg, Kl = spf_q_g.size(0), spf_q_l.size(0)
        assert lines_g.size(0)==Kg and lines_l.size(0)==Kl
        gdf_g_cats = torch.cat((gdf_q_g[:, 0:3], gdf_q_g[:, 3:6]), dim=0)
        gdf_l1_cats = torch.cat((gdf_q_l1[:, 0:3], gdf_q_l1[:, 3:6]), dim=0)
        gdf_l2_cats = torch.cat((gdf_q_l2[:, 0:3], gdf_q_l2[:, 3:6]), dim=0)
        spf_g_cats = torch.cat((spf_q_g[:, 0:3], spf_q_g[:, 3:6]), dim=0)
        spf_l_cats = torch.cat((spf_q_l[:, 0:3], spf_q_l[:, 3:6]), dim=0)
        all_cats = torch.cat((sdf_q, gdf_g_cats, gdf_l1_cats, gdf_l2_cats, spf_g_cats, spf_l_cats), dim=0)
        ftr__all = self.lifting(all_cats)
        ftr__sdf_q = ftr__all[0:N]
        ftr__gdf_g_cats = ftr__all[(N):(N+2*Mg)]
        ftr__gdf_l1_cats = ftr__all[(N+2*Mg):(N+2*Mg+2*Ml1)]
        ftr__gdf_l2_cats = ftr__all[(N+2*Mg+2*Ml1):(N+2*Mg+2*Ml1+2*Ml2)]
        ftr__spf_g_cats = ftr__all[(N+2*Mg+2*Ml1+2*Ml2):(N+2*Mg+2*Ml1+2*Ml2+2*Kg)]
        ftr__spf_l_cats = ftr__all[(N+2*Mg+2*Ml1+2*Ml2+2*Kg):(N+2*Mg+2*Ml1+2*Ml2+2*Kg+2*Kl)]
        sdist_out = self.sdf_head(ftr__sdf_q).squeeze(-1)
        diff_g = (ftr__gdf_g_cats[:Mg] - ftr__gdf_g_cats[Mg:]).abs()
        diff_l1 = (ftr__gdf_l1_cats[:Ml1] - ftr__gdf_l1_cats[Ml1:]).abs()
        diff_l2 = (ftr__gdf_l2_cats[:Ml2] - ftr__gdf_l2_cats[Ml2:]).abs()
        gdist_out_cats = self.gdf_head(torch.cat((diff_g, diff_l1, diff_l2), dim=0)).squeeze(-1)
        gdist_out_g = gdist_out_cats[0:Mg]
        gdist_out_l1 = gdist_out_cats[(Mg):(Mg+Ml1)]
        gdist_out_l2 = gdist_out_cats[(Mg+Ml1):(Mg+Ml1+Ml2)]
        a = torch.cat((ftr__spf_g_cats[:Kg], ftr__spf_g_cats[Kg:]), dim=-1)
        b = torch.cat((ftr__spf_l_cats[:Kl], ftr__spf_l_cats[Kl:]), dim=-1)
        cdw__spf_cats = self.cdw_fc(torch.cat((a, b), dim=0))
        cdw_g_rep = cdw__spf_cats[:Kg].unsqueeze(1).repeat(1, self.num_pp_g, 1)
        cdw_l_rep = cdw__spf_cats[Kg:].unsqueeze(1).repeat(1, self.num_pp_l, 1)
        spath_out_g = self.spf_head(torch.cat((lines_g, cdw_g_rep), dim=-1))
        spath_out_l = self.spf_head(torch.cat((lines_l, cdw_l_rep), dim=-1))
        return sdist_out, [gdist_out_g, gdist_out_l1, gdist_out_l2], [spath_out_g, spath_out_l]


class NeuroGF_Online_Querier_GDistOnly(nn.Module):
    # their distance querier: the trained network's lifting + distance head (the other heads are kept so
    # the trainer's state dict loads as in their scripts)
    def __init__(self, D):
        super(NeuroGF_Online_Querier_GDistOnly, self).__init__()
        self.lifting = nn.Sequential(FC(3, D//4, True, 'relu'), FC(D//4, D//2, True, 'relu'), FC(D//2, D, True, 'relu'))
        self.sdf_head = nn.Sequential(FC(D, D//4, True, 'relu'), FC(D//4, 64, True, 'relu'), FC(64, 1, False, 'none'))
        self.gdf_head = nn.Sequential(FC(D, D//4, True, 'relu'), FC(D//4, 64, True, 'relu'), FC(64, 1, False, 'none'))
        self.cdw_fc = FC(2*D, D, True, 'relu')
        self.spf_head = nn.Sequential(MLP(D+3, 128, True, 'relu'), MLP(128, 64, True, 'relu'), MLP(64, 32, True, 'relu'), MLP(32, 3, False, 'none'))
    def forward(self, gdf_q):
        N = gdf_q.size(0)
        with torch.no_grad():
            ftr__all = self.lifting(torch.cat((gdf_q[:, 0:3], gdf_q[:, 3:6]), dim=0))
        diff = (ftr__all[:N] - ftr__all[N:]).abs()
        return self.gdf_head(diff).squeeze(-1)


# their helpers (cdbs/custom.py, cdbs/general.py)
def approx_path_len_batched(path_pts_batched):
    B, N, _ = path_pts_batched.size()
    x = path_pts_batched[:, 0:N-1, :]
    y = path_pts_batched[:, 1:N, :]
    return ((((x-y)**2).sum(dim=-1))**0.5).sum(dim=-1)


def generate_lines_from_end_points(ps, pe, n):
    t = np.linspace(0, 1, n, dtype=np.float32)[None, :, None]
    return (ps[:, None, :] * (1 - t) + pe[:, None, :] * t).astype(np.float32)  # = their per-row np.linspace


def random_jittering(pc, sigma, bound):
    return pc + np.clip(np.random.normal(0, sigma, size=pc.shape), -bound, bound).astype(np.float32)


# ---------------------------------------------------------------- data
ap = argparse.ArgumentParser()
ap.add_argument("--mesh", required=True); ap.add_argument("--data", required=True)
ap.add_argument("--queries", required=True); ap.add_argument("--out", required=True)
ap.add_argument("--epochs", type=int, default=400); ap.add_argument("--sdf-epochs", type=int, default=100)
ap.add_argument("--iters", type=int, default=200)
ap.add_argument("--surface-samples", type=int, default=2_000_000)
ap.add_argument("--test-scale", type=int, default=1)  # local smoke tests only: divides the sources per epoch
a = ap.parse_args()
os.makedirs(a.out, exist_ok=True)
dev = torch.device("cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu")
result = {"tool": "neurogf", "argv": " ".join(sys.argv), "mesh": a.mesh, "device": str(dev),
          "gpu": torch.cuda.get_device_name(0) if dev.type == "cuda" else None, "torch": torch.__version__, "cpus": NCPU,
          "epochs": a.epochs, "sdf_epochs": a.sdf_epochs, "iters": a.iters, "test_scale": a.test_scale, "complete": False}
def flush():
    json.dump(result, open(os.path.join(a.out, "neurogf.json"), "w"), indent=1)

V, F = read_off(a.mesh)
V = np.asarray(V, dtype=np.float64); F = np.asarray(F, dtype=np.int64)
lo, hi = V.min(0), V.max(0)
ctr = 0.5 * (lo + hi); scale = 0.5 * float((hi - lo).max())
norm = lambda x: ((np.asarray(x, dtype=np.float64) - ctr) / scale).astype(np.float32)
mv = norm(V)
result.update(n=int(len(V)), faces=int(len(F)), scale=scale)
# dense surface sample for the (unsigned) distance to the surface
P = mv[F]
area = 0.5 * np.linalg.norm(np.cross(P[:, 1] - P[:, 0], P[:, 2] - P[:, 0]), axis=1)
rng = np.random.default_rng(0)
fi = rng.choice(len(F), a.surface_samples, p=area / area.sum())
r1 = np.sqrt(rng.random(a.surface_samples, dtype=np.float32))[:, None]; r2 = rng.random(a.surface_samples, dtype=np.float32)[:, None]
sp = (P[fi, 0] * (1 - r1) + P[fi, 1] * r1 * (1 - r2) + P[fi, 2] * r1 * r2).astype(np.float32)
tree = cKDTree(sp)

def fetch_train_data_sdf(num_fetch):
    k = num_fetch // 3
    q1 = ((np.random.rand(k, 3) - 0.5) * 2 * 1.01).astype(np.float32)
    q2 = sp[np.random.choice(len(sp), k, replace=False)]
    q3 = random_jittering(sp[np.random.choice(len(sp), k, replace=False)], 0.03, 0.03)
    q = np.concatenate((q1, q2, q3), axis=0).astype(np.float32)
    return q, tree.query(q, workers=NCPU)[0].astype(np.float32)

ld = lambda f: np.load(os.path.join(a.data, f + ".npy"))
def nrows(x):  # (..., 7) rows: coordinates and the distance normalised
    x = x.astype(np.float32).copy()
    x[..., 0:3] = norm(x[..., 0:3]); x[..., 3:6] = norm(x[..., 3:6]); x[..., 6] /= scale
    return x
G = nrows(ld("train_global"))                                  # (S, G, 7)
L1, O1 = nrows(ld("train_l1")), ld("train_l1_off")
L2, O2 = nrows(ld("train_l2")), ld("train_l2_off")
PG = norm(ld("train_paths_global"))                            # (S, P, 128, 3)
PL, OP = norm(ld("train_paths_local")), ld("train_paths_local_off")
VAL = nrows(ld("val"))
S = G.shape[0]
meta = json.load(open(os.path.join(a.data, "meta.json")))
result["data"] = meta

def pick(n, k):
    return np.random.choice(n, k, replace=n < k)

def fetch_train_data_gdf(num_s, nt_g, nt_l1, nt_l2):
    srcs = np.random.choice(S, num_s)
    g = np.stack([G[s][pick(G.shape[1], nt_g)] for s in srcs]).reshape(-1, 7)
    l1 = np.stack([L1[O1[s] + pick(O1[s + 1] - O1[s], nt_l1)] for s in srcs]).reshape(-1, 7)
    l2 = np.stack([L2[O2[s] + pick(O2[s + 1] - O2[s], nt_l2)] for s in srcs]).reshape(-1, 7)
    return g[:, 0:6], l1[:, 0:6], l2[:, 0:6], g[:, 6], l1[:, 6], l2[:, 6]

def fetch_train_data_spf(num_s, nt_g, nt_l):
    srcs = np.random.choice(S, num_s)
    g = np.stack([PG[s][pick(PG.shape[1], nt_g)] for s in srcs]).reshape(-1, PG.shape[2], 3)
    l = np.stack([PL[OP[s] + pick(OP[s + 1] - OP[s], nt_l)] for s in srcs]).reshape(-1, PL.shape[1], 3)
    return g, l

T = lambda x: torch.tensor(np.ascontiguousarray(x)).float().to(dev)

# ---------------------------------------------------------------- 1. distance-to-surface querier (pretraining_sdf_querier.py)
t0 = time.time()
sdist_querier = SDist_Querier_Offline_Pretrainer().to(dev)
optimizer = optim.AdamW(sdist_querier.parameters(), lr=1e-2, weight_decay=1e-8)
scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=a.sdf_epochs, eta_min=1e-4)
criterion = nn.MSELoss()
best_loss, best_state = 1e8, None
prep_s = 0.0
for epc in range(1, a.sdf_epochs + 1):
    sdist_querier.train()
    tp = time.time()
    q, gt = fetch_train_data_sdf(90000); q, gt = T(q), T(gt)
    prep_s += time.time() - tp
    for _ in range(a.iters):
        optimizer.zero_grad()
        loss = criterion(sdist_querier(q), gt)
        loss.backward(); optimizer.step()
    scheduler.step()
    if loss.item() < best_loss:
        best_loss = loss.item(); best_state = {k: v.detach().cpu().clone() for k, v in sdist_querier.state_dict().items()}
sdist_querier.load_state_dict(best_state); sdist_querier.eval()
for p in sdist_querier.parameters(): p.requires_grad_(False)
torch.save(best_state, os.path.join(a.out, "sdf_querier.pth"))
result["sdf_pretrain_s"] = time.time() - t0; result["sdf_prep_s"] = prep_s; flush()
print(f"distance-to-surface querier: {result['sdf_pretrain_s']:.0f} s ({prep_s:.0f} s preparing data), loss {best_loss:.3g}", flush=True)

# ---------------------------------------------------------------- 2. NeuroGF (training_neurogf.py)
D, num_pp_g, num_pp_l = 256, PG.shape[2], PL.shape[1]
net = NeuroGF_Offline_Trainer(D, num_pp_g, num_pp_l).to(dev)
optimizer = optim.AdamW(net.parameters(), lr=1e-2, weight_decay=1e-8)
scheduler = optim.lr_scheduler.CosineAnnealingLR(optimizer, T_max=a.epochs, eta_min=1e-4)
criterion = nn.L1Loss()

def querier_from(state):
    q = NeuroGF_Online_Querier_GDistOnly(D); q.load_state_dict(state); return q.to(dev).eval()

def predict(q, X, bs=100000):
    out = []
    with torch.no_grad():
        for i in range(0, len(X), bs): out.append(q(T(X[i:i + bs])).cpu().numpy())
    return np.concatenate(out) if out else np.zeros(0, np.float32)

t1 = time.time()
best_mre, best_state, log = 1e8, None, []
prep_s = 0.0
for epc in range(1, a.epochs + 1):
    net.train()
    tp = time.time()
    sdf_q, sdist_gt = fetch_train_data_sdf(30000)
    gq_g, gq_l1, gq_l2, gd_g, gd_l1, gd_l2 = fetch_train_data_gdf(256 // a.test_scale, 256, 32, 64)
    pg_uni, pl_uni = fetch_train_data_spf(128 // a.test_scale, 32, 32)
    spath_gt_g = np.concatenate((pg_uni, np.flip(pg_uni, axis=1)), axis=0)
    spath_gt_l = np.concatenate((pl_uni, np.flip(pl_uni, axis=1)), axis=0)
    spf_q_g = np.concatenate((spath_gt_g[:, 0, :], spath_gt_g[:, -1, :]), axis=-1)
    spf_q_l = np.concatenate((spath_gt_l[:, 0, :], spath_gt_l[:, -1, :]), axis=-1)
    lines_g = generate_lines_from_end_points(spf_q_g[:, 0:3], spf_q_g[:, 3:6], num_pp_g)
    lines_l = generate_lines_from_end_points(spf_q_l[:, 0:3], spf_q_l[:, 3:6], num_pp_l)
    sdf_q, sdist_gt = T(sdf_q), T(sdist_gt)
    gq_g, gq_l1, gq_l2, gd_g, gd_l1, gd_l2 = map(T, (gq_g, gq_l1, gq_l2, gd_g, gd_l1, gd_l2))
    spf_q_g, spf_q_l, spath_gt_g, spath_gt_l, lines_g, lines_l = map(T, (spf_q_g, spf_q_l, spath_gt_g, spath_gt_l, lines_g, lines_l))
    prep_s += time.time() - tp
    for _ in range(a.iters):
        optimizer.zero_grad()
        sdist_out, (go_g, go_l1, go_l2), (so_g, so_l) = net(sdf_q, gq_g, gq_l1, gq_l2, spf_q_g, spf_q_l, lines_g, lines_l)
        loss_sdist = criterion(sdist_out, sdist_gt)
        loss_gdist = criterion(go_g, gd_g) + criterion(go_l1, gd_l1) + criterion(go_l2, gd_l2)
        loss_spath = criterion(so_g, spath_gt_g) + criterion(so_l, spath_gt_l)
        cstt_ccl = criterion(approx_path_len_batched(so_g), approx_path_len_batched(spath_gt_g))
        idx = np.random.choice(so_g.shape[0], min(1024, so_g.shape[0]), replace=False)
        cstt_dcp = sdist_querier(so_g[idx].reshape(-1, 3)).abs().mean()
        loss = 0.1*loss_sdist + 1.0*loss_gdist + 0.1*loss_spath + 0.1*cstt_ccl + 0.1*cstt_dcp
        loss.backward(); optimizer.step()
    scheduler.step()
    if epc % 50 == 0 or (epc >= a.epochs // 2 and epc % 10 == 0) or epc == a.epochs:
        state = {k: v.detach().cpu().clone() for k, v in net.state_dict().items()}
        pv = predict(querier_from(state), VAL[:, 0:6])
        mre = float((np.abs(pv - VAL[:, 6]) / VAL[:, 6]).mean())
        log.append({"epoch": epc, "val_mre": mre, "loss_gdist": float(loss_gdist.item()), "s": time.time() - t1, "prep_s": prep_s})
        if mre < best_mre:
            best_mre, best_state = mre, state
            torch.save(best_state, os.path.join(a.out, "neurogf_best.pth"))
        print(f"epoch {epc}: val MRE {100 * mre:.3f}% (best {100 * best_mre:.3f}%), {time.time() - t1:.0f} s ({prep_s:.0f} s preparing data)", flush=True)
        result.update(train_s=time.time() - t1, val_mre=best_mre, log=log); flush()
result.update(train_s=time.time() - t1, val_mre=best_mre, log=log)

# ---------------------------------------------------------------- 3. evaluation on our query sets
q = querier_from(best_state)
nparams = lambda m: sum(p.numel() for p in m.parameters())
result["model_bytes"] = 4 * nparams(net)                                  # the whole trained network
result["index_bytes"] = 4 * (nparams(q.lifting) + nparams(q.gdf_head))    # what a distance query needs

def stats(pred, truth, tag):
    e = np.abs(pred - truth) / truth
    s = (pred - truth) / truth
    return {"set": tag, "pairs": int(len(e)), "err_mean": float(e.mean()), "err_p50": float(np.median(e)),
            "err_p99": float(np.quantile(e, 0.99)), "err_max": float(e.max()), "err_mean_signed": float(s.mean()),
            "over_eps_0.2": float((e > 0.2).mean())}

sets = []
pairs_f = a.queries + ".V2V.pairs"
if os.path.exists(pairs_f):
    pr = np.loadtxt(pairs_f, dtype=np.int64).reshape(-1, 2); tv = np.loadtxt(a.queries + ".V2V.truth").reshape(-1)
    Xv = np.concatenate((mv[pr[:, 0]], mv[pr[:, 1]]), axis=1)
    pv = predict(q, Xv) * scale
    np.save(os.path.join(a.out, "pred_V2V.npy"), pv)
    sets.append(stats(pv, tv, "V2V"))
    if os.path.exists(a.queries + ".V2V.detour"):
        r = np.loadtxt(a.queries + ".V2V.detour").reshape(-1)
        cut = [1.0, 1.01, 1.05, 1.2, 1e9]
        for b in range(4):
            m = (r >= cut[b]) & (r < cut[b + 1])
            if m.any(): sets.append(dict(stats(pv[m], tv[m], f"detour{b}"), lo=cut[b], hi=cut[b + 1]))
    # query time: GPU batched (all pairs in one call), and single queries on one CPU thread
    Xt = T(Xv)
    with torch.no_grad():
        q(Xt)
        if dev.type == "cuda": torch.cuda.synchronize()
        t = time.time(); reps = 20
        for _ in range(reps): q(Xt)
        if dev.type == "cuda": torch.cuda.synchronize()
        result["query_us_batch"] = 1e6 * (time.time() - t) / (reps * len(Xv))
    qc = querier_from(best_state).to("cpu").eval()
    torch.set_num_threads(1)
    Xc = torch.tensor(Xv[:200]).float()
    with torch.no_grad():
        for i in range(10): qc(Xc[i:i + 1])
        t = time.time()
        for i in range(len(Xc)): qc(Xc[i:i + 1])
        result["query_us"] = 1e6 * (time.time() - t) / len(Xc)
        t = time.time(); qc(torch.tensor(Xv).float()); result["query_us_cpu_batch"] = 1e6 * (time.time() - t) / len(Xv)
if os.path.exists(a.queries + ".A2A.query"):
    Q = np.loadtxt(a.queries + ".A2A.query").reshape(-1, 8); ta = np.loadtxt(a.queries + ".A2A.truth").reshape(-1)
    Xa = np.concatenate((norm(Q[:, 0:3]), norm(Q[:, 4:7])), axis=1)
    pa = predict(q, Xa) * scale
    np.save(os.path.join(a.out, "pred_A2A.npy"), pa)
    sets.append(stats(pa, ta, "A2A"))
result["sets"] = sets
result["build_s"] = result["sdf_pretrain_s"] + result["train_s"]   # GPU time; the CPU data generation is result["data"]["cpu_s"]
result["complete"] = True
flush()
for s_ in sets:
    print(f"{s_['set']}: {s_['pairs']} pairs, mean {100 * s_['err_mean']:.3f}% p99 {100 * s_['err_p99']:.2f}% max {100 * s_['err_max']:.2f}%")
print(f"query {result.get('query_us', float('nan')):.1f} us (1 CPU thread), {result.get('query_us_batch', float('nan')):.3f} us batched on {dev}; "
      f"model {result['index_bytes'] / 1e6:.2f} MB for distances; training {result['build_s']:.0f} s")
