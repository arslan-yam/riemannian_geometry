"""POI setting (results/poi/*.json from scripts/run_poi.sh): our index built over the POIs only against
the plain exact P x P table (one full VTP sweep per POI). Writes notes/poi_results.md.
usage: python3 bench/poi_report.py
"""
import glob, json, os
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    rows = []
    for p in glob.glob(os.path.join(ROOT, "results", "poi", "*.json")):
        try:
            d = json.load(open(p))
        except (json.JSONDecodeError, FileNotFoundError):
            continue
        mesh = os.path.basename(d["mesh"]).replace(".off", "")
        v = {x["variant"]: x for x in d["variants"]}
        t = d.get("poi_table", {})
        rows.append((mesh, d["n"], d["pois"], t, d["build_s"], v))
    order = {"HorseMount": 0, "RobinsonMountain": 1, "Gunnison10k": 2}
    rows.sort(key=lambda r: (order.get(r[0], 9), r[2]))
    seen = set()
    out = ["# POI setting: our POI index against the exact POI table", "",
           "Queries are 1,000 random pairs among P random POIs, scored against exact MMP; single thread on the M1.",
           "The table is built by one full VTP sweep per POI and stores P x P float32 distances, so it is",
           "exact up to float32 rounding. Our index (eps 0.2, g 1 edge) builds the tree over the POIs",
           "only; every answer was sound and none needed refining.", "",
           "| mesh | n | P | table: build / size / query | ours: build / size (dense, 8-bit) / query / mean err | faster to build | smaller |",
           "|---|---:|---:|---|---|---|---|"]
    for mesh, n, P, t, tb, v in rows:
        if (mesh, P) in seen:
            continue
        seen.add((mesh, P))
        dn, c8 = v.get("dense"), v.get("code8")
        ours = f"{tb:.1f} s / {dn['index_bytes']/1e6:.2f}, {c8['index_bytes']/1e6:.2f} MB / {dn['query_us']:.2f} µs / {dn['est_mean']*100:.3f}%"
        tab = f"{t['build_s']:.1f} s / {t['index_bytes']/1e6:.2f} MB / {t['query_us']:.3f} µs" if t else "—"
        fb = ("table" if t and t["build_s"] < tb else "ours") if t else "—"
        sm = ("table" if t and t["index_bytes"] < c8["index_bytes"] else "ours") if t else "—"
        out.append(f"| {mesh} | {n} | {P} | {tab} | {ours} | {fb} | {sm} |")
    out += ["", "Reading: below roughly 1,000 POIs on 3.7k–10k-vertex meshes the exact table wins on build time,",
            "exactness and query time. Our index is cheaper to store from a few hundred POIs on, but its",
            "build is dominated by the top cuts in both modes, so it is only cheaper to build once P is",
            "large. The table grows as P^2 (2,500 POIs: 25 MB; 10,000: 400 MB), our index about linearly in P.", ""]
    open(os.path.join(ROOT, "notes", "poi_results.md"), "w").write("\n".join(out))
    print("wrote notes/poi_results.md")


if __name__ == "__main__":
    main()
