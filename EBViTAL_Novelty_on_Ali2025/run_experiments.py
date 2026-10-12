"""Fairness check + all EB-ViTAL experiments -> results/fairness.json, results/runs.json."""
import json
from multiprocessing import Pool

import numpy as np

from wsn_sim import simulate

NS = (100, 200, 300, 400, 500)
SEEDS = range(30)
KEYS = ("FND", "HND", "LND", "packets", "res500", "res1000")
ABL = {"EB-ViTAL": {}, "EB -M1": {"m1": False}, "EB -M2": {"m2": False}, "EB -M3": {"m3": False}}


def init():
    global DEP
    DEP = dict(np.load("results/deployments.npz"))


def job(a):
    label, proto, kind, N, s, kw = a
    o = simulate(proto, N, s, xy=DEP[f"{kind}_{N}_{s}"], **kw)
    return dict(label=label, protocol=proto, deploy=kind, N=N, seed=s,
                L=kw.get("L"), mobile=kw.get("mobile", False)) | o


def old_job(r):
    kw = {k: r[k] for k in ("L", "mobile", "dist_inverted") if k in r}
    if "w" in r:
        kw["w"] = tuple(r["w"])
    o = simulate(r["protocol"], r["N"], r["seed"], xy=DEP[f"{r['deploy']}_{r['N']}_{r['seed']}"], **kw)
    return all(o[k] == r[k] for k in KEYS)


if __name__ == "__main__":
    old = json.load(open("results/baseline_runs_from_reconstruction.json"))
    with Pool(64, initializer=init) as pool:
        same = pool.map(old_job, old, chunksize=1)
        json.dump(dict(runs=len(old), identical=int(sum(same))), open("results/fairness.json", "w"))
        print("fairness:", sum(same), "/", len(old), "identical")
        assert all(same), "copied simulator no longer reproduces the reconstruction"

        jobs = []
        for L in (800, 4000):
            base = [("LEACH (random)", "LEACH", "random", {}), ("LEACH (HSA)", "LEACH", "hsa", {}),
                    ("ViTAL", "ViTAL", "hsa", {}), ("ViTAL-noCDS", "4F", "hsa", {}),
                    ("EB-ViTAL (random)", "EB-ViTAL", "random", {})]
            base += [(lab, "EB-ViTAL", "hsa", {"eb": e}) for lab, e in ABL.items()]
            jobs += [(lab, p, k, N, s, dict(L=L) | kw) for lab, p, k, kw in base for N in NS for s in SEEDS]
            mob = [("LEACH (HSA)", "LEACH", "hsa"), ("ViTAL", "ViTAL", "hsa"), ("EB-ViTAL", "EB-ViTAL", "hsa")]
            jobs += [(lab, p, k, N, s, dict(L=L, mobile=True)) for lab, p, k in mob for N in NS for s in SEEDS]
        rows = pool.map(job, jobs, chunksize=1)
    json.dump(rows, open("results/runs.json", "w"))
    print("done", len(rows))
