"""All runs for the ViTAL reproduction -> results/*.json, results/deployments.npz."""
import json
import os
from multiprocessing import Pool

import numpy as np
from scipy.spatial import cKDTree

from vital_repro import coverage, deploy, simulate

NS = (100, 200, 300, 400, 500)
SEEDS = range(30)
OUT = "results"
DEP = f"{OUT}/deployments.npz"
STATIC = [("LEACH", "random"), ("LEACH", "hsa"), ("EE-LEACH*", "random"), ("EE-LEACH*", "hsa"),
          ("I-LEACH*", "random"), ("I-LEACH*", "hsa"), ("4F", "random"), ("ViTAL", "hsa")]
MOBILE = [("LEACH", "hsa"), ("EE-LEACH*", "hsa"), ("I-LEACH*", "hsa"), ("ViTAL", "hsa")]
W_VARIANTS = {"equal": (0.25, 0.25, 0.25, 0.25), "energy_only": (1, 0, 0, 0), "neigh_only": (0, 1, 0, 0),
              "dist_only": (0, 0, 1, 0), "snr_only": (0, 0, 0, 1), "energy_heavy": (0.4, 0.2, 0.2, 0.2)}
_dep = {}


def make_dep(a):
    kind, N, s = a
    xy = deploy(kind, N, s)
    dmin = cKDTree(xy).query(xy, k=2)[0][:, 1].min()
    return f"{kind}_{N}_{s}", xy, coverage(xy), dmin


def init():
    global _dep
    _dep = dict(np.load(DEP))


def job(a):
    proto, kind, N, s, kw = a
    o = simulate(proto, N, s, xy=_dep[f"{kind}_{N}_{s}"], **{k: v for k, v in kw.items() if k != "wname"})
    meta = {k: (list(v) if k == "w" else v) for k, v in kw.items()}
    return dict(protocol=proto, deploy=kind, N=N, seed=s) | meta | o


if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    with Pool(64) as pool:
        if not os.path.exists(DEP):
            deps = pool.map(make_dep, [(k, N, s) for k in ("random", "hsa") for N in NS for s in SEEDS])
            np.savez_compressed(DEP, **{n: xy for n, xy, _, _ in deps})
            json.dump([dict(name=n, cr=c, dmin=d) for n, _, c, d in deps], open(f"{OUT}/deployment_quality.json", "w"))
    with Pool(64, initializer=init) as pool:
        jobs = [(p, k, N, s, dict(L=L)) for L in (800, 4000) for p, k in STATIC for N in NS for s in SEEDS]
        jobs += [(p, k, N, s, dict(L=L, mobile=True)) for L in (800, 4000) for p, k in MOBILE for N in NS for s in SEEDS]
        jobs += [(p, k, N, s, dict(L=800, w=w, wname=name)) for p, k in (("4F", "random"), ("ViTAL", "hsa"))
                 for name, w in W_VARIANTS.items() for N in (100, 300, 500) for s in range(10)]
        jobs += [(p, k, N, s, dict(L=800, dist_inverted=True, wname="equal_inverted_dist"))
                 for p, k in (("4F", "random"), ("ViTAL", "hsa")) for N in (100, 300, 500) for s in range(10)]
        rows = pool.map(job, jobs, chunksize=1)
    json.dump(rows, open(f"{OUT}/runs.json", "w"))
    print("done", len(rows))
