"""Fairness gate + test-set deployments + all CEB-ViTAL runs (PREREGISTRATION_v2.md)."""
import json
import os
import sys
from multiprocessing import Pool

import numpy as np

from ceb_sim import deploy, simulate

NS = (100, 200, 300, 400, 500)
SEEDS = range(1000, 1030)
KEYS = ("FND", "HND", "LND", "packets", "res500", "res1000")
EBFLAGS = {"EB -M1": {"m1": False}, "EB -M2": {"m2": False}, "EB -M3": {"m3": False}}
SCEN = {"S1": dict(side=100, bs0=(5, 5), L=800), "S2": dict(side=100, bs0=(50, 50), L=800),
        "S3": dict(side=200, bs0=(5, 5), L=800), "S4": dict(side=200, bs0=(100, 100), L=800),
        "S5": dict(side=100, bs0=(5, 5), L=800, mobile=True), "S6": dict(side=100, bs0=(5, 5), L=4000)}
FULL = [("LEACH (random)", "LEACH", "random", {}), ("LEACH (HSA)", "LEACH", "hsa", {}),
        ("DEEC", "DEEC", "hsa", {}), ("ViTAL", "ViTAL", "hsa", {}), ("EB-ViTAL", "EB-ViTAL", "hsa", {}),
        ("LEACH + CovSched", "LEACH", "hsa", {"sched": "cov"}),
        ("EB-ViTAL + RandSleep", "EB-ViTAL", "hsa", {"sched": "rand"}),
        ("CEB-ViTAL", "EB-ViTAL", "hsa", {"sched": "cov"}),
        ("CEB-ViTAL -M3", "EB-ViTAL", "hsa", {"sched": "cov", "eb": {"m3": False}})]
CORE = ["LEACH (HSA)", "DEEC", "ViTAL", "EB-ViTAL", "CEB-ViTAL"]


def init(path):
    global DEP
    DEP = dict(np.load(path))


def regress(r):
    kw = {k: r[k] for k in ("L", "mobile", "dist_inverted") if k in r and r[k] is not None}
    if "w" in r:
        kw["w"] = tuple(r["w"])
    if "label" in r and r["label"] in EBFLAGS:
        kw["eb"] = EBFLAGS[r["label"]]
    o = simulate(r["protocol"], r["N"], r["seed"], xy=DEP[f"{r['deploy']}_{r['N']}_{r['seed']}"], **kw)
    return all(o[k] == r[k] for k in KEYS)


def make_dep(a):
    kind, side, N, s = a
    return f"{kind}_{side}_{N}_{s}", deploy(kind, N, s, side)


def job(a):
    sc, label, proto, kind, N, s, extra = a
    cfg = SCEN[sc]
    o = simulate(proto, N, s, xy=DEP[f"{kind}_{cfg['side']}_{N}_{s}"], rmax=100000, **cfg, **extra)
    return dict(scenario=sc, label=label, protocol=proto, deploy=kind, N=N, seed=s) | o


if __name__ == "__main__":
    os.makedirs("results", exist_ok=True)
    if not os.path.exists("results/fairness.json"):
        old = json.load(open("results/regress_reconstruction.json")) + json.load(open("results/regress_ebvital_v1.json"))
        with Pool(64, initializer=init, initargs=("results/dev_deployments.npz",)) as pool:
            same = pool.map(regress, old, chunksize=4)
        json.dump(dict(runs=len(old), identical=int(sum(same))), open("results/fairness.json", "w"))
        print("fairness:", sum(same), "/", len(old), flush=True)
        if not all(same):
            sys.exit("fairness gate failed")
    if not os.path.exists("results/test_deployments.npz"):
        todo = [(k, side, N, s) for side in (100, 200) for k in ("hsa", "random") for N in NS[::-1] for s in SEEDS]
        with Pool(64) as pool:
            deps = pool.map(make_dep, todo, chunksize=1)
        np.savez_compressed("results/test_deployments.npz", **dict(deps))
        print("deployments:", len(deps), flush=True)
    jobs = []
    for sc in SCEN:
        for label, proto, kind, extra in FULL:
            if sc in ("S5", "S6") and label not in CORE:
                continue
            jobs += [(sc, label, proto, kind, N, s, extra) for N in NS for s in SEEDS]
    jobs.sort(key=lambda j: (-j[4], "sched" not in j[6]))          # largest first
    with Pool(64, initializer=init, initargs=("results/test_deployments.npz",)) as pool:
        rows = []
        for i, row in enumerate(pool.imap_unordered(job, jobs, chunksize=1)):
            rows.append(row)
            if (i + 1) % 500 == 0:
                print(f"{i + 1}/{len(jobs)}", flush=True)
    json.dump(rows, open("results/runs.json", "w"))
    print("done", len(rows), flush=True)
