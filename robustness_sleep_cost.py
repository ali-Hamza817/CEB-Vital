"""POST-HOC robustness check (not pre-registered): sleeping nodes pay sleep_frac * L * Eelec per round.
Always-on protocols have no sleeping nodes, so their results are unchanged and are read from runs.json."""
import json
from multiprocessing import Pool

import numpy as np

from ceb_sim import simulate
from run_experiments import SCEN


def init():
    global DEP
    DEP = dict(np.load("results/test_deployments.npz"))


def job(a):
    sc, label, proto, extra, N, s, f = a
    cfg = SCEN[sc]
    o = simulate(proto, N, s, xy=DEP[f"hsa_{cfg['side']}_{N}_{s}"], rmax=100000, sleep_frac=f, **cfg, **extra)
    o.pop("deaths")
    return dict(scenario=sc, label=label, N=N, seed=s, sleep_frac=f) | o


if __name__ == "__main__":
    P = [("CEB-ViTAL", "EB-ViTAL", {"sched": "cov"}), ("LEACH + CovSched", "LEACH", {"sched": "cov"})]
    jobs = [(sc, lab, p, ex, N, s, f) for f in (0.05, 0.25, 1.0) for sc in ("S1", "S3") for lab, p, ex in P
            for N in (500, 300, 100) for s in range(1000, 1030)]
    with Pool(64, initializer=init) as pool:
        rows = pool.map(job, jobs, chunksize=1)
    json.dump(rows, open("results/robustness_sleep_cost.json", "w"))
    print("done", len(rows))
