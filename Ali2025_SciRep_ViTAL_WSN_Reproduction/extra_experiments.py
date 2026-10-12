"""Robustness checks added after the main runs (both are labelled post hoc in the report).

(a) ViTAL-noCDS: 4F selection on the HSA deployment, with CH candidates NOT restricted
    to the CDS (alternative reading of P10).
(b) Feasibility scan over packet size L and Eelec (Eelec is not given in the paper):
    LEACH (random deployment, N=100) FND/LND versus the per-node lifetime bound
    E0/(L*Eelec) that caps ViTAL if every alive node sends one packet per round.
"""
import json
from multiprocessing import Pool

import numpy as np

import vital_repro as V

NS = (100, 200, 300, 400, 500)


def init():
    global DEP
    DEP = dict(np.load("results/deployments.npz"))


def run_nocds(a):
    N, s, L = a
    o = V.simulate("4F", N, s, xy=DEP[f"hsa_{N}_{s}"], L=L)
    return dict(protocol="ViTAL-noCDS", deploy="hsa", N=N, seed=s, L=L) | o


def run_scan(a):
    L, eelec, s = a
    V.EELEC = eelec                       # module global read by tx() and CH receive cost
    o = V.simulate("LEACH", 100, s, xy=DEP[f"random_100_{s}"], L=L, rmax=200000)
    return dict(L=L, Eelec=eelec, seed=s, FND=o["FND"], LND=o["LND"], bound=V.E0 / (L * eelec))


if __name__ == "__main__":
    with Pool(64, initializer=init) as pool:
        nocds = pool.map(run_nocds, [(N, s, L) for L in (800, 4000) for N in NS for s in range(30)], chunksize=1)
        scan = pool.map(run_scan, [(L, e, s) for L in (800, 2000, 4000, 8000, 16000, 32000)
                                   for e in (5e-9, 10e-9, 25e-9, 50e-9) for s in range(5)], chunksize=1)
    json.dump(dict(nocds=nocds, scan=scan), open("results/extra.json", "w"))
    print("done", len(nocds), len(scan))
