"""Shared data loading and statistics helpers for the paper builders."""
import json
import os

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
R = json.load(open(f"{ROOT}/results/runs.json"))
Q = json.load(open(f"{ROOT}/results/robustness_sleep_cost.json"))
V1 = json.load(open(f"{ROOT}/../EBViTAL_Novelty_on_Ali2025/results/runs.json"))
NS, N3 = [100, 200, 300, 400, 500], [100, 300, 500]
COMP = ["LEACH (HSA)", "DEEC", "ViTAL", "EB-ViTAL"]
IDX = {}
for r in R:
    IDX.setdefault((r["scenario"], r["label"], r["N"]), []).append(r)


def get(sc, lab, N, k):
    S = sorted(IDX[(sc, lab, N)], key=lambda r: r["seed"])
    return np.array([(r["rounds"] if r[k] is None else r[k]) if k in ("T90", "T50") else
                     (np.nan if r[k] is None else r[k]) for r in S], float)


def mci(v):
    v = v[~np.isnan(v)]
    return v.mean(), 1.96 * v.std(ddof=1) / np.sqrt(len(v))


def pm(v, scale=1.0, dec=0):
    m, c = mci(v)
    return f"{m/scale:.{dec}f} ± {c/scale:.{dec}f}"


def gain(a, b):
    return 100 * (np.nanmean(a) / np.nanmean(b) - 1)


def rng_(xs):
    lo, hi = min(xs), max(xs)
    return f"{lo:+.0f}%" if round(lo) == round(hi) else f"{lo:+.0f}% to {hi:+.0f}%"


def best(sc, N):
    return max(COMP, key=lambda c: np.nanmean(get(sc, c, N, "T90")))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj, run = np.empty_like(p), 0.0
    for rank, i in enumerate(o):
        run = max(run, (len(p) - rank) * p[i])
        adj[i] = min(1.0, run)
    return adj


