"""Builds the repository README.md from paper/README_template.md and the run data (no hand-typed results)."""
import json
import os
import re
import sys

import numpy as np

from build_paper import F
from common import COMP, IDX, N3, NS, Q, ROOT, V1, best, gain, get, rng_

sys.path.insert(0, ROOT)
from ceb_sim import greedy_cds  # noqa: E402

facts = json.load(open(f"{ROOT}/results/figure_facts.json"))
X = dict(F)
X["VITAL_VS_LEACH"] = rng_([gain(get("S1", "ViTAL", N, "FND"), get("S1", "LEACH (HSA)", N, "FND")) for N in NS])
dep = dict(np.load(f"{ROOT}/results/test_deployments.npz"))
X["CDS_SIZE"] = f"{np.mean([greedy_cds(dep[f'hsa_100_{N}_{s}'], np.ones(N, bool)).sum() for N in (100, 300, 500) for s in range(1000, 1010)]):.0f}"

# ---- EB-ViTAL v1 (development seeds)
v1 = lambda lab, N, k: np.nanmean([np.nan if r[k] is None else r[k] for r in V1 if r["label"] == lab and r["N"] == N
                                   and r["L"] == 800 and not r["mobile"]])
X["V1_FND"] = rng_([100 * (v1("EB-ViTAL", N, "FND") / v1("LEACH (HSA)", N, "FND") - 1) for N in NS])
X["V1_LND"] = rng_([100 * (v1("EB-ViTAL", N, "LND") / v1("LEACH (HSA)", N, "LND") - 1) for N in NS])
X["V1_READ"] = rng_([100 * (v1("EB-ViTAL", N, "readings") / v1("LEACH (HSA)", N, "readings") - 1) for N in NS])

# ---- CEB-ViTAL extras
X["NRUNS"] = f"{sum(len(v) for v in IDX.values()):,}"
X["NROB"] = f"{len(Q):,}"
X["C0_S1_100"] = f"{100*np.mean([r['C0'] for r in IDX[('S1', 'CEB-ViTAL', 100)]]):.0f}%"
X["C0_S1_500"] = f"{100*np.mean([r['C0'] for r in IDX[('S1', 'CEB-ViTAL', 500)]]):.0f}%"
X["LEACHCOV_S1"] = rng_([gain(get("S1", "CEB-ViTAL", N, "T90"), get("S1", "LEACH + CovSched", N, "T90")) for N in NS])
X["LEACHCOV_S3"] = rng_([gain(get("S3", "CEB-ViTAL", N, "T90"), get("S3", "LEACH + CovSched", N, "T90")) for N in NS])
X["SCHED_S1"] = rng_([gain(get("S1", "CEB-ViTAL", N, "T90"), get("S1", "EB-ViTAL", N, "T90")) for N in NS])
X["SCHED_S3"] = rng_([gain(get("S3", "CEB-ViTAL", N, "T90"), get("S3", "EB-ViTAL", N, "T90")) for N in NS])
X["M3_S3"] = rng_([gain(get("S3", "CEB-ViTAL", N, "T90"), get("S3", "CEB-ViTAL -M3", N, "T90")) for N in NS])


def rob(sc, f):
    out = []
    for N in N3:
        b = np.nanmean(get(sc, best(sc, N), N, "T90"))
        S = [r for r in Q if r["scenario"] == sc and r["label"] == "CEB-ViTAL" and r["N"] == N and r["sleep_frac"] == f]
        out.append(100 * (np.mean([(r["T90"] if r["T90"] is not None else r["rounds"]) for r in S]) / b - 1))
    return out


X["ROB_S1_F1"], X["ROB_S3_F1"] = rng_(rob("S1", 1.0)), rng_(rob("S3", 1.0))
X["HEAT_MIN"], X["HEAT_MAX"] = f"+{facts['heat_min']:.0f}%", f"+{facts['heat_max']:.0f}%"
for key, lab in (("LEACH", "LEACH (HSA)"), ("VITAL", "ViTAL"), ("EB", "EB-ViTAL"), ("CEB", "CEB-ViTAL")):
    d = facts[lab]
    X[f"SNAP_{key}"] = f"{d['alive']} alive, {d['sensing']} sensing, {d['chs']} CHs, {100*d['cov']:.0f}% coverage, mean residual energy {d['meanE']:.2f} J"
X["SENS_CEB"] = f"{facts['sensing_mid']['CEB-ViTAL']:.0f}"
X["T90_CEB300"] = f"{np.nanmean(get('S1', 'CEB-ViTAL', 300, 'T90')):.0f}"
X["T90_BEST300"] = f"{np.nanmean(get('S1', best('S1', 300), 300, 'T90')):.0f}"

m, e = F["ABS_PMAX"].split("e")
X["PMAX_PRETTY"] = f"{m} × 10<sup>{int(e)}</sup>"
X = {k: re.sub(r"(?<![\w.])-(?=\d)", "−", v).replace("−0%", "0%") for k, v in X.items()}   # typographic minus
t = open(os.path.join(os.path.dirname(__file__), "README_template.md")).read()
for k, v in X.items():
    t = t.replace("{{" + k + "}}", v)
left = re.findall(r"\{\{([A-Z0-9_]+)\}\}", t)
assert not left, left
open(f"{ROOT}/README.md", "w").write(t)
print({k: v for k, v in X.items() if not k.startswith("TABLE") and k != "LND_NOTE"})
