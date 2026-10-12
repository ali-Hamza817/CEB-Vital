"""Fills paper/PAPER_DRAFT_template.md with numbers computed from the run data -> paper/PAPER_DRAFT.md."""
import os

import numpy as np
from scipy.stats import wilcoxon

from common import COMP, HERE, IDX, N3, NS, Q, V1, best, gain, get, holm, mci, pm, rng_

F = {}
F["ABS_T90_BEST"] = rng_([gain(get("S1", "CEB-ViTAL", N, "T90"), get("S1", best("S1", N), N, "T90")) for N in NS])
F["ABS_FND_LEACH"] = rng_([gain(get("S1", "CEB-ViTAL", N, "FND"), get("S1", "LEACH (HSA)", N, "FND")) for N in NS])
F["ABS_FND_VITAL"] = rng_([gain(get("S1", "CEB-ViTAL", N, "FND"), get("S1", "ViTAL", N, "FND")) for N in NS])
F["ABS_COV_BEST"] = rng_([gain(get("S1", "CEB-ViTAL", N, "cov_int"), get("S1", best("S1", N), N, "cov_int")) for N in NS])
F["ABS_READ"] = rng_([gain(get("S1", "CEB-ViTAL", N, "readings"), get("S1", best("S1", N), N, "readings")) for N in NS])
F["ABS_READ_ALL"] = rng_([gain(get("S1", "CEB-ViTAL", N, "readings"), get("S1", c, N, "readings")) for N in NS for c in COMP])
lnd = lambda lab, N: np.nanmean(get("S1", lab, N, "LND"))
F["LND_NOTE"] = (f"Last-node death is not uniformly better. At N = 100, CEB-ViTAL's LND ({lnd('CEB-ViTAL', 100):.0f} rounds) is below "
                 f"LEACH (HSA) ({lnd('LEACH (HSA)', 100):.0f}) and ViTAL ({lnd('ViTAL', 100):.0f}). At N = 500 it exceeds every "
                 f"always-on baseline ({lnd('CEB-ViTAL', 500):.0f} vs at most {max(lnd(c, 500) for c in COMP):.0f}).")
F["ABS_S3_T90"] = rng_([gain(get("S3", "CEB-ViTAL", N, "T90"), get("S3", best("S3", N), N, "T90")) for N in NS])
F["C0_S3_N100"] = f"{100*np.mean([r['C0'] for r in IDX[('S3', 'CEB-ViTAL', 100)]]):.0f}%"
v1 = lambda lab, N: np.mean([r["FND"] for r in V1 if r["label"] == lab and r["N"] == N and r["L"] == 800 and not r["mobile"]])
drops = [100 * (1 - v1("EB -M1", N) / v1("EB-ViTAL", N)) for N in NS]
F["ABL_M1"] = f"{min(drops):.0f}–{max(drops):.0f}% (v1 study, development seeds)"

# hypotheses
hyp, pmax = [], 0.0
for hid, k in (("H1", "T90"), ("H2", "cov_int"), ("H3", "FND")):
    ps, gs = [], []
    for c in COMP:
        for N in NS:
            a, b = get("S1", "CEB-ViTAL", N, k), get("S1", c, N, k)
            ps.append(wilcoxon(a, b, alternative="greater").pvalue)
            gs.append((c, gain(a, b)))
    ph = holm(ps)
    pmax = max(pmax, ph.max())
    for c in COMP:
        sel = [i for i, (cc, _) in enumerate(gs) if cc == c]
        hyp.append(f"| {hid} ({k}) | {c} | {rng_([gs[i][1] for i in sel])} | {ph[sel].max():.1e} | {(ph[sel] < 0.05).sum()}/5 |")
F["ABS_PMAX"] = f"{pmax:.1e}"
F["TABLE_HYP"] = ("**Table 2. Pre-registered hypotheses in S1.** One-sided paired Wilcoxon over 30 held-out seeds, with Holm correction over 20 tests per hypothesis.\n\n"
                  "| Hypothesis | Comparator | CEB-ViTAL gain (N = 100–500) | Largest Holm p | Significant N |\n|---|---|---|---|---|\n" + "\n".join(hyp))

F["TABLE_SCEN_DEF"] = ("**Table S. Scenarios.** Every scenario uses E0 = 0.5 J, Rs = 5 m, N = 100–500 and 30 seeds per point.\n\n"
                       "| ID | Field | Sink | Packet |\n|---|---|---|---|\n"
                       "| S1 | 100 × 100 m | (5,5), static (base paper) | 800 bits |\n| S2 | 100 × 100 m | (50,50), static | 800 bits |\n"
                       "| S3 | 200 × 200 m | (5,5), static | 800 bits |\n| S4 | 200 × 200 m | (100,100), static | 800 bits |\n"
                       "| S5 | 100 × 100 m | random waypoint, 1–10 m/s, 10 s pause | 800 bits |\n| S6 | 100 × 100 m | (5,5), static | 4000 bits |")

LABS = ["LEACH (random)", "LEACH (HSA)", "DEEC", "ViTAL", "EB-ViTAL", "LEACH + CovSched", "CEB-ViTAL"]
t1a = ["**Table 1a. S1 lifetimes in rounds**, mean ± 95% CI over 30 held-out seeds.\n",
       "| Protocol | T90, N=100 | T90, N=300 | T90, N=500 | FND, N=100 | FND, N=300 | FND, N=500 |", "|---|---|---|---|---|---|---|"]
t1b = ["\n**Table 1b. S1 last node death (rounds), coverage integral (coverage-rounds) and readings delivered (×10³).**\n",
       "| Protocol | LND, N=100 | LND, N=500 | Cov. integral, N=100 | Cov. integral, N=500 | Readings, N=100 | Readings, N=500 |", "|---|---|---|---|---|---|---|"]
for lab in LABS:
    b = "**" if lab == "CEB-ViTAL" else ""
    t1a.append(f"| {b}{lab}{b} | " + " | ".join(pm(get("S1", lab, N, "T90")) for N in N3) + " | " +
               " | ".join(pm(get("S1", lab, N, "FND")) for N in N3) + " |")
    t1b.append(f"| {b}{lab}{b} | " + " | ".join(pm(get("S1", lab, N, "LND")) for N in (100, 500)) + " | " +
               " | ".join(pm(get("S1", lab, N, "cov_int")) for N in (100, 500)) + " | " +
               " | ".join(pm(get("S1", lab, N, "readings"), 1e3, 1) for N in (100, 500)) + " |")
F["TABLE_S1"] = "\n".join(t1a + t1b)

cross = ["**Table 3. CEB-ViTAL against the strongest always-on baseline** (highest mean T90 among LEACH (HSA), DEEC, ViTAL and EB-ViTAL).\n",
         "| Scenario | N | Strongest baseline | T90 | FND | Coverage integral | LND | Readings |", "|---|---|---|---|---|---|---|---|"]
for sc in ("S1", "S2", "S3", "S4", "S5", "S6"):
    for N in N3:
        bb = best(sc, N)
        g = [gain(get(sc, "CEB-ViTAL", N, k), get(sc, bb, N, k)) for k in ("T90", "FND", "cov_int", "LND", "readings")]
        cross.append(f"| {sc} | {N} | {bb} | " + " | ".join(f"{x:+.0f}%" for x in g) + " |")
F["TABLE_CROSS"] = "\n".join(cross)

abl = ["**Table 4. Ablation.** Change in CEB-ViTAL relative to each variant, shown as T90 / FND / coverage integral. "
       "'≈0†' means the variant's mean T90 is below 100 rounds, so a ratio is meaningless.\n",
       "| Scenario | Variant | N = 100 | N = 300 | N = 500 |", "|---|---|---|---|---|"]
VAR = [("EB-ViTAL", "no CovSched (always on)"), ("EB-ViTAL + RandSleep", "random sleepers, same count"),
       ("LEACH + CovSched", "LEACH clustering instead of EB"), ("CEB-ViTAL -M3", "no multi-hop (M3)")]
for sc in ("S1", "S3"):
    for lab, what in VAR:
        cells = []
        for N in N3:
            parts = []
            for k in ("T90", "FND", "cov_int"):
                b = get(sc, lab, N, k)
                parts.append("≈0†" if np.nanmean(b) < 100 else f"{gain(get(sc, 'CEB-ViTAL', N, k), b):+.0f}%")
            cells.append(" / ".join(parts))
        abl.append(f"| {sc} | {what} | " + " | ".join(cells) + " |")
F["TABLE_ABL"] = "\n".join(abl)

rob = ["**Table 5. Sleep-energy robustness (post hoc).** CEB-ViTAL against the strongest always-on baseline, shown as T90 / FND / coverage integral. "
       "Each sleeping node pays f · L · Eelec per round.\n",
       "| Scenario | N | f = 0 | f = 0.05 | f = 0.25 | f = 1 |", "|---|---|---|---|---|---|"]
for sc in ("S1", "S3"):
    for N in N3:
        bb = best(sc, N)
        base = [np.nanmean(get(sc, bb, N, k)) for k in ("T90", "FND", "cov_int")]
        cells = []
        for f in (0.0, 0.05, 0.25, 1.0):
            if f == 0:
                v = [np.nanmean(get(sc, "CEB-ViTAL", N, k)) for k in ("T90", "FND", "cov_int")]
            else:
                S = [r for r in Q if r["scenario"] == sc and r["label"] == "CEB-ViTAL" and r["N"] == N and r["sleep_frac"] == f]
                v = [np.mean([(r[k] if r[k] is not None else r["rounds"]) for r in S]) for k in ("T90", "FND", "cov_int")]
            cells.append(" / ".join(f"{100*(x/b-1):+.0f}%" for x, b in zip(v, base)))
        rob.append(f"| {sc} | {N} | " + " | ".join(cells) + " |")
F["TABLE_ROB"] = "\n".join(rob)

t = open(f"{HERE}/PAPER_DRAFT_template.md").read()
for k, v in F.items():
    t = t.replace("{{" + k + "}}", v)
assert "{{" not in t, [x for x in t.split("{{")[1:]][:3]
open(f"{HERE}/PAPER_DRAFT.md", "w").write(t)
print({k: v for k, v in F.items() if not k.startswith("TABLE")})
