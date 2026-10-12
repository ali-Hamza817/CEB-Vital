# EB-ViTAL: an energy-balanced extension of ViTAL for clustered WSNs

**Base paper.** S. J. Ali, J. Iqbal, F. Khan, Y.-K. Moon, "Vigorous technique for augmented lifetime in WSNs", *Scientific Reports* 15:15400 (2025). It is referred to below as ViTAL.

## Summary

ViTAL has three phases: harmony-search (HSA) deployment, a connected dominating set (CDS), and a four-factor (4F) cluster-head score. We re-implemented it from the paper and found three weaknesses in its CH selection:

- **Energy-blind scoring.** The energy term is divided by the *total* network energy, so it is about 1/N and has almost no effect.
- **Wrong distance direction.** The distance term favours nodes *far* from the sink.
- **CDS bottleneck.** Restricting CHs to the CDS concentrates the load on about 28 nodes.

**EB-ViTAL** keeps ViTAL's deployment phase and its four factors. It adds:

- **M1.** A normalised, direction-corrected score whose energy weight rises as the network drains.
- **M2.** Deterministic election of the analytically optimal number of CHs, kept apart by a minimum separation.
- **M3.** Energy-aware multi-hop forwarding from CHs to the sink.

**Results.** The evaluation used 30 seeds per point, N = 100–500 and the sink in a corner. All pre-registered hypotheses were tested with a one-sided paired Wilcoxon test and Holm correction.

| Comparison | First node dead (stability period) | Half nodes dead | Last node dead | Readings delivered |
|---|---|---|---|---|
| vs ViTAL (reconstructed) | **+170% to +1059%** | ±0% | −41% to −46% | +2% to +4% |
| vs LEACH, same HSA deployment | **+15% to +22%** | −0.5% to +1.4% | −26% to −31% | −1% to −2% |
| vs LEACH, random deployment | +17% to +24% | −0% to +2% | −37% to −42% | −1% to −2% |

- **H1 and H2 hold at every N.** EB-ViTAL's FND beats ViTAL and LEACH, with Holm-corrected p < 10⁻⁵.
- **H3 holds only at N = 100–200.** At larger N, half-node death is level with LEACH.
- **The trade-off is real.** All nodes stay alive much longer, then die within a few hundred rounds of each other. Last-node death and total delivered readings are therefore lower.
- **Same direction elsewhere.** The pattern holds for 4000-bit packets (FND +17% to +21% vs LEACH) and for a mobile sink (FND +14% to +16%).

## 1. Evaluation framework

The baseline simulator is the independent reconstruction of ViTAL in `../Ali2025_SciRep_ViTAL_WSN_Reproduction`. That reconstruction did **not** reproduce the paper's absolute numbers; see its REPORT.md. All comparisons here are therefore *within one simulator*, never against numbers printed in the paper.

**Shared settings.** Every protocol uses:

- the same first-order radio model: Eelec 50 nJ/bit, εfs 10 pJ/bit/m², εmp 0.0013 pJ/bit/m⁴, EDA 5 nJ/bit;
- a 100×100 m field with the sink at (5,5), E0 = 0.5 J and an 800-bit packet (primary), with 4000 bits as a check;
- the same deployments and random streams, so seed s gives the same node positions and sink trajectory for every protocol.

**Fairness check.** After EB-ViTAL was added, the simulator reproduced all 4,020 earlier reconstruction runs bit-for-bit (`results/fairness.json`).

## 2. Gaps in ViTAL

| ID | Gap | Evidence |
|---|---|---|
| G1 | Eq. 20 divides residual energy by the network total. The energy term is about 1/N, while the other terms are about 1, so CH choice ignores energy in practice. | Without the M1 correction, EB's deterministic election makes FND collapse to 156–424 rounds, because the same nodes are chosen again and again (§4.3). |
| G2 | The term Dist / Avg_dist rewards nodes far from the sink, against the paper's stated intent. With a corner sink, those CHs pay the d⁴ cost. | Eq. 20 as printed. |
| G3 | CH candidates are limited to the CDS, about 28 nodes for any N. | ViTAL's FND is 50–91% below LEACH and falls as N grows. |
| G4 | The number of CHs is uncontrolled. NV < 1 shrinks the LEACH threshold. | 1.4–1.8 CHs per round on average, and 40% of rounds have no CH (N = 100). |
| G5 | CHs reach the sink in one hop over up to about 127 m, which is the d⁴ regime. | Field geometry. |
| G6 | The weights are static. The paper lists dynamic weights as future work. | Conclusion of the base paper. |

## 3. EB-ViTAL

For each alive node i in round t, with Ē(t) the mean residual energy of alive nodes:

**M1. Score.**

score_i = w1(t)·Eres_i / max Eres + w(t)·[ Neigh_i / max Neigh + (1 − Dist_i / max Dist) + SNR_i / max SNR ]

- w1(t) = 0.25 + 0.5·(1 − Ē(t)/E0)
- w(t) = (1 − w1(t)) / 3

Each term lies in [0, 1]. Nodes closer to the sink score higher, and the energy weight grows from 0.25 to 0.75 as the network drains.

**M2. Election.**

- **Cluster count.** k = round( sqrt(n_alive/2π) · sqrt(εfs/εmp) · M / d̄² ), the optimum of Heinzelman et al. (2002), kept between 1 and n_alive.
- **Selection.** Nodes are scanned in descending score order. A node becomes CH if it lies at least r_sep = sqrt(M² / (π·k)) from every CH already chosen.

**M3. Routing.**

- **Tree.** Build a shortest-path tree to the sink over {CHs + sink}. Link cost is E_tx(L, d_ij) + L·Eelec·Ē/E_j for a relay CH j, and E_tx(L, d_i,sink) for the sink.
- **Forwarding.** Relays forward packets without re-aggregation.

Members join the nearest CH, or the sink if it is closer.

```
Algorithm EB-ViTAL (one round, computed at the sink as in LEACH-C)
1  deploy nodes once with ViTAL's HSA phase
2  for each round t:
3      compute Ē, max-normalised factors, w1(t)           # M1
4      score ← w1·e + w·(n + d + s)
5      k ← k_opt(n_alive, mean distance to sink)           # M2
6      CH ← ∅
7      for node i in descending score:
8          if |CH| = k: break
9          if dist(i, c) ≥ r_sep for all c in CH: add i to CH
10     members join nearest CH (or sink if nearer); CHs aggregate
11     build min-cost tree over CH ∪ {sink}; forward along tree   # M3
```

**Cost per round.** Sorting is O(N log N), the separation test is O(k·N) and Dijkstra over k + 1 nodes is O(k²). ViTAL already needs network-wide quantities (total energy, mean distance, mean SNR), so EB-ViTAL uses the same centralised information model.

## 4. Results

The full tables, p-values and every variant are in [`results/summary.md`](results/summary.md).

### 4.1 Pre-registered hypotheses (static, 800 bits)

| Hypothesis | Result |
|---|---|
| H1: EB-ViTAL FND > ViTAL | **Supported at all 5 N.** Mean 5633–5680 vs 486–2100 rounds; Holm p = 4.3×10⁻⁶. |
| H2: EB-ViTAL FND > LEACH (HSA) | **Supported at all 5 N.** +14.7% to +22.0%; Holm p ≤ 3.5×10⁻⁶. |
| H3: EB-ViTAL HND > LEACH (HSA) | **Partly supported.** Holds at N = 100 and 200 (+1.4%, +0.4%). Not at N = 300–500, where the change is −0.0% to −0.5%. |

### 4.2 Figures

- `figures/F1_static_L800.png`: FND, HND, LND and readings per node vs N, mean ± 95% CI.
- `figures/F5_alive_nodes_N100.png`: alive nodes over time. EB-ViTAL keeps all 100 nodes alive until about round 5600, while LEACH starts losing nodes about 600 rounds earlier and ViTAL about 3500 rounds earlier. EB-ViTAL's nodes then die close together.
- `figures/F2_static_L4000.png`, `figures/F3_mobile_L800.png`: the same pattern under 4000-bit packets and a mobile sink.
- `figures/F4_ablation_L800.png`: ablation (§4.3).

### 4.3 Ablation (static, 800 bits)

Each row shows how the full EB-ViTAL compares with that variant.

| Variant | Full EB-ViTAL FND vs variant | Full EB-ViTAL LND vs variant | Reading |
|---|---|---|---|
| −M1 (Eq. 20 as printed) | +1238% to +3512% | −42% to −45% | M1 is essential. Without proper energy normalisation, a deterministic election drains the same nodes. |
| −M2 (probabilistic threshold instead) | +10% to +14% | −21% to −25% | M2 adds about 10–14% to the stability period. |
| −M3 (one-hop CH to sink) | **−1% to −2%** | +0% to +6% | **M3 gives no stability gain in a 100 m field.** It trades 1–2% of FND for up to 6% more LND. |
| HSA vs random deployment | +2% to +4% | −19% to −27% | The HSA phase adds a small FND gain. |

The M3 result is reported as found. Multi-hop forwarding may pay off more in larger fields, where more CHs fall beyond d0; that was not tested.

### 4.4 Mobile sink and packet size

| Setting | EB-ViTAL vs LEACH (HSA): FND / HND / LND | EB-ViTAL vs ViTAL: FND / HND / LND |
|---|---|---|
| Mobile, 800 bits | +14–16% / −1–+1% / −16–−21% | +130–830% / −10–−14% / −39–−42% |
| Mobile, 4000 bits | +16–18% / −1–0% / −15–−18% | +134–882% / −10–−14% / −37–−39% |
| Static, 4000 bits | +17–21% / 0–+1% / −26–−31% | see summary §5 |

With a mobile sink, the reconstructed ViTAL keeps half its nodes alive longer than EB-ViTAL. EB-ViTAL's advantage is specifically the stability period.

## 5. Discussion

- **Where EB-ViTAL helps.** It suits applications that need *every* sensor working, such as coverage-critical monitoring, healthcare and intrusion detection. It holds full coverage 15–22% longer than LEACH and up to 11.6× longer than the ViTAL reconstruction.
- **What it gives up.** The long tail of a few surviving nodes, and 1–2% of total readings. Applications that value data volume over full coverage may prefer LEACH.
- **Natural follow-up (not tested).** Switch from M2 to probabilistic rotation once the first node dies, to recover part of the tail.

## 6. Novelty and prior art: check before claiming

The building blocks are established ideas:

- energy-normalised multi-criteria CH scores (e.g. DEEC-style residual and average energy);
- the k_opt cluster count (Heinzelman et al. 2002, LEACH-C);
- CH separation radii (HEED and EEUC-style competition ranges);
- multi-hop CH routing (EEM-LEACH, which ViTAL itself cites as ref. 29).

The contribution here is the **diagnosis of ViTAL's normalisation, distance-direction and CDS-restriction problems**, together with an integrated fix whose components are each tested by ablation. Before presenting this as novel, run a literature search on "dynamic weight cluster head selection", "k_opt deterministic cluster head", "LEACH-C separation radius" and "ViTAL WSN extension" for 2023–2026.

## 7. Limitations

- **Reconstruction, not original code.** Results come from the reconstruction framework, not the authors' code. Absolute numbers depend on the assumptions in the reconstruction's PREREGISTRATION.md. The relative comparison is fair because every protocol shares the simulator.
- **No radio-range limit.** All protocols may transmit over any distance. The paper's 15 m radio range is not enforced for data, as in LEACH studies.
- **Untested settings.** One field size (100 m) and one sink position for static runs were tested. Packet loss, MAC and control-message costs are not modelled.
- **Design fixed in advance.** EB-ViTAL's constants (0.25/0.5 weight schedule, r_sep, relay cost) were chosen before running and not tuned. Other choices might do better or worse.

## 8. Re-running

```bash
python3 wsn_sim.py          # self-checks, including M2 and M3 unit tests
python3 run_experiments.py  # fairness check (4,020 runs) + 3,600 EB-ViTAL runs, ~6 min on 64 cores
python3 analyze.py          # tests, tables, figures
```

The files are DESIGN_AND_PREREGISTRATION.md (frozen before runs), wsn_sim.py, run_experiments.py, analyze.py, results/ and figures/.

## Appendix: draft message to the base-paper author

> Dear Dr. Iqbal,
>
> I have been studying your 2025 Scientific Reports paper on ViTAL and built an independent simulation of it from the paper. While doing so I noticed that the energy term in Eq. 20 (Eres/Erestotal) becomes very small relative to the other three terms. I also noticed that restricting cluster heads to the CDS concentrates load on a small backbone.
>
> I designed an extension, EB-ViTAL, with three parts: a normalised score with a dynamic energy weight, election of the optimal number of well-separated cluster heads, and multi-hop routing. In my simulator it extends first-node lifetime by 15–22% over LEACH and well beyond my ViTAL reconstruction, with a trade-off in last-node lifetime.
>
> Because some parameters were not stated in the paper (Eelec, p, w1–w4, the SNR model), my baseline numbers differ from yours. Would you be willing to share your simulation code or those parameters, so I can validate EB-ViTAL against the original ViTAL implementation? I would be glad to work on this under your guidance.
