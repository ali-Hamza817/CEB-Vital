# CEB-ViTAL: Coverage-Aware, Energy-Balanced Clustering for Lifetime Extension of Wireless Sensor Networks

*Draft manuscript. Every number in this file is filled in by `build_paper.py` from `results/runs.json` and `results/robustness_sleep_cost.json`; edit the template, not the output.*

## Abstract

ViTAL (Ali et al., *Sci. Rep.* 2025) extends wireless-sensor-network (WSN) lifetime in three phases: harmony-search deployment for coverage, a connected dominating set (CDS), and a four-factor cluster-head (CH) score. We identify three weaknesses in its CH selection:

- its energy term is normalised by total network energy, so it has almost no effect;
- its distance term favours CHs far from the sink;
- restricting CHs to the CDS overloads a small backbone.

ViTAL's operating phases also ignore the coverage that its deployment phase optimises.

We propose **CEB-ViTAL**. It keeps ViTAL's deployment phase and adds:

- a coverage-preserving, energy-weighted sleep scheduler (CovSched);
- a normalised four-factor score with a dynamic energy weight;
- deterministic election of the analytically optimal number of spatially separated CHs;
- energy-aware multi-hop CH-to-sink forwarding.

All protocols were simulated in one shared simulator: 6,900 runs on 30 held-out random networks, 6 scenarios and 100–500 nodes. All design constants and hypotheses were pre-registered before the test runs.

In the base-paper setting (100 m field, corner sink), CEB-ViTAL extends the coverage lifetime T90 by +21% to +130% over the strongest baseline. It extends first-node death by +40% to +184% over LEACH and +224% to +2519% over ViTAL, and the coverage integral by +15% to +148%. All pre-registered tests are significant at a Holm-corrected p ≤ 6.1e-06. Because sleeping nodes do not report, total readings delivered differ from the always-on baselines by -4% to +1%.

Gains are largest in dense networks and smallest in a sparse 200 m field (+2% to +17% T90). They persist when sleeping nodes are charged a non-zero sleep energy, except in the sparse 200 m field under the most pessimistic setting.

**Keywords:** wireless sensor networks, clustering, LEACH, coverage, sleep scheduling, energy balance, network lifetime

## 1. Introduction

Clustered WSN protocols descended from LEACH [2] rotate the cluster-head role to spread energy use. ViTAL [1] combines three ideas:

- coverage-driven deployment by harmony search (HSA);
- a connected-dominating-set (CDS) backbone;
- a four-factor CH score built from residual energy, neighbour count, distance to the sink and SNR.

ViTAL's operating phases, however, keep every node sensing and reporting every round, even where neighbours already cover its area. Its own conclusion names scheduling as future work.

This paper makes four contributions:

1. **Analysis of ViTAL's CH selection.** We show that it is almost energy-blind (§3.2).
2. **CEB-ViTAL.** The protocol couples coverage-preserving sleep scheduling with corrected, energy-balanced clustering (§4).
3. **Pre-registered evaluation.** Hypotheses were fixed before testing, on held-out random networks, against LEACH [2], DEEC [4], ViTAL [1] and our earlier extension EB-ViTAL, with ablations (§5–6).
4. **Open simulator and data.** Every number in this paper can be regenerated with one command.

## 2. Related work (to be completed by the authors)

The literature for this section still has to be read and cited; the points below only outline it.

- **LEACH and its centralised variant.** LEACH [2] introduced randomised CH rotation, and LEACH-C [3] derived the optimal number of clusters used in our election step.
- **Energy-aware CH probabilities.** DEEC [4] weights CH probability by residual energy relative to the network average.
- **Competition radius.** HEED [5] uses a competition radius to spread CHs.
- **Coverage-preserving scheduling.** Greedy weighted set cover [6] is the standard approximation for selecting a covering subset.
- **Still to survey.** Coverage-preserving clustering and duty-cycling protocols (2023–2026) must be searched. The novelty claim of §4 is limited to the integration with ViTAL and the corrections in §3.2.

## 3. System model and base protocol

### 3.1 Network and energy model

- **Field and sink.** N sensors are deployed in an M × M field, with a sink at a fixed position or following a random-waypoint path.
- **Energy.** Each node starts with E0 = 0.5 J. The first-order radio model [2] uses Eelec = 50 nJ/bit, εfs = 10 pJ/bit/m², εmp = 0.0013 pJ/bit/m⁴, d0 = 87.7 m and aggregation cost EDA = 5 nJ/bit.
- **Traffic.** Each sensing node produces one L-bit reading per round.
- **Coverage.** A node covers 1 m demand cells within Rs = 5 m (ViTAL's binary disk model). Coverage C(r) is the covered fraction of the field in round r, and C0 = C(0) with all nodes active.

### 3.2 ViTAL and its weaknesses

ViTAL's score (Eq. 20 of [1]) is

NV = w1·Eres/Erestotal + w2·Neigh/N + w3·Dist/Avg_dist + w4·SNR/SNR_avg

and a node is elected if rand ≤ NV·T(n).

- **G1, near energy-blindness.** Erestotal is the network total, so the energy term is about 1/N while the others are about 1. In our ablation, a deterministic election using this score drains the same nodes repeatedly: first-node death falls by 93–97% (v1 study, development seeds).
- **G2, wrong distance direction.** Dist/Avg_dist increases with distance, favouring far CHs, which pay the d⁴ amplifier cost when the sink is in a corner.
- **G3, CDS bottleneck.** Restricting CHs to the CDS (about 28 nodes for any N) concentrates load. ViTAL's first-node death is the lowest of all protocols and falls as N grows (Table 1).
- **G4, no coverage awareness.** Coverage is optimised at deployment and then ignored.

## 4. CEB-ViTAL

**Phase 1, deployment.** Unchanged from ViTAL: HSA maximising coverage × minimum spacing.

**Phase 2, CovSched.** Every r_s = 10 rounds, and whenever an active node dies, the sink computes the active set:

- Let U be the cells coverable by all alive nodes.
- Greedily add the node maximising (new cells covered) × (E_i / E_max) until τ·|U| cells are covered, with τ = 0.95.
- The remaining nodes sleep: they neither sense nor transmit.
- Every alive node pays the energy to receive a 200-bit schedule at each (re)schedule.

**Phase 3, clustering among active nodes.**

- *(M1) Score:* s_i = w1(t)·E_i/max E + w(t)·[Neigh_i/max Neigh + (1 − Dist_i/max Dist) + SNR_i/max SNR], with w1(t) = 0.25 + 0.5·(1 − Ē/E0) and w(t) = (1 − w1)/3.
- *(M2) Election:* k = round( sqrt(n/2π) · sqrt(εfs/εmp) · M / d̄² ) [3]. Nodes are taken in descending score; a node is accepted if it lies at least sqrt(M²/(πk)) from every accepted CH.
- *(M3) Routing:* a shortest-path tree to the sink over CHs. Link cost is E_tx(L, d) plus L·Eelec·Ē/E_j for a relay CH j.

**Complexity per schedule.** The greedy cover costs O(|A|·nnz), where nnz ≈ 78·N incidences and |A| is the active-set size. Per round, the election costs O(N log N + kN) and routing costs O(k²) (Dijkstra on k + 1 nodes). Like LEACH-C and ViTAL, the method is centralised at the sink.

## 5. Evaluation methodology

- **Pre-registration.** Constants, scenarios, hypotheses and tests were frozen before the test runs (`PREREGISTRATION_v2.md`).
- **Held-out test set.** Results use 30 held-out seeds (1000–1029) that were never used in development. A seed fixes the deployment and the random streams for every protocol.
- **Fairness check.** Before testing, the simulator reproduced all 7,620 earlier runs bit-for-bit.
- **Scenarios.** Table S lists S1–S6.
- **Metrics.** The primary metric is T90, the first round in which C(r) < 0.9·C0. Secondary metrics are T50, the coverage integral Σ_r C(r), first/half/last node death (FND/HND/LND) and readings delivered.
- **Tests.** One-sided paired Wilcoxon signed-rank, with Holm correction over 20 tests per hypothesis.
- **Compute.** All runs took under one hour on a 64-core server.

**Table S. Scenarios.** Every scenario uses E0 = 0.5 J, Rs = 5 m, N = 100–500 and 30 seeds per point.

| ID | Field | Sink | Packet |
|---|---|---|---|
| S1 | 100 × 100 m | (5,5), static (base paper) | 800 bits |
| S2 | 100 × 100 m | (50,50), static | 800 bits |
| S3 | 200 × 200 m | (5,5), static | 800 bits |
| S4 | 200 × 200 m | (100,100), static | 800 bits |
| S5 | 100 × 100 m | random waypoint, 1–10 m/s, 10 s pause | 800 bits |
| S6 | 100 × 100 m | (5,5), static | 4000 bits |

## 6. Results

### 6.1 Base-paper setting (S1)

**Table 1a. S1 lifetimes in rounds**, mean ± 95% CI over 30 held-out seeds.

| Protocol | T90, N=100 | T90, N=300 | T90, N=500 | FND, N=100 | FND, N=300 | FND, N=500 |
|---|---|---|---|---|---|---|
| LEACH (random) | 5170 ± 27 | 5386 ± 18 | 5505 ± 16 | 4825 ± 46 | 4665 ± 46 | 4549 ± 53 |
| LEACH (HSA) | 5230 ± 21 | 5508 ± 18 | 5652 ± 13 | 4931 ± 39 | 4750 ± 42 | 4650 ± 37 |
| DEEC | 5502 ± 18 | 5639 ± 11 | 5674 ± 10 | 5378 ± 22 | 5330 ± 24 | 5257 ± 17 |
| ViTAL | 3048 ± 60 | 4529 ± 86 | 5193 ± 61 | 2127 ± 38 | 818 ± 20 | 504 ± 17 |
| EB-ViTAL | 5716 ± 17 | 5758 ± 10 | 5767 ± 6 | 5676 ± 37 | 5667 ± 33 | 5650 ± 24 |
| LEACH + CovSched | 5632 ± 22 | 6645 ± 54 | 8411 ± 148 | 5621 ± 22 | 6621 ± 55 | 8324 ± 148 |
| **CEB-ViTAL** | 6900 ± 58 | 10686 ± 55 | 13290 ± 181 | 6889 ± 59 | 10642 ± 57 | 13195 ± 180 |

**Table 1b. S1 last node death (rounds), coverage integral (coverage-rounds) and readings delivered (×10³).**

| Protocol | LND, N=100 | LND, N=500 | Cov. integral, N=100 | Cov. integral, N=500 | Readings, N=100 | Readings, N=500 |
|---|---|---|---|---|---|---|
| LEACH (random) | 10087 ± 543 | 11187 ± 324 | 3075 ± 37 | 6025 ± 11 | 586.1 ± 2.8 | 2944.6 ± 2.9 |
| LEACH (HSA) | 8954 ± 330 | 9272 ± 368 | 3346 ± 35 | 5814 ± 21 | 584.8 ± 1.3 | 2944.5 ± 2.4 |
| DEEC | 7447 ± 375 | 7949 ± 448 | 3287 ± 35 | 5560 ± 20 | 580.3 ± 1.1 | 2948.5 ± 2.2 |
| ViTAL | 11376 ± 289 | 11696 ± 216 | 3375 ± 37 | 6695 ± 37 | 565.5 ± 3.5 | 2816.3 ± 10.5 |
| EB-ViTAL | 6039 ± 368 | 6760 ± 459 | 3240 ± 38 | 5426 ± 25 | 575.0 ± 1.5 | 2919.3 ± 2.8 |
| LEACH + CovSched | 9226 ± 324 | 18222 ± 442 | 3502 ± 27 | 10744 ± 69 | 573.8 ± 1.3 | 2858.7 ± 3.1 |
| **CEB-ViTAL** | 7294 ± 309 | 16415 ± 193 | 3723 ± 24 | 13440 ± 62 | 566.8 ± 1.7 | 2820.0 ± 2.9 |

**Table 2. Pre-registered hypotheses in S1.** One-sided paired Wilcoxon over 30 held-out seeds, with Holm correction over 20 tests per hypothesis.

| Hypothesis | Comparator | CEB-ViTAL gain (N = 100–500) | Largest Holm p | Significant N |
|---|---|---|---|---|
| H1 (T90) | LEACH (HSA) | +32% to +135% | 5.2e-06 | 5/5 |
| H1 (T90) | DEEC | +25% to +134% | 5.2e-06 | 5/5 |
| H1 (T90) | ViTAL | +126% to +156% | 5.2e-06 | 5/5 |
| H1 (T90) | EB-ViTAL | +21% to +130% | 5.2e-06 | 5/5 |
| H2 (cov_int) | LEACH (HSA) | +11% to +131% | 1.9e-08 | 5/5 |
| H2 (cov_int) | DEEC | +13% to +142% | 1.9e-08 | 5/5 |
| H2 (cov_int) | ViTAL | +10% to +101% | 1.9e-08 | 5/5 |
| H2 (cov_int) | EB-ViTAL | +15% to +148% | 1.9e-08 | 5/5 |
| H3 (FND) | LEACH (HSA) | +40% to +184% | 6.1e-06 | 5/5 |
| H3 (FND) | DEEC | +28% to +151% | 6.1e-06 | 5/5 |
| H3 (FND) | ViTAL | +224% to +2519% | 6.1e-06 | 5/5 |
| H3 (FND) | EB-ViTAL | +21% to +134% | 6.1e-06 | 5/5 |

Last-node death is not uniformly better. At N = 100, CEB-ViTAL's LND (7294 rounds) is below LEACH (HSA) (8954) and ViTAL (11376). At N = 500 it exceeds every always-on baseline (16415 vs at most 11696).

Figure F1 plots the S1 metrics against N, and Figure F9 shows coverage over time at N = 300. CEB-ViTAL holds coverage at about 95% of the always-on level, as designed, and keeps it roughly twice as long. Always-on protocols lose coverage once nodes start to die.

### 6.2 Other scenarios

**Table 3. CEB-ViTAL against the strongest always-on baseline** (highest mean T90 among LEACH (HSA), DEEC, ViTAL and EB-ViTAL).

| Scenario | N | Strongest baseline | T90 | FND | Coverage integral | LND | Readings |
|---|---|---|---|---|---|---|---|
| S1 | 100 | EB-ViTAL | +21% | +21% | +15% | +21% | -1% |
| S1 | 300 | EB-ViTAL | +86% | +88% | +80% | +89% | -2% |
| S1 | 500 | EB-ViTAL | +130% | +134% | +148% | +143% | -3% |
| S2 | 100 | EB-ViTAL | +14% | +17% | +11% | +7% | -1% |
| S2 | 300 | ViTAL | +67% | +696% | +27% | +24% | -15% |
| S2 | 500 | ViTAL | +90% | +1561% | +68% | +75% | -16% |
| S3 | 100 | EB-ViTAL | +12% | +141% | -0% | +1% | -1% |
| S3 | 300 | EB-ViTAL | +11% | +184% | +9% | +11% | -1% |
| S3 | 500 | EB-ViTAL | +17% | +251% | +20% | +18% | -1% |
| S4 | 100 | EB-ViTAL | +4% | +7% | -0% | +5% | -2% |
| S4 | 300 | EB-ViTAL | +16% | +18% | +11% | +19% | -3% |
| S4 | 500 | EB-ViTAL | +29% | +30% | +24% | +37% | -3% |
| S5 | 100 | EB-ViTAL | +20% | +20% | +15% | +24% | -1% |
| S5 | 300 | DEEC | +82% | +89% | +76% | +93% | -2% |
| S5 | 500 | ViTAL | +114% | +2036% | +82% | +66% | -9% |
| S6 | 100 | EB-ViTAL | +22% | +22% | +16% | +22% | -0% |
| S6 | 300 | EB-ViTAL | +89% | +91% | +82% | +92% | -1% |
| S6 | 500 | EB-ViTAL | +140% | +141% | +148% | +147% | -1% |

- **Sparse fields.** Gains shrink when there is little redundant coverage to exploit: N = 100 in a 200 m field covers only about 19% of the area.
- **Dense fields.** Gains grow with node density.
- **Mobile sink and larger packets.** The S5 and S6 results show the same pattern as S1.

### 6.3 Ablation

**Table 4. Ablation.** Change in CEB-ViTAL relative to each variant, shown as T90 / FND / coverage integral. '≈0†' means the variant's mean T90 is below 100 rounds, so a ratio is meaningless.

| Scenario | Variant | N = 100 | N = 300 | N = 500 |
|---|---|---|---|---|
| S1 | no CovSched (always on) | +21% / +21% / +15% | +86% / +88% / +80% | +130% / +134% / +148% |
| S1 | random sleepers, same count | ≈0† / -1% / +10% | ≈0† / -9% / +17% | ≈0† / -26% / +12% |
| S1 | LEACH clustering instead of EB | +23% / +23% / +6% | +61% / +61% / +23% | +58% / +59% / +25% |
| S1 | no multi-hop (M3) | +0% / +0% / -0% | +0% / +0% / +0% | +0% / +0% / +0% |
| S3 | no CovSched (always on) | +12% / +141% / -0% | +11% / +184% / +9% | +17% / +251% / +20% |
| S3 | random sleepers, same count | +35% / +129% / +1% | ≈0† / +137% / +7% | ≈0† / +146% / +10% |
| S3 | LEACH clustering instead of EB | +41% / +41% / +8% | +110% / +108% / +24% | +140% / +138% / +36% |
| S3 | no multi-hop (M3) | +3% / +3% / +3% | -0% / +1% / +2% | +2% / +4% / +3% |

- **CovSched** gives the largest gain.
- **Coverage-aware selection is essential.** Sleeping the same number of nodes chosen at random drops coverage below 0.9·C0 within the first rounds in every setting except the sparsest ones (N = 100 in the 200 m fields).
- **EB clustering beats LEACH clustering** under the same scheduler.
- **Multi-hop forwarding (M3)** matters only in the 200 m field with a corner sink. With a centre sink no relay is ever cheaper than a direct link, so the M3 variant gives identical results.

### 6.4 Robustness to sleep energy (post hoc)

The model above assumes that sleeping is free. This check was not pre-registered. Each sleeping node is charged f·L·Eelec per round, which penalises only the scheduled protocols.

**Table 5. Sleep-energy robustness (post hoc).** CEB-ViTAL against the strongest always-on baseline, shown as T90 / FND / coverage integral. Each sleeping node pays f · L · Eelec per round.

| Scenario | N | f = 0 | f = 0.05 | f = 0.25 | f = 1 |
|---|---|---|---|---|---|
| S1 | 100 | +21% / +21% / +15% | +20% / +21% / +14% | +18% / +18% / +12% | +10% / +10% / +4% |
| S1 | 300 | +86% / +88% / +80% | +83% / +86% / +78% | +75% / +75% / +66% | +37% / +35% / +28% |
| S1 | 500 | +130% / +134% / +148% | +125% / +128% / +141% | +108% / +105% / +115% | +56% / +49% / +45% |
| S3 | 100 | +12% / +141% / -0% | +11% / +140% / -0% | +10% / +137% / -0% | +6% / +129% / -1% |
| S3 | 300 | +11% / +184% / +9% | +10% / +182% / +9% | +6% / +174% / +7% | -3% / +144% / +3% |
| S3 | 500 | +17% / +251% / +20% | +15% / +248% / +20% | +10% / +234% / +17% | +0% / +157% / +9% |

## 7. Discussion and threats to validity

- **Coverage level.** CEB-ViTAL trades a fixed ~5% of instantaneous coverage (τ = 0.95) for a much longer coverage lifetime. Applications needing every possible covered cell should use a larger τ.
- **Readings delivered.** These differ by -4% to +1% from the always-on baselines in S1, because sleeping nodes do not report. Redundant readings from overlapping sensors are what is dropped.
- **Model limits.** The model is the first-order radio model with no packet loss, MAC contention or idle listening. The sleep-energy check bounds one of these effects.
- **Baselines.** ViTAL is our reconstruction from its paper. Some parameters are unpublished (Eelec, p, w1–w4, SNR model), and its absolute numbers differ from those reported in [1]. LEACH and DEEC are standard and fully specified.
- **Centralised control.** The scheduler and election are centralised, as in LEACH-C and ViTAL. A distributed version is future work.
- **Untuned constants.** τ, r_s, L_CTRL and the M1 weight schedule were fixed a priori and not tuned. Better values may exist.

## 8. Conclusion

Coupling ViTAL's coverage-driven deployment with coverage-preserving scheduling and corrected, energy-balanced clustering extends coverage lifetime substantially in dense networks, with no meaningful loss of delivered data. The pre-registered hypotheses hold in all 20 comparisons in the base-paper setting. Future work:

- distributed scheduling;
- validation against the original ViTAL implementation;
- packet-level simulation;
- hardware tests.

## Data and code availability

The simulator, pre-registration, raw run data and analysis scripts are in this folder. Run `python3 run_experiments.py && python3 robustness_sleep_cost.py && python3 analyze.py && python3 paper/build_paper.py` to regenerate every table and figure.

## References (verify details before submission)

1. S. J. Ali, J. Iqbal, F. Khan, Y.-K. Moon, "Vigorous technique for augmented lifetime in WSNs," *Scientific Reports*, 15, 15400, 2025. doi:10.1038/s41598-025-96830-w
2. W. R. Heinzelman, A. Chandrakasan, H. Balakrishnan, "Energy-efficient communication protocol for wireless microsensor networks," *Proc. HICSS-33*, 2000.
3. W. B. Heinzelman, A. P. Chandrakasan, H. Balakrishnan, "An application-specific protocol architecture for wireless microsensor networks," *IEEE Trans. Wireless Commun.*, 1(4), 660–670, 2002.
4. L. Qing, Q. Zhu, M. Wang, "Design of a distributed energy-efficient clustering algorithm for heterogeneous wireless sensor networks," *Computer Communications*, 29(12), 2230–2237, 2006.
5. O. Younis, S. Fahmy, "HEED: a hybrid, energy-efficient, distributed clustering approach for ad hoc sensor networks," *IEEE Trans. Mobile Comput.*, 3(4), 366–379, 2004.
6. V. Chvátal, "A greedy heuristic for the set-covering problem," *Mathematics of Operations Research*, 4(3), 233–235, 1979.
