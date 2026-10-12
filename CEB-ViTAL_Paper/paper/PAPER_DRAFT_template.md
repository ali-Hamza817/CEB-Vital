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

In the base-paper setting (100 m field, corner sink), CEB-ViTAL extends the coverage lifetime T90 by {{ABS_T90_BEST}} over the strongest baseline. It extends first-node death by {{ABS_FND_LEACH}} over LEACH and {{ABS_FND_VITAL}} over ViTAL, and the coverage integral by {{ABS_COV_BEST}}. All pre-registered tests are significant at a Holm-corrected p ≤ {{ABS_PMAX}}. Because sleeping nodes do not report, total readings delivered differ from the always-on baselines by {{ABS_READ_ALL}}.

Gains are largest in dense networks and smallest in a sparse 200 m field ({{ABS_S3_T90}} T90). They persist when sleeping nodes are charged a non-zero sleep energy, except in the sparse 200 m field under the most pessimistic setting.

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

- **G1, near energy-blindness.** Erestotal is the network total, so the energy term is about 1/N while the others are about 1. In our ablation, a deterministic election using this score drains the same nodes repeatedly: first-node death falls by {{ABL_M1}}.
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

{{TABLE_SCEN_DEF}}

## 6. Results

### 6.1 Base-paper setting (S1)

{{TABLE_S1}}

{{TABLE_HYP}}

{{LND_NOTE}}

Figure F1 plots the S1 metrics against N, and Figure F9 shows coverage over time at N = 300. CEB-ViTAL holds coverage at about 95% of the always-on level, as designed, and keeps it roughly twice as long. Always-on protocols lose coverage once nodes start to die.

### 6.2 Other scenarios

{{TABLE_CROSS}}

- **Sparse fields.** Gains shrink when there is little redundant coverage to exploit: N = 100 in a 200 m field covers only about {{C0_S3_N100}} of the area.
- **Dense fields.** Gains grow with node density.
- **Mobile sink and larger packets.** The S5 and S6 results show the same pattern as S1.

### 6.3 Ablation

{{TABLE_ABL}}

- **CovSched** gives the largest gain.
- **Coverage-aware selection is essential.** Sleeping the same number of nodes chosen at random drops coverage below 0.9·C0 within the first rounds in every setting except the sparsest ones (N = 100 in the 200 m fields).
- **EB clustering beats LEACH clustering** under the same scheduler.
- **Multi-hop forwarding (M3)** matters only in the 200 m field with a corner sink. With a centre sink no relay is ever cheaper than a direct link, so the M3 variant gives identical results.

### 6.4 Robustness to sleep energy (post hoc)

The model above assumes that sleeping is free. This check was not pre-registered. Each sleeping node is charged f·L·Eelec per round, which penalises only the scheduled protocols.

{{TABLE_ROB}}

## 7. Discussion and threats to validity

- **Coverage level.** CEB-ViTAL trades a fixed ~5% of instantaneous coverage (τ = 0.95) for a much longer coverage lifetime. Applications needing every possible covered cell should use a larger τ.
- **Readings delivered.** These differ by {{ABS_READ_ALL}} from the always-on baselines in S1, because sleeping nodes do not report. Redundant readings from overlapping sensors are what is dropped.
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
