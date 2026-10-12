# CEB-ViTAL study: pre-registration

Frozen 2026-10-12 before any test-set run.

**Before freezing, I had seen:**
- the EB-ViTAL v1 results on development seeds 0–29 (`../EBViTAL_Novelty_on_Ali2025`);
- the `ceb_sim.py` self-check output on development seed 0, N = 100, which is debugging only.

No constant below was changed after seeing any run. All reported results use **held-out test seeds 1000–1029** with freshly generated deployments.

## 1. Motivation (from the v1 study)

- **Strength.** EB-ViTAL v1 extended the stability period, with FND +15–22% vs LEACH.
- **Weakness 1.** It shortened the network tail: LND −26–31%.
- **Weakness 2.** It delivered 1–2% fewer readings.
- **Weakness 3.** Every alive node still senses and transmits every round, even where its area is already covered by neighbours.

ViTAL's phase 1 is built around coverage (HSA, coverage ratio cr, threshold ct), but its operating phases ignore coverage. The base paper also lists scheduling as future work.

## 2. Proposed method: CEB-ViTAL (Coverage-aware, Energy-Balanced ViTAL)

| Phase | ViTAL | CEB-ViTAL |
|---|---|---|
| 1 Deployment | HSA, Eq. 9–19 | unchanged |
| 2 Topology/activity | CDS; CHs restricted to it | **CovSched**: every r_sched rounds, and whenever an active node dies, the sink computes an active set by energy-weighted greedy set cover. The rest sleep. |
| 3 Clustering | 4F score, Eq. 20–22 | EB-ViTAL among active nodes: M1 score, M2 k_opt with separation, M3 multi-hop |

**CovSched in detail.**
- **Target.** Let U be the cells coverable by all alive nodes, using ViTAL's binary disk model with Rs = 5 m and 1 m cells. Greedily add the node that maximises (new cells covered) × (E_i / E_max) until τ·|U| cells are covered.
- **Sleep.** Sleeping nodes neither sense nor transmit.
- **Control cost.** Every alive node pays L_CTRL·Eelec per (re)schedule to receive the schedule.

| Constant | Value |
|---|---|
| τ (coverage target, fraction of currently coverable area) | 0.95 |
| r_sched (rounds between schedules) | 10 |
| L_CTRL (schedule message, bits) | 200 |
| M1/M2/M3 constants | unchanged from DESIGN_AND_PREREGISTRATION.md (v1) |

## 3. Scenarios

All scenarios use E0 = 0.5 J, the ViTAL radio model, N = 100, 200, 300, 400, 500, and seeds 1000–1029. A seed gives the same deployment for all protocols within a field size and deployment type.

| ID | Field | Sink | Packet | Notes |
|---|---|---|---|---|
| S1 | 100×100 m | (5,5) static | 800 bits | **primary**; the base paper's setting |
| S2 | 100×100 m | (50,50) static | 800 | centre sink |
| S3 | 200×200 m | (5,5) static | 800 | large field, corner sink |
| S4 | 200×200 m | (100,100) static | 800 | large field, centre sink |
| S5 | 100×100 m | mobile, random waypoint from (5,5) | 800 | the paper's mobile scenario |
| S6 | 100×100 m | (5,5) static | 4000 | packet-size check |

## 4. Protocols

All use HSA deployment except LEACH (random).

| Label | Definition | Run in |
|---|---|---|
| LEACH (random) | LEACH, random deployment | S1–S4 |
| LEACH (HSA) | LEACH | S1–S6 |
| DEEC | DEEC, homogeneous (Qing et al. 2006) | S1–S6 |
| ViTAL | reconstruction, with CDS | S1–S6 |
| EB-ViTAL | v1, every node always active | S1–S6 |
| LEACH + CovSched | isolates the clustering contribution | S1–S4 |
| EB-ViTAL + RandSleep | same number of active nodes as CovSched, chosen at random; isolates node selection | S1–S4 |
| **CEB-ViTAL** | proposed | S1–S6 |
| CEB-ViTAL −M3 | single-hop CHs; tests whether multi-hop matters in the 200 m fields | S1–S4 |

## 5. Metrics

- **T90 (primary).** First round in which the area covered by sensing nodes falls below 90% of the initial full-deployment coverage C0. If coverage is still above that level when the run ends, the run is recorded as censored.
- **T50.** The same, with 50% of C0.
- **Coverage integral.** Σ_r coverage(r), counted in coverage-rounds.
- **FND, HND, LND.** First, half and last node death.
- **Readings delivered.** Σ sensing nodes per round.
- **Packets at the sink, control energy.**
- **Run length.** rmax = 100000 rounds.

## 6. Hypotheses (S1; one-sided paired Wilcoxon over the 30 test seeds)

Each hypothesis is tested against four comparators: LEACH (HSA), DEEC, ViTAL and EB-ViTAL. Holm correction runs over 5 N × 4 comparators = 20 tests per hypothesis, at family-wise α = 0.05.

- **H1.** CEB-ViTAL T90 > comparator.
- **H2.** CEB-ViTAL coverage integral > comparator.
- **H3.** CEB-ViTAL FND > comparator.

## 7. Analyses reported without hypotheses

- **Scenarios S2–S6.** All metrics, mean ± 95% CI.
- **Scheduling effect.** CEB-ViTAL vs EB-ViTAL.
- **Node-selection effect.** CEB-ViTAL vs EB-ViTAL + RandSleep.
- **Clustering effect.** CEB-ViTAL vs LEACH + CovSched.
- **Multi-hop effect.** CEB-ViTAL vs CEB-ViTAL −M3, especially S3.
- **Readings trade-off.** Readings delivered and LND are reported for every protocol, whatever the outcome.

## 8. Integrity rules

- **Fairness gate.** Before the test runs, `ceb_sim.py` must reproduce all 7,620 earlier runs bit-for-bit: 4,020 from the reconstruction and 3,600 from v1.
- **No re-tuning.** If a hypothesis fails, it is reported as failed. Any new variant is labelled post hoc and run on new seeds.
