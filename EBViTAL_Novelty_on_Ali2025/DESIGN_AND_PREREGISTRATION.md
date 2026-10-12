# EB-ViTAL: design and pre-registration

Written 2026-10-12T14:55Z, before any EB-ViTAL simulation was run. The parameters, hypotheses and tests below are frozen. Anything added later is labelled post hoc.

**Base paper.** S. J. Ali, J. Iqbal, F. Khan, Y.-K. Moon, "Vigorous technique for augmented lifetime in WSNs", *Sci. Rep.* 15:15400 (2025). The paper is referred to below as ViTAL.

**Evaluation framework.** The simulator and assumptions P1–P18 are taken unchanged from the reconstruction in `../Ali2025_SciRep_ViTAL_WSN_Reproduction` (see its PREREGISTRATION.md). All protocols share one round loop, one radio model, the same deployments and the same random streams.

## 1. Gaps in ViTAL

| ID | Gap | Evidence |
|---|---|---|
| G1 | The energy term in Eq. 20 is Eres / Erestotal. Its typical value is about 1/N, while the other three terms are about 1. CH choice is therefore almost blind to energy. | At N = 100 with equal weights, the energy term contributes about 0.0025, against about 0.25 for each other term. |
| G2 | The distance term Dist / Avg_dist *raises* the score of nodes far from the sink. With the sink at the corner (5,5), those CHs pay the d⁴ amplifier cost. The paper's text says the opposite preference. | Eq. 20 as printed. |
| G3 | CH candidates are restricted to the CDS backbone, about 28 nodes for any N. Those nodes are drained first. | Reconstruction: ViTAL FND is 50–91% below LEACH, and it falls as N grows. |
| G4 | The number of CHs per round is not controlled. NV × T(n) shrinks the LEACH threshold (NV ≈ 0.5), and nothing ties it to the optimal cluster count. | Reconstruction: no-CH rounds are 40% for ViTAL vs 31% for LEACH (N = 100). |
| G5 | CHs send to the sink in one hop, although the sink sits in a corner. Distances reach about 127 m, beyond d0 = 87.7 m (the d⁴ regime), and well beyond the paper's stated 15 m radio range. | Field geometry. |
| G6 | The four weights are fixed. The paper names dynamic weighting as future work. | Paper, Conclusion. |

## 2. Proposed protocol: EB-ViTAL (Energy-Balanced ViTAL)

EB-ViTAL keeps ViTAL's phase 1 (HSA deployment) and its four factors: energy, neighbours, distance and SNR. It changes three things.

**M1. Normalised, direction-corrected score with a dynamic energy weight (fixes G1, G2, G6).** For each alive node i:

- e_i = Eres_i / max(Eres over alive nodes)
- n_i = Neigh_i / max(Neigh)
- d_i = 1 − Dist_i / max(Dist), so nodes closer to the sink score higher
- s_i = SNR_i / max(SNR)
- score_i = w1(t)·e_i + w(t)·(n_i + d_i + s_i)
- w1(t) = 0.25 + 0.5·(1 − Ē(t)/E0), where Ē(t) is the mean residual energy of alive nodes
- w(t) = (1 − w1(t)) / 3

The energy weight rises from 0.25 to 0.75 as the network drains.

**M2. Deterministic, spatially separated election of k_opt CHs (fixes G3, G4).**

- **Cluster count.** k_opt = sqrt(n_alive / 2π) · sqrt(εfs/εmp) · M / d̄², the classical optimum of Heinzelman et al. (2002), where d̄ is the mean node-to-sink distance of alive nodes. k_opt is rounded, kept at least 1 and at most n_alive.
- **Selection.** Alive nodes are taken in descending score order. A node is accepted as CH if it is at least r_sep = sqrt(A / (π·k_opt)) from every CH already accepted. Selection stops at k_opt CHs or when the list is exhausted.
- **No CDS restriction.** CH candidates are no longer limited to the CDS. The separation rule gives the well-spread coverage the CDS was meant to give.

**M3. Energy-aware multi-hop CH-to-sink routing (fixes G5).**

- **Routing tree.** A shortest-path tree to the sink is built over {CHs + sink}. Link cost is c(i→j) = E_tx(L, d_ij) + L·Eelec·(Ē / E_j) for a relay CH j, and E_tx(L, d_i,sink) for the sink.
- **Forwarding.** Relays forward packets without further aggregation. Each hop charges the sender E_tx and the receiver L·Eelec.

Members join the nearest CH, or the sink if it is nearer, exactly as in all other protocols.

## 3. Ablations

Each ablation removes exactly one component.

| Variant | M1 | M2 | M3 | Notes |
|---|---|---|---|---|
| EB-ViTAL | yes | yes | yes | full proposal |
| EB −M1 | ViTAL's Eq. 20 as printed, equal weights | yes | yes | |
| EB −M2 | yes | no | yes | M1 score used as a LEACH threshold: rand ≤ score·T(n), over all nodes |
| EB −M3 | yes | yes | no | CHs send to the sink in one hop |
| EB-ViTAL (random deployment) | yes | yes | yes | isolates the HSA phase |

## 4. Comparison set

All protocols use HSA deployment unless stated otherwise:

- LEACH (random deployment)
- LEACH (HSA)
- ViTAL, as reconstructed with the CDS restriction
- ViTAL-noCDS
- EB-ViTAL and its ablations

Packet sizes L = 800 (primary) and 4000 bits. Static and mobile-sink scenarios. N = 100–500, with 30 seeds; the same seed means the same deployment and sink trajectory for every protocol.

## 5. Hypotheses and tests (primary scenario: static, L = 800)

- **H1.** EB-ViTAL FND > ViTAL (reconstructed) FND at every N.
- **H2.** EB-ViTAL FND > LEACH (HSA) FND at every N.
- **H3.** EB-ViTAL HND > LEACH (HSA) HND at every N.

**Test.** One-sided paired Wilcoxon signed-rank over the 30 seeds, with Holm correction across the five N values, at family-wise α = 0.05 per hypothesis.

**Secondary, reported without hypothesis.** LND, sensor readings delivered (Σ alive nodes over rounds), packets at the sink, residual energy at round 1000, the L = 4000 runs, the mobile sink, and the ablation deltas.

**Fairness check.** Before any EB run, the copied simulator must reproduce `baseline_runs_from_reconstruction.json` bit-for-bit for every non-EB configuration.
