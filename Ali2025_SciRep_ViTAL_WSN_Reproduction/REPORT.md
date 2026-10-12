# Independent reconstruction of ViTAL (Ali, Iqbal, Khan & Moon, *Scientific Reports* 15:15400, 2025)

**Target paper.** S. J. Ali, J. Iqbal, F. Khan, Y.-K. Moon, "Vigorous technique for augmented lifetime in WSNs", *Scientific Reports* 15, 15400 (2025). https://doi.org/10.1038/s41598-025-96830-w. The paper is open access under CC BY 4.0. It was received on 16 Sep 2024 and accepted on 1 Apr 2025.

**Status of original artefacts.** No code or data are public. The paper says datasets are "available from the corresponding author on reasonable request". GitHub searches found nothing; the three hits for "ViTAL WSN" are unrelated projects.

**Bottom line.**
- **Method is underspecified.** The paper does not specify enough of ViTAL to reproduce it exactly. Eleven modelling choices had to be assumed. All were fixed in [`PREREGISTRATION.md`](PREREGISTRATION.md) before any ViTAL run.
- **Main result does not reproduce.** Under those assumptions ViTAL does **not** reproduce the paper's main result. Its first-node death is shorter than plain LEACH's, not 2.4–5× longer.
- **Two other gaps cannot be resolved here.** The paper's LEACH baseline matches a packet size of about 4000–8000 bits, not the stated 800. And in a scan of packet size and electronics energy, no setting fits the paper's LEACH numbers while still allowing ViTAL's reported lifetimes.
- **Headline percentages.** These are computed as (ViTAL − baseline) / ViTAL. In the usual convention, (ViTAL − baseline) / baseline, the figures imply much larger gains, for example +437% instead of 80.9%.

**Hardware.** The machine has 80 CPU cores, 4× NVIDIA L40S GPUs and 503 GB RAM, far more than needed. The 4,440 simulations plus 300 harmony-search deployments ran in about 3 minutes on 64 cores. No GPU is required.

---

## 1. What the paper specifies and what it leaves out

| Component | Given in paper | Missing, so assumed (pre-registered ID) |
|---|---|---|
| Field, sink, N, E0 | 100×100 m, sink (5,5), N = 100–500, 0.5 J | — |
| Radio constants | εfs 10 pJ, εmp 0.0013 pJ (labelled "energy receive"), 800-bit packets | Eelec, EDA, CH fraction p (P3, P5) |
| 4F score, Eq. 20 | four ratio terms | weights w1–w4 (P6); SNR model parameters (P8) |
| Eq. 22 and Algorithm 2 | ViTAL = NV × T(n); "select the candidate with the highest ViTAL value" | whether this is one CH per round or a LEACH-style threshold (P9) |
| CDS phase | "each node builds its CDS"; "CHs selected from the CDS"; "integrated after 10% of rounds" | construction algorithm and update rule (P10) |
| HSA deployment | Eq. 9–19, Algorithm 1 | HMS, HMCR, PAR, iterations, Ru (P11); how "optimal" deployment relates to fixed N |
| Baselines | one sentence each for EE-LEACH and I-LEACH | equations. The originals could not be retrieved, so EE-LEACH\* and I-LEACH\* are reconstructions (P13) |
| Metrics | FND, LND, packets to BS, average residual energy | how packets are counted (P15); when residual energy is measured (P16) |

## 2. Results

Each configuration has 30 runs per N. The paper's values were digitised from its Fig. 6 and Fig. 7 and are stored in [`paper_digitized/paper_values.json`](paper_digitized/paper_values.json), with an uncertainty of about ±30 rounds. Full tables are in [`results/summary.md`](results/summary.md).

### 2.1 The LEACH baseline decides which packet size is plausible

All values are for N = 100.

| Setting | LEACH (random) first node dead | LEACH last node dead |
|---|---|---|
| Paper Fig. 6 | 560 | 1180 |
| Reproduction, 800 bits (as stated) | 4844 ± 112 | 10174 ± 1483 |
| Reproduction, 4000 bits | 942 ± 27 | 2043 ± 291 |
| Reproduction, 8000 bits (scan, 5 runs) | 457 | 1026 |

With the stated 800-bit packets, every protocol lives 5–8× longer than the paper's baselines. The paper's LEACH numbers fit packets of roughly 4000–8000 bits.

### 2.2 ViTAL against the paper (static scenario, N = 100)

| Series | Paper FND | Repro FND, 800 bits | Repro FND, 4000 bits | Paper LND | Repro LND, 800 / 4000 bits |
|---|---|---|---|---|---|
| LEACH, random deployment | 560 | 4844 | 942 | 1180 | 10174 / 2043 |
| 4F selection, random deployment | 2800 | 4378 | 852 | 3660 | 9793 / 1971 |
| ViTAL (HSA + CDS + 4F) | 3000 | 2100 | 417 | 4050 | 11187 / 2220 |
| ViTAL without CDS restriction (post hoc) | — | 4404 | 853 | — | 8158 / 1660 |

The pattern holds for every N from 100 to 500, in both static and mobile-sink scenarios (figures R1–R4):
- **First node dies earlier.** ViTAL's first-node death is 50–91% *below* LEACH's.
- **Last node dies slightly later.** ViTAL's last-node death is 2–37% above LEACH's.
- **More packets counted.** ViTAL counts 69–252% more packets at the sink.
- **More residual energy only with 4000-bit packets.** There, residual energy at round 1000 is 44–67% above LEACH. With 800-bit packets the difference is within ±1%.

### 2.3 Headline improvement claims

The paper's percentages match (V − B)/V computed from its own figures to within about 1 point for first-node death. Examples: 81.3% against the stated 80.9%, 61.6% against 62%, and 58.0% against 58%. The usual convention, (V − B)/B, gives these values from the same figures:

| Static, vs randomly deployed baselines | Paper text | Same figure, (V−B)/B | Reproduction, (V−B)/B, 800 / 4000 bits |
|---|---|---|---|
| FND vs LEACH | 80.9% | +437% | −78% / −79% |
| LND vs LEACH | 75% | +293% | +6% / +6% |
| Packets vs LEACH | 83% | +494% | +78% / +75% |
| Residual energy vs LEACH | 68% | +148% | −1% / +48% |

Mobile-sink results follow the same pattern; see §4–5 of `results/summary.md`.

## 3. Why the reconstruction diverges

1. **The CDS restriction overloads a few nodes.** A connected dominating set at a 15 m radio range holds only about 28 nodes on average, whether N is 100 or 500. If cluster heads must come from it, those nodes serve repeatedly and die first. FND falls as N grows, which is the opposite of the paper's trend.
2. **4F scoring alone does not beat LEACH here.** As printed, Eq. 20's distance term w3·Dist/Avg_dist *raises* the score of nodes far from the sink. With the sink in a corner, those cluster heads pay the d⁴ amplifier cost.
   - Of the six weight settings tried, only distance-only (w3 = 1) brings 4F level with LEACH: 4795 vs 4844 at N = 100, and up to 4% above LEACH at N = 300–500. That is far from the paper's +400%.
   - Inverting the distance term does not help (§6 of the summary).
3. **Packet counts depend on the counting rule.** Under ViTAL, more rounds have no eligible cluster head: 40% vs 31% of the lifetime at N = 100. In those rounds every node sends directly to the sink, and each such packet is counted.
4. **No energy setting fits all reported numbers.** A 24-point scan covered packet size (800–32,000 bits) and Eelec (5–50 nJ/bit). No setting gave both the paper's LEACH FND and LND within ±25% while keeping the per-node lifetime cap E0/(L·Eelec) above ViTAL's reported 3000 rounds (§8 of the summary). A finer grid or undisclosed mechanisms could change this, for example duty cycling, suppression of transmissions through the CDS, or a different packet count.

## 4. Inconsistencies noticed in the paper

These are reported neutrally; each may be a typesetting or reporting slip.

1. **Coverage target.** With Rs = 5 m, 100 sensors can cover at most 78.5% of the 100×100 m field, so the 90% coverage target is infeasible at N = 100. Measured HSA coverage was 0.56.
2. **Fig. 5 field.** It shows a 50×50 field with roughly 10 m sensing disks, not the stated set-up.
3. **Run length.** The paper says "each simulation run lasts for 1000 seconds", yet it reports lifetimes up to 7600 rounds and gives no round duration.
4. **Static deployment.** The text says the static scenario uses grid-based deployment. Fig. 6 compares random and HSA deployment.
5. **Packet counts too low.** Fig. 6c gives LEACH about 950 packets to the BS at N = 100 over about 1180 rounds. That is under one packet per round, while LEACH with any usual p sends several cluster-head packets per round.
6. **Radio constants mislabelled.** Table 3 labels the amplifier constants "energy receive" and gives them in pJ/bit without the /m² and /m⁴ units. Eelec is absent.
7. **Distance term direction.** Eq. 20's distance term favours nodes far from the BS, while the text says better-positioned nodes are preferred.
8. **Single vs multiple CHs.** Algorithm 2 selects the single highest-scoring candidate per round. Eq. 22 has the form of a LEACH threshold.
9. **Percentage convention.** The improvement percentages use (V − B)/V (§2.3).
10. **Plot style.** The figures are in Python/matplotlib style, while the paper states MATLAB 2019a.

## 5. Limitations of this reconstruction

- **Choices outside the paper.** Eleven pre-registered choices (P3, P5, P6, P8–P11, P13–P16) are not in the paper. Different choices could change the absolute numbers. The robustness checks bound some of them: weights, distance term, CDS restriction, packet size and Eelec.
- **EE-LEACH\* and I-LEACH\*.** These are reconstructions from one-sentence descriptions. Comparisons with the paper's EE-LEACH and I-LEACH curves are therefore weak. LEACH itself is exact.
- **Packet-level effects.** The model is the first-order radio model, as in LEACH studies. The paper mentions log-normal shadowing and an IEEE 802.15.4 MAC; those packet-level effects are not modelled beyond SNR entering the 4F score.
- **Post-hoc checks.** §2.2's no-CDS row and the energy scan were added after seeing the main results. They are labelled as such and do not replace the pre-registered primary runs.
- **Digitised values.** Paper values were read from the figures, with an uncertainty of about ±30 rounds, ±150 packets and ±0.005 J.
- **IoT scenario not attempted.** The in-home health-monitoring scenario (Figs. 8–9) gives no parameters.

## 6. Files and re-running

```bash
python3 vital_repro.py          # self-checks: coverage, CDS domination/connectivity, run sanity
python3 run_experiments.py      # HSA deployments + 4,020 runs  -> results/runs.json  (~3 min, 64 cores)
python3 extra_experiments.py    # post-hoc checks             -> results/extra.json (~20 s)
python3 analyze.py              # tables + figures            -> results/summary.md, figures/R1..R4
```

Requires numpy, scipy and matplotlib. All randomness is seeded.

## 7. Use and attribution

- **Label.** Present this as an *independent reconstruction under stated assumptions*, citing the paper above.
- **Paper's numbers.** Quote them only as "reported by Ali et al.".
- **Gaps.** Describe the gaps as non-reproduced results, not as errors by the authors.
- **Next step.** Javed Iqbal is a co-author, and the paper offers data on request. Asking for the MATLAB code or the exact parameters (Eelec, p, w1–w4, SNR model, CDS rule, packet counting) is the step that would settle the open points.
