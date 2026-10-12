# CEB-ViTAL

**Coverage-aware, energy-balanced clustering for lifetime extension in wireless sensor networks.**

This repository extends ViTAL ([Ali, Iqbal, Khan & Moon, *Scientific Reports* 15:15400, 2025](https://doi.org/10.1038/s41598-025-96830-w)) in three steps:

- an independent re-implementation of ViTAL;
- an analysis of its cluster-head selection;
- a new protocol, **CEB-ViTAL**, which adds coverage-preserving sleep scheduling and corrected, energy-balanced clustering on top of ViTAL's coverage-driven deployment.

The repository holds the simulator, all raw run data, the analysis scripts and the LaTeX manuscript.

![Network state at round 5500](CEB-ViTAL_Paper/figures/fig_snapshot.png)

*Network state at round 5500, from one held-out network with 300 nodes. CEB-ViTAL keeps every node alive while only about half of them sense at a time.*

## Key results

These results are for the base setting of ViTAL: a 100 m field with the sink at (5, 5). They come from 30 held-out random networks per point and 100–500 nodes, and the hypotheses were pre-registered.

| CEB-ViTAL versus | Change |
|---|---|
| Strongest baseline, coverage lifetime T90 | +21% to +130% |
| LEACH, first node dead | +40% to +184% |
| ViTAL (re-implementation), first node dead | +224% to +2519% |
| Strongest baseline, coverage integral | +15% to +148% |
| Always-on baselines, readings delivered | −4% to +1% |

- **Significance.** All pre-registered tests are significant, with Holm-corrected p ≤ 6.1 × 10⁻⁶.
- **Scale of the study.** The study comprises 6,900 test runs over 6 scenarios.
- **Range of gains.** Across all scenarios the T90 gain ranges from +2% (a sparse 200 m field) to +140%.

![T90 gain across scenarios](CEB-ViTAL_Paper/figures/fig_gain_heatmap.png)

![S1 metrics](CEB-ViTAL_Paper/figures/fig_s1_metrics.png)

![Network dynamics](CEB-ViTAL_Paper/figures/fig_dynamics.png)

## Repository layout

| Folder | Contents |
|---|---|
| `Ali2025_SciRep_ViTAL_WSN_Reproduction/` | Independent re-implementation of ViTAL from its paper: assumptions (`PREREGISTRATION.md`), simulator, results and report |
| `EBViTAL_Novelty_on_Ali2025/` | EB-ViTAL v1: gap analysis, normalised score, k-opt separated cluster heads, multi-hop routing, and ablations |
| `CEB-ViTAL_Paper/` | CEB-ViTAL: simulator (`ceb_sim.py`), pre-registration, experiments, analysis, figures and the LaTeX manuscript (`latex/`) |

## Reproducing the results

You need Python 3 with numpy, scipy and matplotlib.

```bash
cd CEB-ViTAL_Paper
python3 run_experiments.py          # fairness gate + 6,900 test runs (about 1 h on 64 cores)
python3 robustness_sleep_cost.py    # post-hoc sleep-energy check
python3 analyze.py                  # statistics, tables, figures
cd paper && python3 make_latex_assets.py && python3 make_paper_figures.py && python3 make_extra_figures.py
cd ../latex && pdflatex main && bibtex main && pdflatex main && pdflatex main
```

All randomness is seeded. The test networks use seeds 1000–1029 and were never used during development.

## Notes

- **ViTAL baseline.** ViTAL is re-implemented from its paper. Several of its parameters are unpublished, so the absolute ViTAL numbers differ from those printed in the original paper.
- **Model.** The study uses the first-order radio model, without packet loss or MAC effects.

## Author

Ali Hamza
