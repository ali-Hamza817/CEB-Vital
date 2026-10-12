# CEB-ViTAL study

This study continues ViTAL (Ali et al., Sci. Rep. 15:15400, 2025).

Read the files in this order:

1. `PREREGISTRATION_v2.md`: design, hypotheses and tests, frozen before the test runs.
2. `paper/PAPER_DRAFT.md`: manuscript draft. Every number in it is generated from the data by `paper/build_paper.py`.
3. `results/summary.md`: every table, scenario and variant.
4. `figures/`: F1–F9.
5. `results/tables/`: CSV and LaTeX tables for the manuscript.

Code files:

- `ceb_sim.py`: the simulator. It covers LEACH, DEEC, ViTAL (reconstruction), EB-ViTAL and CEB-ViTAL, plus CovSched and RandSleep.
- `run_experiments.py`: the fairness gate (7,620 earlier runs reproduced bit-for-bit) and the 6,900 test runs on seeds 1000–1029.
- `robustness_sleep_cost.py`: a post-hoc sleep-energy check, 1,080 runs.
- `analyze.py`: the statistics, tables and figures.

To regenerate everything:

```bash
python3 run_experiments.py && python3 robustness_sleep_cost.py && python3 analyze.py && python3 paper/build_paper.py
```

The run takes about one hour on 64 cores.
