"""Journal-sized vector figures for latex/ (IEEE: column 3.5 in, text width 7.16 in)."""
import os
import sys
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import NS, ROOT, get, mci

sys.path.insert(0, ROOT)
from ceb_sim import simulate  # noqa: E402

OUT = os.path.join(ROOT, "latex", "figures")
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True, "grid.color": "#e4e3df", "grid.linewidth": 0.5,
                     "pdf.fonttype": 42})
STY = {"LEACH (HSA)": ("#2a78d6", "-", "o"), "DEEC": ("#6250d6", "-", "s"), "ViTAL": ("#eb6834", "-", "^"),
       "EB-ViTAL": ("#1baf7a", "--", "v"), "CEB-ViTAL": ("#1baf7a", "-", "D")}
KEYS = [("Coverage lifetime $T_{90}$ (rounds)", "T90", 1), ("First node dead (rounds)", "FND", 1),
        ("Coverage integral (coverage-rounds)", "cov_int", 1), ("Readings delivered ($\\times10^6$)", "readings", 1e6)]


def metrics(sc, fname):
    fig, axs = plt.subplots(2, 2, figsize=(7.16, 4.4))
    for ax, (title, k, sc_) in zip(axs.ravel(), KEYS):
        for lab, (c, ls, mk) in STY.items():
            m = [mci(get(sc, lab, N, k)) for N in NS]
            ax.errorbar(NS, [a / sc_ for a, _ in m], yerr=[b / sc_ for _, b in m], color=c, ls=ls, marker=mk,
                        ms=3.5, lw=1.4, capsize=2, label=lab)
        ax.set_title(title)
        ax.set_xticks(NS)
    for ax in axs[1]:
        ax.set_xlabel("Number of nodes $N$")
    fig.legend(*axs[0, 0].get_legend_handles_labels(), loc="upper center", ncol=5, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig(f"{OUT}/{fname}")
    fig.savefig(f"{ROOT}/figures/{fname[:-4]}.png", dpi=220)
    plt.close(fig)


DEP = None
CFG = {"LEACH (HSA)": ("LEACH", {}), "DEEC": ("DEEC", {}), "ViTAL": ("ViTAL", {}),
       "EB-ViTAL": ("EB-ViTAL", {}), "CEB-ViTAL": ("EB-ViTAL", {"sched": "cov"})}


def init():
    global DEP
    DEP = dict(np.load(f"{ROOT}/results/test_deployments.npz"))


def trace(a):
    lab, s = a
    proto, kw = CFG[lab]
    t = []
    simulate(proto, 300, s, xy=DEP[f"hsa_100_300_{s}"], rmax=100000, side=100, bs0=(5, 5), L=800, cov_trace=t, **kw)
    return lab, np.array(t)


if __name__ == "__main__":
    metrics("S1", "fig_s1_metrics.pdf")
    metrics("S3", "fig_s3_metrics.pdf")
    with Pool(50, initializer=init) as pool:
        tr = pool.map(trace, [(lab, s) for lab in CFG for s in range(1000, 1010)])
    fig, ax = plt.subplots(figsize=(3.5, 2.5))
    for lab, (c, ls, _) in STY.items():
        T = [t for l, t in tr if l == lab]
        n = max(len(t) for t in T)
        ax.plot(np.arange(n), np.mean([np.pad(t, (0, n - len(t))) for t in T], 0), ls, color=c, lw=1.4, label=lab)
    ax.set_xlabel("Round")
    ax.set_ylabel("Covered fraction of field")
    ax.legend(frameon=False, loc="lower left")
    fig.tight_layout()
    fig.savefig(f"{OUT}/fig_cov_time.pdf")
    fig.savefig(f"{ROOT}/figures/fig_cov_time.png", dpi=220)
    print("figures written to", OUT)
