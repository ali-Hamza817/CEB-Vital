"""Additional paper figures, all computed from simulations on the held-out test networks.

  fig_snapshot     network state of four protocols at the same round (S1, N = 300, seed 1000)
  fig_gain_heatmap T90 gain of CEB-ViTAL over the strongest always-on baseline, scenario x N
  fig_dynamics     residual energy, energy imbalance and sensing nodes over time (S1, N = 300, 10 seeds)
PDFs go to latex/figures/, PNGs to figures/ (for the repository README). Facts quoted in the
paper are written to latex/numbers_extra.tex.
"""
import json
import os
import sys
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

from common import NS, ROOT, best, gain, get

sys.path.insert(0, ROOT)
from ceb_sim import E0, covered_frac, incidence, simulate  # noqa: E402

PDF = os.path.join(ROOT, "latex", "figures")
PNG = os.path.join(ROOT, "figures")
plt.rcParams.update({"font.size": 8, "axes.titlesize": 8.5, "axes.labelsize": 8, "legend.fontsize": 7,
                     "xtick.labelsize": 7, "ytick.labelsize": 7, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42})
BLUE = LinearSegmentedColormap.from_list("blue_seq", ["#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
INK, MUTED, GRID = "#0b0b0b", "#8a8984", "#e4e3df"
CFG = {"LEACH (HSA)": ("LEACH", {}), "DEEC": ("DEEC", {}), "ViTAL": ("ViTAL", {}),
       "EB-ViTAL": ("EB-ViTAL", {}), "CEB-ViTAL": ("EB-ViTAL", {"sched": "cov"})}
STY = {"LEACH (HSA)": ("#2a78d6", "-"), "DEEC": ("#6250d6", "-"), "ViTAL": ("#eb6834", "-"),
       "EB-ViTAL": ("#1baf7a", "--"), "CEB-ViTAL": ("#1baf7a", "-")}
SNAP_N, SNAP_SEED, SNAP_ROUND = 300, 1000, 5500
DEP = None


def init():
    global DEP
    DEP = dict(np.load(f"{ROOT}/results/test_deployments.npz"))


def run(a):
    lab, N, s, snap = a
    proto, kw = CFG[lab]
    et = []
    o = simulate(proto, N, s, xy=DEP[f"hsa_100_{N}_{s}"], rmax=100000, side=100, bs0=(5, 5), L=800,
                 snap_round=snap, e_trace=et if snap is None else None, **kw)
    return lab, s, o.get("snapshot"), np.array(et) if et else None


def save(fig, name):
    fig.savefig(f"{PDF}/{name}.pdf")
    fig.savefig(f"{PNG}/{name}.png", dpi=220)
    plt.close(fig)


def snapshot(snaps, facts):
    fig, axs = plt.subplots(2, 2, figsize=(7.16, 7.0))
    for ax, lab in zip(axs.ravel(), ["LEACH (HSA)", "ViTAL", "EB-ViTAL", "CEB-ViTAL"]):
        sn = snaps[lab]
        xy, E, alive, part, ch = sn["xy"], sn["E"], sn["alive"], sn["part"], sn["ch"]
        cov = covered_frac(incidence(xy, 100), part)
        for m, h in zip(sn["mem"], sn["head"]):                     # member -> CH links
            if h >= 0:
                ax.plot(*zip(xy[m], xy[h]), color=GRID, lw=0.5, zorder=1)
        dead = ~alive
        sleep = alive & ~part
        ax.scatter(*xy[dead].T, marker="x", s=9, lw=0.6, color=MUTED, zorder=2)
        ax.scatter(*xy[sleep].T, s=13, facecolor="white", edgecolor=MUTED, lw=0.6, zorder=2)
        awake = part.copy()
        awake[ch] = False
        sc = ax.scatter(*xy[awake].T, c=np.clip(E[awake], 0, None), cmap=BLUE, vmin=0, vmax=E0, s=16,
                        edgecolor="white", lw=0.3, zorder=3)
        ax.scatter(*xy[ch].T, c=np.clip(E[ch], 0, None), cmap=BLUE, vmin=0, vmax=E0, marker="*", s=120,
                   edgecolor=INK, lw=0.6, zorder=4)
        ax.scatter(*sn["bs"], marker="s", s=60, color=INK, zorder=5)
        ax.set_xlim(-2, 102)
        ax.set_ylim(-2, 102)
        ax.set_aspect("equal")
        ax.set_xticks([0, 50, 100])
        ax.set_yticks([0, 50, 100])
        ax.grid(False)
        ax.set_title(lab, fontweight="bold", pad=15)
        ax.text(0.5, 1.015, f"alive {alive.sum()}  ·  sensing {part.sum()}  ·  CHs {len(ch)}  ·  coverage {100*cov:.0f}%",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=7, color=INK)
        facts[lab] = dict(alive=int(alive.sum()), sensing=int(part.sum()), chs=len(ch), cov=cov,
                          meanE=float(np.clip(E[alive], 0, None).mean()))
    handles = [Line2D([], [], marker="o", ls="", color="#3987e5", label="sensing node (fill = residual energy)"),
               Line2D([], [], marker="*", ls="", color="#3987e5", mec=INK, ms=10, label="cluster head"),
               Line2D([], [], marker="o", ls="", mfc="white", mec=MUTED, label="sleeping node"),
               Line2D([], [], marker="x", ls="", color=MUTED, label="dead node"),
               Line2D([], [], marker="s", ls="", color=INK, label="sink")]
    fig.legend(handles=handles, loc="upper center", ncol=5, frameon=False, fontsize=7)
    cax = fig.add_axes([0.25, 0.06, 0.5, 0.012])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.set_label("Residual energy (J)")
    cb.outline.set_visible(False)
    fig.subplots_adjust(left=0.06, right=0.98, top=0.90, bottom=0.14, hspace=0.30, wspace=0.12)
    save(fig, "fig_snapshot")


def heatmap(facts):
    scen = ["S1", "S2", "S3", "S4", "S5", "S6"]
    names = ["S1  100 m, corner sink", "S2  100 m, centre sink", "S3  200 m, corner sink",
             "S4  200 m, centre sink", "S5  100 m, mobile sink", "S6  100 m, 4000-bit packets"]
    G = np.array([[gain(get(sc, "CEB-ViTAL", N, "T90"), get(sc, best(sc, N), N, "T90")) for N in NS] for sc in scen])
    fig, ax = plt.subplots(figsize=(7.16, 2.5))
    im = ax.imshow(G, cmap=BLUE, vmin=0, vmax=max(150, G.max()), aspect="auto")
    for i in range(G.shape[0]):
        for j in range(G.shape[1]):
            ax.text(j, i, f"+{G[i, j]:.0f}%", ha="center", va="center", fontsize=8,
                    color="white" if G[i, j] > 75 else INK)
    ax.set_xticks(range(len(NS)), NS)
    ax.set_yticks(range(len(scen)), names)
    ax.set_xlabel("Number of nodes $N$")
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("$T_{90}$ gain over strongest baseline (%)")
    cb.outline.set_visible(False)
    fig.tight_layout()
    save(fig, "fig_gain_heatmap")
    facts["heat_min"], facts["heat_max"] = float(G.min()), float(G.max())


def dynamics(traces, facts):
    fig, axs = plt.subplots(1, 3, figsize=(7.16, 2.4))
    for lab, (c, ls) in STY.items():
        T = traces[lab]
        n = min(len(t) for t in T)
        A = np.stack([t[:n] for t in T])                              # seeds x rounds x 5
        r = A[0, :, 0]
        axs[0].plot(r, A[:, :, 1].mean(0), ls, color=c, lw=1.4, label=lab)
        ok = (A[:, :, 1] > 0) & (A[:, :, 3] >= 0.1 * SNAP_N)          # at least 10% of nodes alive
        cv = np.where(ok, A[:, :, 2] / np.maximum(A[:, :, 1], 1e-12), np.nan)
        axs[1].plot(r, np.nanmean(cv, 0), ls, color=c, lw=1.4)
        axs[2].plot(r, A[:, :, 4].mean(0), ls, color=c, lw=1.4)
        facts.setdefault("sensing_mid", {})[lab] = float(A[:, 1000:3000, 4].mean())
    axs[0].set_title("Mean residual energy (J)")
    axs[1].set_title("Energy imbalance (std / mean)")
    axs[2].set_title("Sensing nodes per round")
    axs[1].set_ylim(0, 1.0)
    for ax in axs:
        ax.set_xlabel("Round")
        ax.grid(True, color=GRID, lw=0.5)
    fig.legend(*axs[0].get_legend_handles_labels(), loc="upper center", ncol=5, frameon=False)
    fig.tight_layout(rect=(0, 0, 1, 0.88))
    save(fig, "fig_dynamics")


if __name__ == "__main__":
    os.makedirs(PDF, exist_ok=True)
    os.makedirs(PNG, exist_ok=True)
    jobs = [(lab, SNAP_N, SNAP_SEED, SNAP_ROUND) for lab in ("LEACH (HSA)", "ViTAL", "EB-ViTAL", "CEB-ViTAL")]
    jobs += [(lab, SNAP_N, s, None) for lab in CFG for s in range(1000, 1010)]
    with Pool(50, initializer=init) as pool:
        out = pool.map(run, jobs)
    snaps = {lab: sn for lab, s, sn, _ in out if sn is not None}
    traces = {}
    for lab, s, _, et in out:
        if et is not None:
            traces.setdefault(lab, []).append(et)
    facts = {}
    snapshot(snaps, facts)
    heatmap(facts)
    dynamics(traces, facts)
    sm = facts["sensing_mid"]
    with open(os.path.join(ROOT, "latex", "numbers_extra.tex"), "w") as f:
        f.write("% Auto-generated by paper/make_extra_figures.py. Do not edit by hand.\n")
        f.write(f"\\newcommand{{\\SnapRound}}{{{SNAP_ROUND}}}\n\\newcommand{{\\SnapN}}{{{SNAP_N}}}\n")
        for key, lab in (("Leach", "LEACH (HSA)"), ("Vital", "ViTAL"), ("Eb", "EB-ViTAL"), ("Ceb", "CEB-ViTAL")):
            d = facts[lab]
            f.write(f"\\newcommand{{\\Snap{key}Alive}}{{{d['alive']}}}\n\\newcommand{{\\Snap{key}Sensing}}{{{d['sensing']}}}\n")
            f.write(f"\\newcommand{{\\Snap{key}Cov}}{{{100*d['cov']:.0f}\\%}}\n\\newcommand{{\\Snap{key}E}}{{{d['meanE']:.2f}}}\n")
        f.write(f"\\newcommand{{\\SensCeb}}{{{sm['CEB-ViTAL']:.0f}}}\n\\newcommand{{\\SensEb}}{{{sm['EB-ViTAL']:.0f}}}\n")
        f.write(f"\\newcommand{{\\HeatMin}}{{+{facts['heat_min']:.0f}\\%}}\n\\newcommand{{\\HeatMax}}{{+{facts['heat_max']:.0f}\\%}}\n")
    json.dump({k: v for k, v in facts.items()}, open(os.path.join(ROOT, "results", "figure_facts.json"), "w"), indent=1, default=float)
    print(json.dumps(facts, indent=1, default=float))
