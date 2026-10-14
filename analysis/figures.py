#!/usr/bin/env python3
"""
Publication figures, built from the per-run CSVs in ``results/``.

Design constraints, all deliberate:

* **Print-safe, colour-blind-safe palette.** The Okabe-Ito order, reordered so
  that no adjacent pair falls below the CVD separation floor. Validated rather
  than eyeballed.
* **Colour is never the only channel.** Every series also carries a distinct
  marker and line style, and is labelled directly on the plot, so the figures
  survive greyscale printing and photocopying.
* **No dual axes anywhere.** Where two quantities of different scale belong
  together (cost and accuracy, cost and compute), they get adjacent panels
  sharing an x-axis instead.
* **Medians with interquartile bands**, not means with standard errors: the
  experiment-count distributions are integer, bounded below and right-skewed.

Writes PDF (for LaTeX) and PNG (for quick viewing) into ``figures/``.
"""

from __future__ import annotations

import pathlib
import sys

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
from stats import load, where  # noqa: E402

FIGURES = pathlib.Path(__file__).resolve().parents[1] / "figures"
RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"

# Okabe-Ito, ordered so adjacent pairs clear the CVD separation floor.
BLUE, VERMILLION, GREEN, ORANGE, PURPLE, SKY = (
    "#0072B2", "#D55E00", "#009E73", "#E69F00", "#CC79A7", "#56B4E9",
)
INK, MUTED, GRID = "#1a1a1a", "#555555", "#d9d9d9"

STYLE = {
    "historical": dict(color=INK, marker="s", ls=(0, (1, 1)), label="Historical (2006)"),
    "adaptive": dict(color=BLUE, marker="o", ls="-", label="Adaptive (info gain)"),
    "informed": dict(color=GREEN, marker="^", ls="-.", label="Adaptive + prior"),
    "random": dict(color=PURPLE, marker="v", ls=(0, (3, 1, 1, 1)), label="Random + Bayes"),
    "half-removal": dict(color=VERMILLION, marker="D", ls="--", label="Half-removal (classical)"),
    "latent": dict(color=ORANGE, marker="P", ls="-", label="Latent size"),
    "two-stage-6": dict(color=SKY, marker="X", ls=":", label="Two-stage (6)"),
    "two-stage-12": dict(color=PURPLE, marker="x", ls=":", label="Two-stage (12)"),
    "fixed-4": dict(color=INK, marker="s", ls=(0, (1, 1)), label="Fixed size 4"),
}


def rc() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "font.size": 8.5,
        "axes.titlesize": 9,
        "axes.labelsize": 8.5,
        "legend.fontsize": 7.5,
        "axes.edgecolor": MUTED,
        "axes.linewidth": 0.7,
        "axes.grid": True,
        "grid.color": GRID,
        "grid.linewidth": 0.5,
        "axes.axisbelow": True,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "text.color": INK,
        "axes.labelcolor": INK,
        "figure.facecolor": "white",
        "savefig.bbox": "tight",
    })


def save(fig, name: str) -> None:
    FIGURES.mkdir(exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(FIGURES / f"{name}.{ext}")
    plt.close(fig)
    print(f"  wrote figures/{name}.pdf and .png")


def tidy(ax) -> None:
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)


def panel(ax, letter: str, note: str = "") -> None:
    """Label a panel A/B/... The interpretation belongs in the caption, not here:
    at journal column width long axes titles collide across panels."""
    ax.set_title(f"{letter}" + (f"   {note}" if note else ""), loc="left",
                 color=INK, fontweight="bold", fontsize=9)


# ------------------------------------------------------------------ figure 1


def figure1_main() -> None:
    """Cost distribution per strategy at 24 candidates, 300 replicates."""
    cols = load("05_scale_main.csv")
    order = ["historical", "random", "adaptive", "informed"]
    fig, ax = plt.subplots(figsize=(5.4, 2.6))

    rng = np.random.default_rng(0)
    for i, key in enumerate(order):
        exp = cols["n_experiments"][where(cols, strategy=key)]
        st = STYLE[key]
        # raw runs, jittered, so the reader sees the distribution not a summary
        ax.scatter(exp + rng.uniform(-0.3, 0.3, exp.size), np.full(exp.size, i)
                   + rng.uniform(-0.16, 0.16, exp.size),
                   s=3, color=st["color"], alpha=0.25, linewidths=0)
        q25, med, q75 = np.percentile(exp, [25, 50, 75])
        ax.plot([q25, q75], [i, i], color=st["color"], lw=2.4, solid_capstyle="round")
        ax.plot([med], [i], marker=st["marker"], color=st["color"], ms=6,
                markeredgecolor="white", markeredgewidth=0.8, zorder=5)
        ok = cols["correct"][where(cols, strategy=key)].mean()
        ax.annotate(f"median {med:.0f}   ·   {100 * ok:.1f}% correct",
                    (med, i + 0.30), color=INK, fontsize=7.5, ha="left", va="bottom")

    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([STYLE[k]["label"] for k in order])
    ax.set_xlabel("wet-lab experiments to identify the recipe")
    ax.set_xlim(0, 70)
    ax.set_ylim(-0.5, len(order) - 0.2)
    ax.grid(axis="y", visible=False)
    tidy(ax)
    panel(ax, "A", "24 candidates, 300 paired replicates")
    save(fig, "fig1_main_comparison")


# ------------------------------------------------------------------ figure 2


def _scaling_rows() -> list[dict]:
    """Experiment 8 (Bayesian arms) and 11 (classical arm) share hidden recipes and
    seeds, so they merge into one table rather than being re-run."""
    import csv as _csv

    rows = []
    for name in ("08_scaling.csv", "11_classical_baseline.csv"):
        path = RESULTS / name
        if not path.exists():
            continue
        for r in _csv.DictReader(path.open(newline="")):
            rows.append({
                "method": r.get("method") or r.get("strategy"),
                "n": int(float(r["n_candidates"])),
                "experiments": float(r["n_experiments"]),
                "correct": float(r["correct"]),
                "seconds": float(r.get("elapsed_s") or 0.0),
            })
    return rows


def figure2_scaling() -> None:
    """The headline: wet-lab cost and compute cost as the candidate list grows."""
    rows = _scaling_rows()
    if not rows:
        raise FileNotFoundError("results/08_scaling.csv and 11_classical_baseline.csv")

    def pick(method, n, key="experiments"):
        return np.array([r[key] for r in rows if r["method"] == method and r["n"] == n])

    sizes = np.array(sorted({r["n"] for r in rows}))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 2.8), layout="constrained")

    for key in ("historical", "half-removal", "adaptive", "informed"):
        xs, meds, los, his = [], [], [], []
        for n in sizes:
            exp = pick(key, n)
            if exp.size == 0:
                continue
            q25, med, q75 = np.percentile(exp, [25, 50, 75])
            xs.append(n), meds.append(med), los.append(q25), his.append(q75)
        st = STYLE[key]
        ax.fill_between(xs, los, his, color=st["color"], alpha=0.12, linewidth=0)
        ax.plot(xs, meds, color=st["color"], ls=st["ls"], marker=st["marker"],
                ms=4.5, lw=1.6, label=st["label"])

    # growth-rate reference: O(g^2 d log n), constant fitted at the smallest size
    ref = 4 * np.log(sizes)
    scale = float(np.median(pick("half-removal", sizes[0]))) / ref[0]
    ax.plot(sizes, scale * ref, color=MUTED, lw=1.0, ls=(0, (4, 2)))
    ax.annotate("$O(d\\log n)$", (sizes[-1], scale * ref[-1]), xytext=(-6, -12),
                textcoords="offset points", color=MUTED, fontsize=7, ha="right")

    ax.set_xlabel("candidate factors $n$")
    ax.set_ylabel("experiments (median, IQR band)")
    ax.set_xticks(sizes)
    ax.set_xlim(sizes[0] - 1, sizes[-1] + 2)
    ax.legend(frameon=False, fontsize=7, loc="upper left", handlelength=2.4,
              borderpad=0.1, labelspacing=0.35)
    tidy(ax)
    panel(ax, "A", "wet-lab cost")

    for key in ("half-removal", "adaptive", "informed"):
        xs, secs = [], []
        for n in sizes:
            s_n = pick(key, n, "seconds")
            if s_n.size == 0:
                continue
            xs.append(n), secs.append(max(float(np.median(s_n)), 1e-4))
        st = STYLE[key]
        bx.plot(xs, secs, color=st["color"], ls=st["ls"],
                marker=st["marker"], ms=4.5, lw=1.6, label=st["label"])
    bx.set_yscale("log")
    bx.set_xlabel("candidate factors $n$")
    bx.set_ylabel("seconds to choose the experiments")
    bx.set_xticks(sizes)
    bx.set_xlim(sizes[0] - 1, sizes[-1] + 2)
    # Upper left is the only region of panel B no curve passes through.
    bx.legend(frameon=False, fontsize=7, loc="upper left", handlelength=2.4,
              borderpad=0.1, labelspacing=0.35)
    tidy(bx)
    panel(bx, "B", "compute cost")

    save(fig, "fig2_scaling")


# ------------------------------------------------------------------ figure 3


def figure3_calibration() -> None:
    """Synergy strength vs calibration error, and how loudly failure announces itself."""
    cols = load("06_synergy_threshold.csv")
    grid = np.array(sorted(set(cols["oracle_partial_3"].tolist())))
    positive = grid[grid > 0]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 2.8), layout="constrained")

    for arm, color, marker, ls, label in (
        ("matched", GREEN, "^", "-", "likelihood knows the truth"),
        ("misspecified", VERMILLION, "D", "--", "likelihood assumes 0.10"),
    ):
        succ = [cols["correct"][where(cols, arm=arm, oracle_partial_3=p)].mean()
                for p in positive]
        ax.plot(positive, 100 * np.array(succ), color=color, marker=marker, ls=ls,
                ms=4.5, lw=1.6)
        offset = (-6, -18) if arm == "matched" else (-6, 14)
        ax.annotate(label, (positive[-1], 100 * succ[-1]), xytext=offset,
                    textcoords="offset points", color=color, fontsize=7, ha="right")
    ax.axvline(0.10, color=MUTED, lw=0.8, ls=(0, (1, 2)))
    ax.annotate("calibrated\nvalue", (0.10, 46), xytext=(3, 0),
                textcoords="offset points", color=MUTED, fontsize=7)
    ax.set_xscale("log")
    ax.set_xlabel("true three-of-four efficiency")
    ax.set_ylabel("recipe recovered (%)")
    ax.set_ylim(30, 104)
    tidy(ax)
    panel(ax, "A", "accuracy")

    # Reported confidence when the answer is wrong -- the silent-failure panel.
    ratios, confs = [], []
    for p in grid:
        mask = where(cols, arm="misspecified", oracle_partial_3=p)
        wrong = cols["correct"][mask] == 0
        if wrong.any():
            ratios.append(p / 0.10)
            confs.append(cols["final_confidence"][mask][wrong].mean())
    bx.scatter(np.maximum(ratios, 1e-2), confs, s=26, color=VERMILLION,
               marker="D", zorder=4, edgecolors="white", linewidths=0.6)
    bx.axhline(0.95, color=MUTED, lw=0.8, ls=(0, (1, 2)))
    bx.annotate("stopping threshold 0.95", (1e-2, 0.952), color=MUTED, fontsize=7,
                va="bottom")
    bx.set_xscale("log")
    bx.set_xlabel("true ÷ assumed synergy  (0.01 marks the zero-synergy case)")
    bx.set_ylabel("confidence reported when wrong")
    bx.set_ylim(0.78, 1.02)
    tidy(bx)
    panel(bx, "B", "confidence when wrong")

    save(fig, "fig3_calibration")


# ------------------------------------------------------------------ figure 4


def figure4_cardinality() -> None:
    """Unknown recipe size: four treatments, three truths."""
    cols = load("12_cardinality_baselines.csv")
    sizes = sorted(set(cols["true_size"].tolist()))
    methods = ["fixed-4", "two-stage-6", "two-stage-12", "latent", "half-removal"]
    methods = [m for m in methods if (cols["method"] == m).any()]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 2.8), layout="constrained")
    width = 0.8 / len(methods)

    for j, m in enumerate(methods):
        st = STYLE[m]
        succ = [100 * cols["correct"][where(cols, method=m, true_size=k)].mean()
                for k in sizes]
        cost = [np.median(cols["n_experiments"][where(cols, method=m, true_size=k)])
                for k in sizes]
        x = np.arange(len(sizes)) + j * width - 0.4 + width / 2
        ax.bar(x, succ, width=width * 0.86, color=st["color"], label=st["label"],
               edgecolor="white", linewidth=0.6)
        bx.bar(x, cost, width=width * 0.86, color=st["color"], edgecolor="white",
               linewidth=0.6)
        for xi, v in zip(x, succ):
            if v < 4:  # make the zeroes visible rather than invisible
                ax.annotate("0", (xi, 1.5), ha="center", fontsize=6.5, color=st["color"])

    for axis, ylabel, title in (
        (ax, "recipe recovered (%)", "A"),
        (bx, "experiments (median)", "B"),
    ):
        axis.set_xticks(range(len(sizes)))
        axis.set_xticklabels([f"{k:.0f} factors" for k in sizes])
        axis.set_xlabel("true recipe size")
        axis.set_ylabel(ylabel)
        axis.grid(axis="x", visible=False)
        tidy(axis)
        panel(axis, title)
    ax.set_ylim(0, 108)
    # One shared legend under both panels: inside panel A it collided with the
    # axis label at journal column width.
    fig.legend(*ax.get_legend_handles_labels(), loc="outside lower center",
               ncol=3, frameon=False, fontsize=7.5)

    save(fig, "fig4_cardinality")


# ------------------------------------------------------------------ figure 5


def figure5_stopping() -> None:
    """The price of certainty, and whether the posterior is honest."""
    cols = load("07_stopping_rule.csv")
    grid = np.array(sorted(set(cols["confidence_target"].tolist())))
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 2.8), layout="constrained")

    for key in ("adaptive", "informed"):
        st = STYLE[key]
        means = [cols["n_experiments"][where(cols, strategy=key, confidence_target=c)].mean()
                 for c in grid]
        succ = [100 * cols["correct"][where(cols, strategy=key, confidence_target=c)].mean()
                for c in grid]
        ax.plot(1 - grid, means, color=st["color"], ls=st["ls"], marker=st["marker"],
                ms=4.5, lw=1.6, label=st["label"])
        bx.plot(1 - grid, succ, color=st["color"], ls=st["ls"], marker=st["marker"],
                ms=4.5, lw=1.6)

    ax.legend(frameon=False, fontsize=7, loc="center right", handlelength=2.4,
              borderpad=0.1, labelspacing=0.35)
    for axis, ylabel, title in (
        (ax, "experiments (mean)", "A"),
        (bx, "recipe recovered (%)", "B"),
    ):
        axis.set_xscale("log")
        axis.invert_xaxis()
        axis.set_xlabel("1 − confidence threshold  (stricter to the right)")
        axis.set_ylabel(ylabel)
        tidy(axis)
        panel(axis, title)

    save(fig, "fig5_stopping_rule")


def figure6_real_data() -> None:
    """The oracle's gap assumption against two measured datasets."""
    import csv as _csv

    def rows(name):
        path = RESULTS / f"real_data_audit_{name}.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        return list(_csv.DictReader(path.open(newline="")))

    norman = [float(r["one_short"]) for r in rows("norman")]
    atlas_rows = rows("tfatlas")
    atlas = [float(r["one_short"]) for r in atlas_rows if r["size"] == "2"]
    triples = [(r["combination"], float(r["one_short"]))
               for r in atlas_rows if r["size"] == "3"]

    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.9, 2.9), layout="constrained",
                                 gridspec_kw={"width_ratios": [1.45, 1]})

    rng = np.random.default_rng(1)
    for i, (vals, color, label) in enumerate((
        (norman, BLUE, f"Norman 2019 · K562 · CRISPRa · {len(norman)} pairs"),
        (atlas, GREEN, f"Joung 2023 · hESC · TF overexpression · {len(atlas)} pairs"),
    )):
        v = np.array(vals)
        ax.scatter(v, np.full(v.size, i) + rng.uniform(-0.17, 0.17, v.size),
                   s=11, color=color, alpha=0.4, linewidths=0)
        q25, med, q75 = np.percentile(v, [25, 50, 75])
        ax.plot([q25, q75], [i, i], color=color, lw=2.6, solid_capstyle="round")
        ax.plot([med], [i], marker="o", ms=7, color=color, markeredgecolor="white",
                markeredgewidth=1, zorder=5)
        ax.annotate(f"median {med:.2f}", (med, i + 0.24), color=INK, fontsize=7.5,
                    ha="center")
        ax.annotate(label, (0.02, i + 0.40), fontsize=7, color=color, va="center")

    for _, value in triples:
        ax.plot([value], [1.62], marker="^", ms=7, color=ORANGE,
                markeredgecolor="white", markeredgewidth=0.8)
    if triples:
        ax.annotate("the two measured 3-factor combinations ("
                    + ", ".join(f"{v:.2f}" for _, v in triples) + ")",
                    (0.02, 1.62), fontsize=7, color=ORANGE, va="center")

    ax.axvline(0.10, color=VERMILLION, lw=1.4, ls=(0, (3, 2)))
    ax.annotate("oracle assumes 0.10", (0.13, -0.42), color=VERMILLION, fontsize=7.5)
    ax.set_xlim(0, 2.0)
    ax.set_ylim(-0.55, 1.80)
    ax.set_yticks([])
    ax.set_xlabel("one factor's effect ÷ the combination's effect")
    ax.grid(axis="y", visible=False)
    tidy(ax)
    panel(ax, "A", "measured, two systems")

    # Anchor sensitivity: the hESC subset has no unperturbed control.
    summary = list(_csv.DictReader((RESULTS / "real_data_summary_tfatlas.csv").open(newline="")))
    names = ["least-squares fit\n(biased to additivity)", "centroid of\nall conditions",
             "centroid of\nsingle-TF conditions"]
    y = np.arange(len(summary))
    med = [float(r["one_short_median"]) for r in summary]
    lo = [float(r["one_short_q25"]) for r in summary]
    hi = [float(r["one_short_q75"]) for r in summary]
    bx.barh(y, med, height=0.5, color=GREEN, alpha=0.85, edgecolor="white")
    bx.errorbar(med, y, xerr=[np.array(med) - lo, np.array(hi) - np.array(med)],
                fmt="none", ecolor=INK, elinewidth=1, capsize=3)
    bx.axvline(0.10, color=VERMILLION, lw=1.4, ls=(0, (3, 2)))
    bx.set_yticks(y)
    bx.set_yticklabels(names[: len(summary)], fontsize=7)
    bx.invert_yaxis()
    bx.set_xlabel("one-short (median, IQR)")
    bx.set_xlim(0, 1.3)
    bx.grid(axis="y", visible=False)
    tidy(bx)
    panel(bx, "B", "reference choice")

    save(fig, "fig6_real_data")


FIGURE_FUNCS = {
    "1": figure1_main,
    "2": figure2_scaling,
    "3": figure3_calibration,
    "4": figure4_cardinality,
    "5": figure5_stopping,
    "6": figure6_real_data,
}


def main() -> None:
    rc()
    wanted = sys.argv[1:] or list(FIGURE_FUNCS)
    for key in wanted:
        try:
            FIGURE_FUNCS[key]()
        except FileNotFoundError as exc:
            print(f"  skipping figure {key}: {exc}")
        except SystemExit as exc:
            print(f"  skipping figure {key}: {exc}")


if __name__ == "__main__":
    main()
