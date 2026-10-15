#!/usr/bin/env python3
"""
Analysis of the large-scale campaign (experiments 05-09).

Run everything, or one section at a time:

    python analysis/analyze.py            # all available sections
    python analysis/analyze.py 05 07      # just those

Each section reads its per-run CSV from ``results/`` and prints a report.
Nothing here re-runs a simulation, so the statistics can be revised without
spending another CPU-hour.
"""

from __future__ import annotations

import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

from stats import (  # noqa: E402
    bootstrap_ci,
    describe,
    load,
    loglog_slope,
    logistic_threshold,
    mcnemar_exact,
    paired_permutation_test,
    where,
    wilson_ci,
)

RESULTS = pathlib.Path(__file__).resolve().parents[1] / "results"
LABELS = {
    "historical": "Historical protocol (2006)",
    "adaptive": "Adaptive design (info gain)",
    "informed": "Adaptive + literature prior",
    "random": "Random pools + Bayes [ablation]",
}


def rule(title: str) -> None:
    print(f"\n{'=' * 78}\n  {title}\n{'=' * 78}")


def pct(x: float) -> str:
    return f"{100 * x:.1f}%"


# ------------------------------------------------------------------ exp 05


def section_05() -> None:
    cols = load("05_scale_main.csv")
    rule("Experiment 5 -- main comparison, 300 paired replicates")

    base = None
    print(f"  {'strategy':<32}{'median':>7}{'95% CI':>11}{'mean':>7}"
          f"{'IQR':>9}{'p90':>6}{'max':>6}   {'success (95% CI)':<24}")
    print("  " + "-" * 96)
    summaries = {}
    for key in ("historical", "adaptive", "informed", "random"):
        mask = where(cols, strategy=key)
        exp = cols["n_experiments"][mask]
        ok = cols["correct"][mask]
        d = describe(exp)
        lo, hi = bootstrap_ci(exp, np.median)
        s_lo, s_hi = wilson_ci(int(ok.sum()), int(ok.size))
        summaries[key] = {"exp": exp, "ok": ok, **d}
        if key == "historical":
            base = d["median"]
        print(
            f"  {LABELS[key]:<32}{d['median']:>7.0f}"
            f"{f'[{lo:.0f},{hi:.0f}]':>11}{d['mean']:>7.1f}"
            f"{f'{d_q(d)}':>9}{d['p90']:>6.0f}{d['max']:>6.0f}   "
            f"{f'{pct(ok.mean())} [{pct(s_lo)}, {pct(s_hi)}]':<24}"
        )

    print(f"\n  Saving vs the 2006 protocol (median {base:.0f} experiments):")
    for key in ("adaptive", "informed", "random"):
        med = summaries[key]["median"]
        print(f"    {LABELS[key]:<32}{100 * (1 - med / base):>6.0f}%"
              f"   ({base / med:.1f}x fewer experiments)")

    rule("Paired tests -- same hidden recipe, same noise seed, different strategy")
    pairs = (
        ("adaptive", "random", "is informative pool choice doing the work?"),
        ("informed", "adaptive", "does the literature prior add anything?"),
        ("adaptive", "historical", "vs the historical protocol"),
    )
    for a, b, question in pairs:
        ea, eb = summaries[a]["exp"], summaries[b]["exp"]
        diff, p = paired_permutation_test(ea, eb)
        wins = int((ea < eb).sum())
        ties = int((ea == eb).sum())
        print(f"\n  {LABELS[a]} vs {LABELS[b]}")
        print(f"    {question}")
        print(f"    mean difference {diff:+.1f} experiments   permutation p = {fmt_p(p)}")
        print(f"    cheaper in {wins}/{ea.size} runs, tied in {ties}")
        oa, ob = summaries[a]["ok"].astype(bool), summaries[b]["ok"].astype(bool)
        both = int((oa & ob).sum())
        only_a = int((oa & ~ob).sum())
        only_b = int((~oa & ob).sum())
        neither = int((~oa & ~ob).sum())
        if only_a or only_b:
            p_mc = mcnemar_exact(both, only_a, only_b, neither)
            print(f"    accuracy: {a} only {only_a}, {b} only {only_b}"
                  f"   McNemar p = {fmt_p(p_mc)}")

    rule("Tail behaviour -- the number a lab would actually plan against")
    for key in ("adaptive", "informed", "random"):
        exp = summaries[key]["exp"]
        ok = summaries[key]["ok"].astype(bool)
        print(f"  {LABELS[key]:<32}"
              f"worst run {exp.max():>3.0f}   "
              f"P(>=61, i.e. no better than 2006) = {pct(np.mean(exp >= base))}   "
              f"failures {int((~ok).sum())}")
        if (~ok).any():
            print(f"      failing runs used {np.median(exp[~ok]):.0f} experiments "
                  f"(median) vs {np.median(exp[ok]):.0f} for successes")

    # Does a failed run announce itself?
    rule("Is failure detectable from the inside?")
    for key in ("adaptive", "informed", "random"):
        mask = where(cols, strategy=key)
        conf = cols["final_confidence"][mask]
        ok = cols["correct"][mask].astype(bool)
        if (~ok).any():
            print(f"  {LABELS[key]:<32}mean reported confidence: "
                  f"{conf[ok].mean():.3f} when right, {conf[~ok].mean():.3f} when wrong")
        else:
            print(f"  {LABELS[key]:<32}no failures to compare")


def d_q(d: dict) -> str:
    return f"{d['q25']:.0f}-{d['q75']:.0f}"


def fmt_p(p: float) -> str:
    return "< 1e-4" if p < 1e-4 else f"{p:.4f}"


# ------------------------------------------------------------------ exp 06


def section_06() -> None:
    cols = load("06_synergy_threshold.csv")
    rule("Experiment 6 -- how much graded synergy does the method need?")
    grid = sorted(set(cols["oracle_partial_3"].tolist()))

    print("  three-of-four efficiency, and what the searcher believes about it")
    print(f"\n  {'partial':>8} │ {'matched: median':>16} {'success':>9} │"
          f" {'misspec: median':>16} {'success':>9}")
    print("  " + "-" * 72)
    table = {}
    for p in grid:
        line = f"  {p:>8.3f} │"
        for arm in ("matched", "misspecified"):
            mask = where(cols, arm=arm, oracle_partial_3=p)
            exp = cols["n_experiments"][mask]
            ok = cols["correct"][mask]
            table[(arm, p)] = (exp, ok)
            line += f" {np.median(exp):>16.0f} {pct(ok.mean()):>9} │"
        print(line)

    print("\n  Expected colonies from a three-of-four pool, for scale:")
    print("    partial 0.10 -> ~4 colonies (the calibrated value, i.e. the 2006 assay)")
    print("    partial 0.02 -> ~1 colony     partial 0.005 -> ~0.2 colonies")

    rule("Where the method breaks, and why")
    for arm in ("matched", "misspecified"):
        xs = np.array(sorted(p for (a, p) in table if a == arm and p > 0))
        success = np.concatenate([table[(arm, p)][1] for p in xs])
        x_rep = np.concatenate([np.full(table[(arm, p)][1].size, p) for p in xs])
        worst = min(table[(arm, p)][1].mean() for p in xs)
        t50, slope = logistic_threshold(x_rep, success, 0.5)
        t90, _ = logistic_threshold(x_rep, success, 0.9)

        def report(level: float, est: float) -> str:
            # A logistic fitted to data that never drops below the level has no
            # crossing point; reporting the extrapolated number would be a lie.
            if np.isnan(est) or not (xs.min() <= est <= xs.max()):
                return f"no {100 * level:.0f}% crossing inside the tested range"
            return f"{100 * level:.0f}% success at partial = {est:.4f} (~{45 * est:.1f} colonies)"

        print(f"  {arm}")
        print(f"    lowest success anywhere in the sweep: {pct(worst)}")
        print(f"    {report(0.5, t50)}")
        print(f"    {report(0.9, t90)}")

    rule("The real axis is not synergy strength -- it is calibration error")
    print("  The misspecified searcher always assumes partial = 0.10. Success as a")
    print("  function of how wrong that assumption is:\n")
    print(f"  {'true/assumed':>13} │ {'true partial':>12} │ {'success':>8} │"
          f" {'median':>7} │ {'reported confidence when wrong':>31}")
    print("  " + "-" * 84)
    for p in grid:
        exp, ok = table[("misspecified", p)]
        ratio = p / 0.10
        mask = where(cols, arm="misspecified", oracle_partial_3=p)
        conf = cols["final_confidence"][mask]
        wrong = ok.astype(bool) == False  # noqa: E712
        conf_wrong = f"{conf[wrong].mean():.3f}" if wrong.any() else "--"
        label = f"{ratio:.2f}x" if p > 0 else "0x"
        print(f"  {label:>13} │ {p:>12.3f} │ {pct(ok.mean()):>8} │"
              f" {np.median(exp):>7.0f} │ {conf_wrong:>31}")

    rule("How much of the failure is a modelling error rather than lost information?")
    print("  Paired by run: same hidden recipe and noise, differing only in whether")
    print("  the searcher's likelihood knows the true partial efficiency.\n")
    for p in grid:
        if p == 0:
            continue
        a_ok = table[("matched", p)][1].astype(bool)
        b_ok = table[("misspecified", p)][1].astype(bool)
        n = min(a_ok.size, b_ok.size)
        a_ok, b_ok = a_ok[:n], b_ok[:n]
        only_a = int((a_ok & ~b_ok).sum())
        only_b = int((~a_ok & b_ok).sum())
        p_mc = mcnemar_exact(int((a_ok & b_ok).sum()), only_a, only_b,
                             int((~a_ok & ~b_ok).sum()))
        gap = a_ok.mean() - b_ok.mean()
        flag = "  <-- calibration matters here" if p_mc < 0.05 else ""
        print(f"    partial {p:>6.3f}:  matched {pct(a_ok.mean()):>6}"
              f"   misspecified {pct(b_ok.mean()):>6}"
              f"   gap {100 * gap:>+5.0f} pts   McNemar p = {fmt_p(p_mc)}{flag}")


# ------------------------------------------------------------------ exp 07


def section_07() -> None:
    cols = load("07_stopping_rule.csv")
    rule("Experiment 7 -- the price of certainty, and is the posterior honest?")
    grid = sorted(set(cols["confidence_target"].tolist()))

    for strategy in ("adaptive", "informed"):
        print(f"\n  {LABELS[strategy]}")
        print(f"  {'target':>8} │ {'median':>7} {'mean':>7} {'p90':>6} │"
              f" {'success':>9} {'95% CI':>16} │ {'mean reported':>14} {'gap':>7}")
        print("  " + "-" * 86)
        for c in grid:
            mask = where(cols, strategy=strategy, confidence_target=c)
            exp = cols["n_experiments"][mask]
            ok = cols["correct"][mask]
            conf = cols["final_confidence"][mask]
            lo, hi = wilson_ci(int(ok.sum()), int(ok.size))
            gap = float(np.nanmean(conf)) - float(ok.mean())
            print(f"  {c:>8.4f} │ {np.median(exp):>7.0f} {exp.mean():>7.1f}"
                  f" {np.percentile(exp, 90):>6.0f} │ {pct(ok.mean()):>9}"
                  f" {f'[{pct(lo)}, {pct(hi)}]':>16} │"
                  f" {float(np.nanmean(conf)):>14.3f} {gap:>+7.3f}")

    rule("Calibration verdict")
    for strategy in ("adaptive", "informed"):
        over = []
        for c in grid:
            mask = where(cols, strategy=strategy, confidence_target=c)
            ok = cols["correct"][mask]
            lo, hi = wilson_ci(int(ok.sum()), int(ok.size))
            # Over-confident iff the stated threshold is above the plausible
            # range of the empirical success rate.
            if c > hi:
                over.append((c, ok.mean(), hi))
        if over:
            print(f"  {LABELS[strategy]}: over-confident at "
                  + ", ".join(f"{c:.4f} (actual {pct(a)})" for c, a, _ in over))
        else:
            print(f"  {LABELS[strategy]}: no threshold where the stated confidence "
                  "exceeds the observed success rate")

    rule("Marginal cost of certainty")
    for strategy in ("adaptive", "informed"):
        meds = {}
        for c in grid:
            mask = where(cols, strategy=strategy, confidence_target=c)
            meds[c] = float(np.median(cols["n_experiments"][mask]))
        steps = []
        keys = sorted(meds)
        for lo_c, hi_c in zip(keys, keys[1:]):
            steps.append(f"{lo_c:.4f}->{hi_c:.4f}: {meds[hi_c] - meds[lo_c]:+.0f}")
        print(f"  {LABELS[strategy]}: " + "   ".join(steps))


# ------------------------------------------------------------------ exp 08


def section_08() -> None:
    cols = load("08_scaling.csv")
    rule("Experiment 8 -- scaling in candidate-list size")
    sizes = sorted(set(cols["n_candidates"].tolist()))

    def n_hypotheses(n: int) -> int:
        from math import comb

        return comb(int(n), 4)

    print(f"  {'cands':>6} {'hypotheses':>12} │ {'historical':>11} │"
          f" {'adaptive':>9} {'success':>8} {'sec/run':>8} │"
          f" {'informed':>9} {'success':>8}")
    print("  " + "-" * 86)
    series = {k: [] for k in ("adaptive", "informed", "historical")}
    secs = []
    for n in sizes:
        row = f"  {n:>6.0f} {n_hypotheses(n):>12,} │"
        for key in ("historical", "adaptive", "informed"):
            mask = where(cols, strategy=key, n_candidates=n)
            exp = cols["n_experiments"][mask]
            ok = cols["correct"][mask]
            series[key].append(float(np.median(exp)))
            if key == "historical":
                row += f" {np.median(exp):>11.0f} │"
            else:
                row += f" {np.median(exp):>9.0f} {pct(ok.mean()):>8}"
                if key == "adaptive":
                    t = float(np.median(cols["elapsed_s"][mask]))
                    secs.append(t)
                    row += f" {t:>8.1f} │"
        print(row)

    rule("Growth rates")
    bits = [np.log2(n_hypotheses(n)) for n in sizes]
    for key in ("historical", "adaptive", "informed"):
        slope = np.polyfit(bits, series[key], 1)[0]
        print(f"  {LABELS[key]:<32}{slope:>6.2f} experiments per bit of "
              f"hypothesis-space entropy")
    print(f"\n  wet-lab cost vs candidates : exponent {loglog_slope(np.array(sizes), np.array(series['adaptive'])):.2f}"
          "   (1.00 would be linear, as the historical protocol is)")
    print(f"  compute cost vs hypotheses : exponent "
          f"{loglog_slope(np.array([n_hypotheses(n) for n in sizes]), np.array(secs)):.2f}"
          "   (design-step cost, seconds per run)")

    print("\n  Saving against the historical protocol at each size:")
    for n, h, a, i in zip(sizes, series["historical"], series["adaptive"], series["informed"]):
        print(f"    {n:>3.0f} candidates: adaptive {100 * (1 - a / h):>3.0f}%"
              f"   informed {100 * (1 - i / h):>3.0f}%")


# ------------------------------------------------------------------ exp 09


def section_09() -> None:
    cols = load("09_prior_quality.csv")
    rule("Experiment 9 -- when does a literature prior stop helping and start hurting?")
    strengths = sorted(set(cols["prior_strength"].tolist()))
    dropped = sorted(set(cols["prior_dropped"].tolist()))

    for k in sorted(set(cols["prior_informed"].tolist())):
        print(f"\n  Prior up-weights {k:.0f} factors "
              f"({k - 4:.0f} of them decoys when coverage is complete)")
        print(f"  {'strength':>9} │ " + " │ ".join(
            f"drop {d:.0f}: {'med':>4} {'succ':>5}" for d in dropped))
        print("  " + "-" * (12 + 17 * len(dropped)))
        for s in strengths:
            line = f"  {s:>9.1f} │ "
            cells = []
            for d in dropped:
                mask = where(cols, prior_strength=s, prior_dropped=d, prior_informed=k)
                exp = cols["n_experiments"][mask]
                ok = cols["correct"][mask]
                cells.append(f"{np.median(exp):>10.0f} {pct(ok.mean()):>5}")
            print(line + " │ ".join(cells))

    rule("The trade: a strong prior is cheap when right and dangerous when wrong")
    print(f"  {'strength':>9} │ {'median when complete':>21} │"
          f" {'median when 1 missing':>22} │ {'success when 1 missing':>23}")
    print("  " + "-" * 82)
    for s in strengths:
        m0 = where(cols, prior_strength=s, prior_dropped=0)
        m1 = where(cols, prior_strength=s, prior_dropped=1)
        e0 = cols["n_experiments"][m0]
        e1 = cols["n_experiments"][m1]
        ok1 = cols["correct"][m1]
        lo, hi = wilson_ci(int(ok1.sum()), int(ok1.size))
        print(f"  {s:>9.1f} │ {np.median(e0):>21.0f} │ {np.median(e1):>22.0f} │"
              f" {pct(ok1.mean())} [{pct(lo)}, {pct(hi)}]")

    print("\n  Failure rate by coverage, pooled over breadth, per strength:")
    for s in strengths:
        parts = []
        for d in dropped:
            mask = where(cols, prior_strength=s, prior_dropped=d)
            ok = cols["correct"][mask]
            parts.append(f"drop {d:.0f}: {pct(1 - ok.mean())}")
        print(f"    strength {s:>6.1f}   " + "   ".join(parts))


def section_10() -> None:
    cols = load("10_unknown_cardinality.csv")
    rule("Experiment 10 -- recipe size as a latent variable instead of an assumption")
    sizes = sorted(set(cols["true_size"].tolist()))

    print(f"  {'true size':>10} │ {'designer':>9} │ {'median':>7} {'mean':>7} │"
          f" {'exact recipe':>13} {'95% CI':>17} │ {'right size':>11} │ {'P(true size)':>13}")
    print("  " + "-" * 96)
    cells = {}
    for k in sizes:
        for designer in ("fixed-4", "latent"):
            mask = where(cols, true_size=k, designer=designer)
            exp = cols["n_experiments"][mask]
            ok = cols["correct"][mask]
            size_ok = cols["size_correct"][mask]
            p_true = cols["p_true_size"][mask]
            cells[(k, designer)] = (exp, ok, p_true)
            lo, hi = wilson_ci(int(ok.sum()), int(ok.size))
            print(f"  {k:>10.0f} │ {designer:>9} │ {np.median(exp):>7.0f} {exp.mean():>7.1f} │"
                  f" {pct(ok.mean()):>13} {f'[{pct(lo)}, {pct(hi)}]':>17} │"
                  f" {pct(size_ok.mean()):>11} │ {np.mean(p_true):>13.3f}")

    rule("Does the fix cost anything on the case it did not need to generalise for?")
    e_fixed, ok_fixed, _ = cells[(4, "fixed-4")]
    e_latent, ok_latent, _ = cells[(4, "latent")]
    diff, p = paired_permutation_test(e_latent, e_fixed)
    print(f"  Four-factor truth, the case the fixed designer was built for:")
    print(f"    fixed-4  median {np.median(e_fixed):.0f}, success {pct(ok_fixed.mean())}")
    print(f"    latent   median {np.median(e_latent):.0f}, success {pct(ok_latent.mean())}")
    print(f"    premium  {diff:+.1f} experiments   permutation p = {fmt_p(p)}")

    rule("Is the failure still silent?")
    print("  The danger was never the failure itself but that the posterior stayed")
    print("  confident through it. Reported confidence, split by whether the answer")
    print("  was right:\n")
    for k in sizes:
        for designer in ("fixed-4", "latent"):
            mask = where(cols, true_size=k, designer=designer)
            conf = cols["final_confidence"][mask]
            ok = cols["correct"][mask].astype(bool)
            right = f"{conf[ok].mean():.3f}" if ok.any() else "--"
            wrong = f"{conf[~ok].mean():.3f}" if (~ok).any() else "--"
            # For the latent designer, the size marginal is the diagnostic the
            # fixed designer cannot produce at all.
            p_true = cols["p_true_size"][mask]
            diag = (
                f"   P(true size) = {p_true[~ok].mean():.3f} when wrong"
                if designer == "latent" and (~ok).any()
                else ""
            )
            print(f"    size {k:.0f}, {designer:<8} confidence {right} when right, "
                  f"{wrong} when wrong{diag}")


def section_12() -> None:
    cols = load("12_cardinality_baselines.csv")
    rule("Experiment 12 -- is joint inference over recipe size worth its machinery?")
    sizes = sorted(set(cols["true_size"].tolist()))
    methods = [m for m in ("fixed-4", "two-stage-6", "two-stage-12", "latent",
                           "half-removal") if (cols["method"] == m).any()]

    print(f"  {'method':<14}" + "".join(f"{f'true size {k:.0f}':>22}" for k in sizes))
    print("  " + "-" * (14 + 22 * len(sizes)))
    for m in methods:
        line = f"  {m:<14}"
        for k in sizes:
            mask = where(cols, method=m, true_size=k)
            exp = cols["n_experiments"][mask]
            ok = cols["correct"][mask]
            lo, hi = wilson_ci(int(ok.sum()), int(ok.size))
            line += f"{f'{np.median(exp):.0f} exp · {pct(ok.mean())}':>22}"
        print(line)

    rule("What each comparator costs on the four-factor case it was built for")
    base = cols["n_experiments"][where(cols, method="fixed-4", true_size=4)]
    for m in methods:
        if m == "fixed-4":
            continue
        other = cols["n_experiments"][where(cols, method=m, true_size=4)]
        n = min(base.size, other.size)
        diff, p = paired_permutation_test(other[:n], base[:n])
        print(f"  {m:<14} {diff:+6.1f} experiments vs fixed-4   p = {fmt_p(p)}")

    rule("Where the two-stage baseline loses")
    print("  It commits to a size estimated from a fixed first-stage budget, and")
    print("  cannot revisit that decision. Share of runs whose committed size was")
    print("  wrong:\n")
    for m in [x for x in methods if x.startswith("two-stage")]:
        for k in sizes:
            mask = where(cols, method=m, true_size=k)
            est = cols["estimated_size"][mask]
            wrong = float(np.mean(est != k))
            print(f"    {m:<13} true size {k:.0f}: committed wrongly in {pct(wrong)} of runs")


# ------------------------------------------------------------------ exp 13


def section_13() -> None:
    cols = load("13_measured_gap.csv")
    rule("Experiment 13 -- the designer under measured gap structure")
    rho = cols["rho"][where(cols, method="adaptive-assumed", shape="gap")]
    print(f"  rho drawn from {rho.size} Norman ratios: median {np.median(rho):.2f},"
          f" IQR {np.percentile(rho, 25):.2f}-{np.percentile(rho, 75):.2f}")
    print(f"  the published oracle assumes 0.10\n")

    order = ("calibrated", "adaptive-assumed", "adaptive-matched",
             "half-removal", "half-removal-matched", "historical")
    names = {
        "calibrated": "Adaptive at the assumed rho=0.10",
        "adaptive-assumed": "Adaptive, published calibration",
        "adaptive-matched": "Adaptive, told the true rho",
        "half-removal": "Half-removal, assumed threshold",
        "half-removal-matched": "Half-removal, matched threshold",
        "historical": "Historical protocol (2006)",
    }
    for shape in ("gap", "geometric"):
        print(f"  fall-off shape: {shape}")
        print(f"    {'method':<34}{'median':>7}{'mean':>8}{'success':>20}{'silent':>9}")
        for method in order:
            mask = where(cols, method=method, shape=shape)
            if not mask.any():
                continue
            exp, ok = cols["n_experiments"][mask], cols["correct"][mask]
            lo, hi = wilson_ci(int(ok.sum()), int(ok.size))
            silent = float(np.mean(cols["silent_failure"][mask]))
            print(f"    {names[method]:<34}{np.median(exp):>7.0f}{exp.mean():>8.1f}"
                  f"{pct(ok.mean()) + f' [{100 * lo:.0f}-{100 * hi:.0f}]':>20}"
                  f"{pct(silent):>9}")
        print()

    # The headline comparison: what a reader of this repository actually gets.
    for shape in ("gap", "geometric"):
        a = where(cols, method="adaptive-assumed", shape=shape)
        b = where(cols, method="adaptive-matched", shape=shape)
        ok_a, ok_b = cols["correct"][a], cols["correct"][b]
        both = int(np.sum((ok_a == 1) & (ok_b == 1)))
        only_a = int(np.sum((ok_a == 1) & (ok_b == 0)))
        only_b = int(np.sum((ok_a == 0) & (ok_b == 1)))
        neither = int(np.sum((ok_a == 0) & (ok_b == 0)))
        p = mcnemar_exact(both, only_a, only_b, neither)
        print(f"  {shape}: published calibration {pct(ok_a.mean())} vs matched"
              f" {pct(ok_b.mean())} success, exact McNemar {fmt_p(p)}")

    worst = cols["final_confidence"][where(cols, method="adaptive-assumed", shape="gap")]
    worst = worst[~np.isnan(worst)]
    if worst.size:
        print(f"\n  A silent failure is not a hedged one: under the published"
              f" calibration the\n  mean final confidence on runs that ended wrong"
              f" was {np.mean(worst):.3f}.")


# ------------------------------------------------------------------ exp 14


def section_14() -> None:
    cols = load("14_norman_retrospective.csv")
    rule("Experiment 14 -- retrospective discovery on measured data")
    print("  Every query returns a measured Norman condition; nothing is simulated.\n")
    print(f"    {'arm':<16}{'found best':>12}{'median tests':>14}{'at 20 tests':>14}"
          f"{'at 60 tests':>14}")
    traces = {}
    for arm in ("random", "singles-first", "surrogate", "exhaustive"):
        mask = where(cols, arm=arm)
        if not mask.any():
            continue
        found = cols["found_best"][mask]
        hits = cols["n_to_best"][mask][found == 1]
        at20 = cols["best_at_20"][mask]
        at60 = cols["best_at_60"][mask]
        at20, at60 = at20[~np.isnan(at20)], at60[~np.isnan(at60)]
        traces[arm] = (hits, at20)
        print(f"    {arm:<16}{f'{int(found.sum())}/{found.size}':>12}"
              f"{np.median(hits) if hits.size else float('nan'):>14.0f}"
              f"{pct(np.median(at20)) if at20.size else '':>14}"
              f"{pct(np.median(at60)) if at60.size else '':>14}")

    if "surrogate" in traces and "singles-first" in traces:
        a, b = traces["surrogate"][1], traces["singles-first"][1]
        n = min(a.size, b.size)
        diff, p = paired_permutation_test(a[:n], b[:n])
        print(f"\n  At a 20-experiment budget the surrogate reaches"
              f" {100 * np.median(a):.0f}% of the best measured combination"
              f" against\n  {100 * np.median(b):.0f}% for the singles-first protocol:"
              f" difference {100 * diff:.0f} points, {fmt_p(p)}.")
        print("  Read this as a statement about the landscape as much as the method:"
              "\n  the measured one-short ratio of 0.77 is exactly what makes a"
              " singles-first\n  protocol nearly sufficient here, and what makes the"
              " simulated assay\n  (where a single factor yields nothing) the easier"
              " problem of the two.")


# ------------------------------------------------------------------ exp 15


def section_15() -> None:
    cols = load("15_proposal_ablation.csv")
    rule("Experiment 15 -- proposal generator ablation")
    print("  Scoring rule, posterior, oracle, recipes and seeds all held fixed;"
          "\n  only the set of pools the designer may score changes.\n")
    sizes = sorted(set(cols["n_candidates"].tolist()))

    classical = {}
    if (RESULTS / "11_classical_baseline.csv").exists():
        cl = load("11_classical_baseline.csv")
        for n in sorted(set(cl["n_candidates"].tolist())):
            mask = where(cl, method="half-removal", n_candidates=n)
            if mask.any():
                classical[n] = (
                    float(np.median(cl["n_experiments"][mask])),
                    float(np.mean(cl["correct"][mask])),
                )

    print(f"    {'cands':>6}│{'original':>10}{'success':>9}{'s/run':>8}"
          f"│{'structured':>12}{'success':>9}{'s/run':>8}│{'half-removal':>14}{'success':>9}")
    print("    " + "-" * 92)
    series = {"original": [], "structured": []}
    for n in sizes:
        row = f"    {n:>6.0f}│"
        for proposals in ("original", "structured"):
            mask = where(cols, proposals=proposals, n_candidates=n)
            exp, ok = cols["n_experiments"][mask], cols["correct"][mask]
            series[proposals].append(float(np.median(exp)))
            width = 10 if proposals == "original" else 12
            row += (f"{np.median(exp):>{width}.0f}{pct(ok.mean()):>9}"
                    f"{np.median(cols['elapsed_s'][mask]):>8.1f}│")
        if n in classical:
            med, ok = classical[n]
            row += f"{med:>14.0f}{pct(ok):>9}"
        print(row)

    rule("Growth rates (exponent of median cost against candidate count)")
    x = np.array(sizes, dtype=float)
    for proposals in ("original", "structured"):
        print(f"  adaptive, {proposals:<11} {loglog_slope(x, np.array(series[proposals])):>6.2f}")
    if classical:
        cx = np.array(sorted(classical), dtype=float)
        cy = np.array([classical[n][0] for n in sorted(classical)])
        print(f"  half-removal (exp 11)   {loglog_slope(cx, cy):>6.2f}")

    print("\n  Paired comparison at each size (same recipe, same seed):")
    for n in sizes:
        a = cols["n_experiments"][where(cols, proposals="original", n_candidates=n)]
        b = cols["n_experiments"][where(cols, proposals="structured", n_candidates=n)]
        k = min(a.size, b.size)
        # a is the original arm, b the structured one, so a positive mean
        # difference is experiments the structured proposals did not spend.
        diff, p = paired_permutation_test(a[:k], b[:k])
        print(f"    {n:>3.0f} candidates: structured saves {diff:>6.1f} experiments"
              f" on average, {fmt_p(p)}")

    # Whether the regime boundary the paper reports has actually moved.
    if classical:
        crossed = [n for n in sizes if n in classical
                   and series["structured"][sizes.index(n)] > classical[n][0]]
        if crossed:
            print(f"\n  Structured proposals still cost more than half-removal from"
                  f" {int(min(crossed))} candidates up.")
        else:
            print("\n  Structured proposals are at or below half-removal's cost at"
                  " every size tested.")


# ------------------------------------------------------------------ exp 16


def section_16() -> None:
    cols = load("16_latent_gap.csv")
    rule("Experiment 16 -- the gap response as a latent variable")
    print("  All arms use structured proposals, the same drawn rho and the same")
    print("  recipes, so neither the proposal generator nor the landscape differs.\n")

    names = {
        "fixed-assumed": "Fixed at the assumed 0.10",
        "fixed-oracle": "Fixed, told the true rho",
        "pilot-then-fixed": "Pilot (5 tests), then fixed",
        "latent-rho": "Latent rho",
        "latent-rho-shape": "Latent rho and shape",
    }
    order = list(names)
    for shape in ("gap", "geometric"):
        print(f"  true fall-off: {shape}")
        print(f"    {'method':<30}{'median':>7}{'mean':>7}{'success':>19}"
              f"{'silent':>8}{'|rho err|':>10}")
        for method in order:
            mask = where(cols, method=method, shape=shape)
            if not mask.any():
                continue
            exp, ok = cols["n_experiments"][mask], cols["correct"][mask]
            lo, hi = wilson_ci(int(ok.sum()), int(ok.size))
            err = cols["rho_abs_error"][mask]
            print(f"    {names[method]:<30}{np.median(exp):>7.0f}{exp.mean():>7.1f}"
                  f"{pct(ok.mean()) + f' [{100 * lo:.0f}-{100 * hi:.0f}]':>19}"
                  f"{pct(np.mean(cols['silent_failure'][mask])):>8}"
                  f"{np.median(err):>10.2f}")
        print()

    rule("Does inferring rho cost more than measuring it first?")
    for shape in ("gap", "geometric"):
        a = where(cols, method="pilot-then-fixed", shape=shape)
        b = where(cols, method="latent-rho", shape=shape)
        cost_a, cost_b = cols["n_experiments"][a], cols["n_experiments"][b]
        k = min(cost_a.size, cost_b.size)
        diff, p = paired_permutation_test(cost_a[:k], cost_b[:k])
        ok_a, ok_b = cols["correct"][a][:k], cols["correct"][b][:k]
        mc = mcnemar_exact(
            int(np.sum((ok_a == 1) & (ok_b == 1))),
            int(np.sum((ok_a == 1) & (ok_b == 0))),
            int(np.sum((ok_a == 0) & (ok_b == 1))),
            int(np.sum((ok_a == 0) & (ok_b == 0))),
        )
        print(f"  {shape}: latent costs {-diff:+.1f} experiments against the pilot"
              f" ({fmt_p(p)}), and recovers\n     {pct(ok_b.mean())} against"
              f" {pct(ok_a.mean())} (exact McNemar {fmt_p(mc)}).")

    rule("What the latent arm pays for not being told rho")
    for shape in ("gap", "geometric"):
        a = where(cols, method="fixed-oracle", shape=shape)
        b = where(cols, method="latent-rho", shape=shape)
        k = min(int(a.sum()), int(b.sum()))
        diff, p = paired_permutation_test(
            cols["n_experiments"][a][:k], cols["n_experiments"][b][:k]
        )
        print(f"  {shape}: {-diff:+.1f} experiments against an oracle that cannot"
              f" exist in practice, {fmt_p(p)}")

    rule("Recovering the assay model itself")
    for shape in ("gap", "geometric"):
        mask = where(cols, method="latent-rho", shape=shape)
        cover = cols["rho_in_ci"][mask]
        cover = cover[~np.isnan(cover)]
        err = cols["rho_abs_error"][mask]
        print(f"  {shape}: rho recovered to a median absolute error of"
              f" {np.median(err):.2f}, against a grid spacing of\n     0.20-0.25"
              f" near the measured range; the 90% credible interval over grid"
              f" cells\n     covers the true value in {pct(np.mean(cover))}"
              f" of runs.")
        shaped = where(cols, method="latent-rho-shape", shape=shape)
        p_shape = cols["p_true_shape"][shaped]
        p_shape = p_shape[~np.isnan(p_shape)]
        if p_shape.size:
            print(f"     With the shape latent too, posterior mass on the true"
                  f" shape is {np.median(p_shape):.2f}"
                  f" (median),\n     above 0.9 in {pct(np.mean(p_shape > 0.9))}"
                  f" of runs.")


SECTIONS = {
    "05": section_05,
    "06": section_06,
    "07": section_07,
    "08": section_08,
    "09": section_09,
    "10": section_10,
    "12": section_12,
    "13": section_13,
    "14": section_14,
    "15": section_15,
    "16": section_16,
}
FILES = {
    "05": "05_scale_main.csv",
    "06": "06_synergy_threshold.csv",
    "07": "07_stopping_rule.csv",
    "08": "08_scaling.csv",
    "09": "09_prior_quality.csv",
    "10": "10_unknown_cardinality.csv",
    "12": "12_cardinality_baselines.csv",
    "13": "13_measured_gap.csv",
    "14": "14_norman_retrospective.csv",
    "15": "15_proposal_ablation.csv",
    "16": "16_latent_gap.csv",
}


def main() -> None:
    wanted = sys.argv[1:] or list(SECTIONS)
    for key in wanted:
        if key not in SECTIONS:
            print(f"  unknown section {key!r}; choose from {', '.join(SECTIONS)}")
            continue
        if not (RESULTS / FILES[key]).exists():
            print(f"  skipping {key}: results/{FILES[key]} not present yet")
            continue
        SECTIONS[key]()
    print()


if __name__ == "__main__":
    main()
