# yamanaka-replay

**Replaying the 2006 discovery of the Yamanaka factors in silico, under its original fitness function, to ask a counterfactual question: how many wet-lab experiments would adaptive experimental design have needed?**

Answer, on a simulator calibrated to the published figures: **a median of 6–10 experiments, against 61 for the historical leave-one-out protocol** — measured over 300 independent replicates with a 95% bootstrap interval of a single experiment.

> **This is a methodological study on a simulator, not a biological result.**
> The oracle is a calibrated model of a selection assay, not an assay. Nothing
> here supports any claim about cellular biology.

The repository also reports where the method fails, because that turned out to
be the more useful half. A campaign of 6,320 replicates (`./campaign.sh`)
**refuted two of this project's own earlier claims** and removed its most
dangerous failure mode. An audit against two measured perturbation datasets
then **refuted the oracle's central calibration** as well. All three
corrections are documented below rather than quietly dropped.

A manuscript drawing these together is in [`paper/`](paper/).

---

## Why cost is measured in experiments

Takahashi & Yamanaka did not optimise a continuous score. Their readout was the
number of drug-resistant colonies surviving Fbx15/G418 selection: an **integer,
stochastic, and very expensive** signal, where a single evaluation is a wet-lab
experiment taking weeks.

Most combinatorial-search benchmarks report compute budget. That is the wrong
currency here. Under this fitness function the binding constraint is the
experiment count, so that is what this repository measures — and what the
historical record lets us compare against, because the 2006 protocol's
experiment count is fully documented.

## Oracle calibration

Colony counts are Poisson draws around a rate set by three effects reported in
the paper:

| Input pool | Expected colonies | Published basis |
|---|---|---|
| All 24 candidate factors | 22 | Figure reported in the paper |
| 10-factor subset | 36 | "The 10 outperformed the 24" |
| The 4-factor recipe | 45 | Clean recipe |
| Any single factor | 0 | No single factor sufficed |
| 3 of the 4 | ~4 | Triplets gave partial colonies |

Two features matter more than the absolute numbers. **Graded synergy**: pools
with three of the four produce a weak but nonzero signal. **A pool-size
penalty**: adding factors lowers per-cell efficiency, which is why the paper's
10-factor subset beat the full 24 — a counter-intuitive result the simulator
reproduces rather than assumes.

## Results

Headline figures below are from the 300-replicate campaign; `./reproduce.sh`
runs a 20-replicate version of the same comparison in under a minute and lands
within noise of these numbers. Raw per-run CSVs are in [`results/`](results/).

### Main comparison — 300 replicates, paired

Every strategy sees the same 300 hidden recipes and the same noise seeds, so
each comparison is paired run by run.

| Strategy | Median | 95% CI | Mean | IQR | p90 | Worst | Success |
|---|---|---|---|---|---|---|---|
| Historical protocol (2006) | 61 | [61, 61] | 61.0 | 61–61 | 61 | 61 | 100% |
| Adaptive design (information gain) | **10** | [10, 10] | 10.3 | 9–11 | 13 | **18** | 99.7% |
| Adaptive + literature prior | **6** | [6, 6] | 6.2 | 6–7 | 7 | **11** | 99.3% |
| Random pools + Bayes *(ablation)* | 44 | [42, 46] | 43.0 | 33–56 | 60 | 60 | 88.0% |

**The ablation is the load-bearing control.** Identical Bayesian machinery with
randomly chosen pools needs 44 experiments and succeeds 88% of the time. The
saving comes from *informative experiment selection*, not from Bayesian
updating: paired over 300 runs the adaptive designer is cheaper by 32.7
experiments on average (permutation *p* < 1e-4) and wins in 298 of 300 runs.

The literature prior buys speed but not accuracy: −4.1 experiments
(*p* < 1e-4), with no detectable difference in success rate (McNemar *p* = 1.0).

Against the historical protocol the adaptive designer is cheaper in **300 out
of 300 runs**, and its worst run over 300 was 18 experiments.

### The stopping rule is nearly free, and the posterior is honest

A Bayesian stopping rule is only usable if its confidence means something, so
the threshold was swept over 1,680 runs:

| Confidence target | Median | Mean | Observed success | Mean reported confidence |
|---|---|---|---|---|
| 0.50 | 10 | 9.9 | 90.8% | 0.924 |
| 0.90 | 10 | 10.2 | 95.8% | 0.983 |
| 0.95 | 10 | 10.4 | 99.2% | 0.995 |
| 0.999 | 10 | 10.6 | 100% | 1.000 |
| 0.9999 | 10 | **10.7** | 100% | 1.000 |

Two things follow. The posterior is **calibrated**: at no threshold does the
stated confidence exceed the observed success rate by more than its Wilson
interval allows (largest gap +0.025). And tightening the threshold by four
orders of magnitude costs **0.8 experiments on average** and never moves the
median, so there is no reason to run this method at a loose threshold.

That same flatness explains an uncomfortable property: confidence does not
climb gradually, it **jumps** at the decisive experiment. There is no
intermediate regime in which the searcher can report "not sure yet".

### Where the method breaks

The designer assumes a likelihood model. Making the oracle disagree with it, one
assumption at a time (`03_misspecification.py`, 20 runs per row):

| Oracle reality vs designer's assumption | Median | Success |
|---|---|---|
| Matches | 11 | 100% |
| Pool-size penalty 3× harsher | 12 | 100% |
| Triplet effect stronger than assumed | 13 | 79% |
| Overall efficiency 3× lower | 18 | 43% |
| Triplets have no effect at all | 25 | 43% |
| **Truth is five factors, not four** | 26 | **0%** |

Experiments 6 and 10 then took each of these apart, and the conclusions changed.

#### Correction 1: absent graded synergy is a calibration error, not a structural one

This README previously listed "no graded synergy" as a structural failure the
method cannot survive, on the strength of the 43% row above. That experiment
conflated two things: it removed the partial signal from the *assay* while
leaving the searcher *believing* in it.

Experiment 6 separates them, pairing by run so that only the searcher's
assumption differs:

| True three-of-four efficiency | Searcher told the truth | Searcher assumes 0.10 |
|---|---|---|
| **0.000** (no synergy at all) | 22 experiments, **100%** | 35 experiments, 52% |
| 0.010 | 19, **100%** | 18, 70% |
| 0.050 | 12, **100%** | 12, 94% |
| 0.100 (calibrated) | 10, 98% | 10, 98% |
| 0.400 | 10, **100%** | 12, **42%** |

With no graded synergy whatsoever and a searcher that knows it, the method
succeeds **100% of the time** in a median of 22 experiments — it degenerates
into pure group testing, pays 2.2× for the lost gradient, and is still 64%
cheaper than the 2006 protocol. The earlier 43% was misspecification, not lost
information.

What does matter is **how wrong the assumed synergy is, and in which
direction**. The tolerance window is asymmetric: under-estimating the partial
effect by 5× still succeeds 70% of the time, while over-estimating it by 4×
drops to 42% — and reports a confidence of 0.989 while doing so. **Calibrate
the partial effect conservatively.**

#### Correction 2: the assumed gap response is far sharper than measurement supports

The oracle assumes a pool one factor short yields 10% of the full response. The
empirical analogue — for a measured combination, how much of its effect one of
its constituents reproduces alone — was computed on two public datasets:

| Dataset | System | Combinations | Median one-short | Below 0.20 |
|---|---|---|---|---|
| [Norman 2019](https://doi.org/10.1126/science.aax4438) | K562, CRISPRa | 131 pairs | **0.77** | 0 of 131 |
| [Joung 2023](https://pubmed.ncbi.nlm.nih.gov/36608654/) | hESC, TF overexpression | 40 pairs | **0.77** | 0 of 40 |

The two measured three-factor combinations give 0.80 and 1.01. Threshold-like
effect measures (counting genes that move past a cutoff rather than a Euclidean
norm) lower the Norman median to 0.48 and put 25% of pairs below 0.20, but no
measure produces a median near 0.10.

Since over-estimating the gap response is the error direction that fails
silently (Correction 1), **our own calibration sits on the dangerous side**. The
same designer informed of the true value succeeds in 98–100% of runs at every
value tested, so the operational rule is to estimate the gap response in a
pilot rather than assume it. Caveat: both datasets read out the transcriptome,
while the 2006 assay reads out a rare discrete phenotype, which may be far
sharper. Reproduce with `python analysis/real_data_audit.py`.

#### Correction 3: the design step does not scale as claimed

This README previously said the per-experiment design cost is independent of
the hypothesis-space size. That is too strong. Experiment 8, over candidate
lists from 24 to 60 factors:

| Candidates | Hypotheses | Historical | Adaptive | Design cost (s/run) | With literature prior |
|---|---|---|---|---|---|
| 24 | 10,626 | 61 | 10 | 0.3 | 6 |
| 32 | 35,960 | 77 | 16 | 2.1 | 6 |
| 40 | 91,390 | 93 | 19 | 5.5 | 6 |
| 48 | 194,580 | 109 | 27 | 17.4 | 6 |
| 60 | 487,635 | 133 | **58** | **159.2** | 6 |

Wet-lab cost grows as roughly `n^1.8` in the candidate count, not
logarithmically, and the saving decays from 84% to 57%. A classical model-free
conjunctive search (`half_removal_search`) settles the cause: it grows as
`n^0.4`, and the two cross at 48 candidates (26 versus 27 experiments). The
bottleneck is not evaluating candidate pools but *proposing* good ones: random
pools of size 2–12 hit all four factors with probability falling as the fourth
power of `n`. That is a fixable defect in the proposal generator rather than in
the theory. Neither method dominates: the classical search's success rate falls
to 83% at 60 candidates, because the pool-size penalty drags its positive
control toward the noise floor, while the Bayesian designer holds 100%
throughout at five orders of magnitude more compute.

Compute cost grows as roughly `H^1.5` in the hypothesis count. The correct
narrow claim is that the design step is independent of **how concentrated the
posterior is** — the five-overlap-class collapse is real — but it remains linear
in the **size of the enumerated space**, so exact enumeration is the binding
constraint past about 60 candidates.

The flat right-hand column is **not** a result. The literature prior up-weights
8 factors regardless of `n`, so the effective space never grows. A longer
candidate list in reality means a *less* informative prior, so treat that column
as an upper bound conditional on prior quality holding, not as a prediction.

#### What survived: fixed recipe cardinality — and it is now fixed

The one genuinely structural failure is assuming the wrong recipe size, and it
is worse than previously reported. Experiment 10, 40 runs per cell, with
`GeneralOracle` hiding recipes of three, four and five factors:

| True recipe size | Designer | Median | Exact recipe recovered | P(true size) inferred |
|---|---|---|---|---|
| 3 | fixed at four | 9 | **0%** | — |
| 3 | **size latent** | 9 | **100%** | 0.999 |
| 4 | fixed at four | 10 | 100% | — |
| 4 | **size latent** | 11 | **100%** | 0.997 |
| 5 | fixed at four | 18 | **0%** | — |
| 5 | **size latent** | 16 | **100%** | 0.995 |

Two findings. First, the failure is **symmetric**: a three-factor truth is as
fatal as a five-factor one. The fixed designer finds the three real factors,
appends a spurious fourth, and reports 0.999 confidence. The defect is fixing
the cardinality at all, in either direction.

Second, [`yamanaka_replay/cardinality.py`](yamanaka_replay/cardinality.py)
removes it. Recipe size becomes a latent variable: the hypothesis space is every
subset of size 3–6 (189,750 of them), and the rate depends on a hypothesis only
through the pair `(|h|, |pool ∩ h|)` — 22 reachable classes — so the same
collapse that makes the fixed-size designer cheap still applies one dimension
up. The premium on the four-factor case it did not need to generalise for is
**+0.9 experiments** (paired, *p* = 0.0036).

The latent designer also returns a posterior over recipe size, which is the
diagnostic the fixed designer cannot produce even in principle: a model that
assumes four factors can never doubt the four.

### What the search actually does

From `experiments/04_trace_single_run.py` (hidden recipe `F03, F04, F11, F17`):

```
 exp │ size │ colonies │ best hypothesis    │ conf
 ────┼──────┼──────────┼────────────────────┼──────
   1 │  12  │     0    │ F00+F01+F02+OCT4   │ 0.000
   4 │  12  │     1    │ F00+F01+F06+F07    │ 0.001
   6 │  12  │     3    │ F00+F01+F06+F07    │ 0.006
   9 │   8  │     6    │ F00+F05+F07+F14    │ 0.067
  10 │  12  │    36    │ OCT4+SOX2+F07+KLF4 │ 0.197
  11 │  11  │    31    │ OCT4+SOX2+KLF4+MYC │ 1.000
```

Three things are visible. Zero-colony pools are not failures — each eliminates
every hypothesis it contains, which is classical group testing. The weak
three-factor signal at experiments 4–9 is the thread the search pulls on. And
the designer never once tests a single factor, where the historical protocol
spent 24 of its 61 experiments.

### Failure is quiet, and efficiency makes it quieter

Across the campaign the same pattern appears in three unrelated designs. When
the adaptive designer is wrong, its mean reported confidence is 0.955; when the
ablation is wrong, it is 0.423. Over-estimated synergy fails at 0.989. A
cardinality mismatch fails at 0.999.

The designer is not miscalibrated — experiment 7 shows it is well calibrated
*inside its model*. The point is narrower and more awkward: **efficient
experiment selection improves cost and degrades the detectability of model
error at the same time**, because a method that gathers consistent evidence
quickly also gathers it quickly for the wrong answer. This is the finding most
likely to generalise beyond this benchmark.

## Related work, and what is actually new here

Adaptive and Bayesian experimental design for biological perturbation selection
is **an established field**, and this repository does not claim to introduce it:

- [Zhang, Uhler & Sapsis, *Nature Machine Intelligence* (2023)](https://www.nature.com/articles/s42256-023-00719-0)
  — active learning over genetic interventions, including transcription-factor
  combinations, evaluated on a simulated cellular-reprogramming task. The
  closest prior work to this one.
- [BATCHIE, *Nature Communications* (2024)](https://www.nature.com/articles/s41467-024-55287-7)
  — Bayesian active learning for combination drug screens; retrospective
  simulation over three public screens, reporting large experiment savings.
- [BioDiscoveryAgent (2024)](https://arxiv.org/html/2405.17631v3) — LLM-driven
  closed-loop design of genetic perturbation experiments.
- [Mogrify, *Nature Genetics* (2016)](https://doi.org/10.1038/ng.3487) —
  predicts reprogramming factors from expression and network data; a ranking
  method, not a sequential design method.

Against that background, the contribution of this repository is narrow and
specific:

1. **A counterfactual replay of one documented historical discovery.** Prior
   retrospective studies replay adaptive data collection over existing datasets
   and report predictive accuracy. This replays the path to a single known
   discovery under its own original assay, and compares against that discovery's
   *actual, published* experiment count.
2. **A reproducible benchmark, and a campaign that stress-tests it.** 24
   candidates, 10,626 hypotheses, an oracle pinned to published numbers, four
   strategies including an ablation, one dependency, under a minute on a laptop
   — plus 6,320 replicates isolating synergy from calibration, pricing the
   stopping rule, and measuring scaling.
3. **A symmetric, silent failure under fixed recipe cardinality, and its
   removal.** Fixing the recipe size at the wrong value fails at 0% success and
   0.999 confidence in *both* directions; treating size as a latent variable
   restores 100% success for a premium of one experiment. We have not seen this
   failure mode reported or addressed elsewhere.
4. **Evidence that search efficiency and error detectability trade off.**
   Three independent experimental designs here show the better searcher failing
   more confidently.

Claims 3 and 4 are the parts most likely to be useful to others. Claims 1 and 2
are infrastructure.

## Reproducing

```bash
git clone https://github.com/mashal-alkhalifah/yamanaka-replay
cd yamanaka-replay
pip install -r requirements.txt
./reproduce.sh      # every README figure from experiments 01-04, under a minute
```

The large-scale campaign behind the intervals, the paired tests and the
corrections above is a separate, much longer run:

```bash
./campaign.sh       # experiments 05-10, 6,320 replicates, ~2 hours on 2 cores
python analysis/analyze.py      # re-analyse the per-run CSVs without re-running
```

**No GPU, no network access, no external data.** Reference environment: Linux
6.18, x86_64, 2 cores, 7.8 GB RAM, Python 3.13.16, numpy 2.5.3.

## Repository layout

```
yamanaka_replay/
  oracle.py        colony-formation oracle, calibrated to the 2006 figures
  inference.py     exact posterior over 10,626 hypotheses; expected info gain
  strategies.py    historical protocol, adaptive design, ablation
  cardinality.py   recipe size as a latent variable; GeneralOracle
experiments/
  01_main_comparison.py     the headline table
  02_prior_robustness.py    degrading the literature prior
  03_misspecification.py    one-at-a-time model error
  04_trace_single_run.py    one run, experiment by experiment
  _harness.py               per-run records, process pool, hypothesis rebuild
  05_scale_main.py          1,200 runs: intervals and paired tests
  06_synergy_threshold.py   1,100 runs: synergy strength vs calibration error
  07_stopping_rule.py       1,680 runs: stopping threshold and calibration
  08_scaling.py               300 runs: 24 to 60 candidates
  09_prior_quality.py       1,800 runs: prior strength x coverage x breadth
  10_unknown_cardinality.py   240 runs: fixed vs latent recipe size
  11_classical_baseline.py    classical threshold-group-testing comparator
  12_cardinality_baselines.py two-stage and model-free comparators
paper/           manuscript source and bibliography
data/
  prepare_norman.py  turns a published .h5ad into a small pseudobulk table
analysis/
  real_data_audit.py  the gap assumption against measured data
  figures.py       the manuscript figures
  stats.py         bootstrap, Wilson, paired permutation, exact McNemar, logistic
  test_stats.py    those five checked against independently known values
  analyze.py       the report; reads CSVs, runs no simulation
results/           per-run CSVs and analysis_report.txt
```

The statistics layer is hand-rolled to keep the dependency list at one package,
so it is tested: `python analysis/test_stats.py` checks each estimator against a
closed-form or planted-truth answer, including the permutation test's
false-positive rate under the null.

## Limitations

- **The oracle is a simulator.** Calibrated to published summary figures, not to
  raw data. It is a model of the assay, not the assay.
- **No real perturbation data is used.** The natural next validation is to
  repeat this against measured combinatorial perturbation datasets, where the
  first question is whether graded synergy of the assumed shape exists at all.
- **The historical protocol is reconstructed from the paper's description**, not
  from lab notebooks. The 61-experiment count follows the documented procedure;
  the real campaign involved further work not captured here.
- **Exact enumeration binds past ~60 candidates**, and the proposal generator
  degrades before that. Both need replacing for genome-scale candidate lists.
- **The latent-cardinality designer is bounded too**: it allows sizes 3–6, so a
  seven-factor truth would reproduce the original failure. The fix widens the
  assumption; it does not eliminate it.

## Citation

See [`CITATION.cff`](CITATION.cff).

## References

- Takahashi K. & Yamanaka S. (2006). *Induction of pluripotent stem cells from
  mouse embryonic and adult fibroblast cultures by defined factors.* Cell
  126(4):663–676. [doi:10.1016/j.cell.2006.07.024](https://doi.org/10.1016/j.cell.2006.07.024)
- Okita K., Ichisaka T. & Yamanaka S. (2007). *Generation of germline-competent
  induced pluripotent stem cells.* Nature 448:313–317.
  [doi:10.1038/nature05934](https://doi.org/10.1038/nature05934)

## License

MIT — see [`LICENSE`](LICENSE).

## Acknowledgement of tooling

The code in this repository was written with the assistance of Claude
(Anthropic), under the direction of the author, who designed the study,
specified the experiments, and is responsible for all claims and
interpretations.
