# yamanaka-replay

**Replaying the 2006 discovery of the Yamanaka factors in silico, under its original fitness function, to ask a counterfactual question: how many wet-lab experiments would adaptive experimental design have needed?**

Answer, on a simulator calibrated to the published figures: **a median of 6–11 experiments, against 61 for the historical leave-one-out protocol.**

> **This is a methodological study on a simulator, not a biological result.**
> The oracle is a calibrated model of a selection assay, not an assay. Nothing
> here supports any claim about cellular biology.

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
with three of the four produce a weak but nonzero signal, which is the only
gradient the search can climb. **A pool-size penalty**: adding factors lowers
per-cell efficiency, which is why the 10-factor subset beat the full 24 — a
counter-intuitive result the simulator reproduces rather than assumes.

## Results

All figures from `./reproduce.sh`; raw CSVs in [`results/`](results/).

### Main comparison — 20 runs, fresh hidden recipe and noise seed each

| Strategy | Median experiments | Mean | Success | Saving |
|---|---|---|---|---|
| Historical protocol (2006) | 61 | 61.0 | 100% | — |
| Adaptive design (information gain) | **11** | 10.3 | 100% | **82%** |
| Adaptive + literature prior | **6** | 6.3 | 100% | **90%** |
| Random pools + Bayes *(ablation)* | 49 | 47.7 | 80% | 20% |

**The ablation is the load-bearing control.** Identical Bayesian machinery with
randomly chosen pools needs 49 experiments and succeeds only 80% of the time.
The saving comes from *informative experiment selection*, not from Bayesian
updating.

### Robustness to an incomplete literature prior

Published candidate lists miss genes, so the prior was degraded by dropping
true factors out of it:

| Prior quality | Median experiments | Success |
|---|---|---|
| Contains all four true factors | 6 | 100% |
| Drops one true factor | 8 | 100% |
| Drops two true factors | 10 | 100% |
| Drops three true factors | 12 | 93% |

Degradation is graceful, because the prior is a **weighting, not a
constraint** — data overturns it within a few rounds. Even a badly misleading
prior stays roughly five times cheaper than the historical protocol.

### Where the method breaks

The designer assumes a likelihood model. Making the oracle disagree with it,
one assumption at a time:

| Oracle reality vs designer's assumption | Kind | Median | Success |
|---|---|---|---|
| Matches | calibration | 11 | 100% |
| Pool-size penalty 3× harsher | calibration | 12 | 100% |
| Triplet effect stronger than assumed | calibration | 13 | 79% |
| Overall efficiency 3× lower | calibration | 18 | 43% |
| **Triplets have no effect at all** | **structural** | 25 | **43%** |
| **Truth is five factors, not four** | **structural** | 26 | **0%** |

Calibration error is survivable. Two structural errors are not:

1. **No graded synergy.** If the phenotype is strictly all-or-nothing at four
   factors, there is no gradient and the search degenerates toward enumeration.
2. **Wrong recipe cardinality.** Searching a fixed size of four never finds a
   five-factor truth — and the failure is **silent**: the posterior still
   concentrates confidently on a wrong answer. Recipe size must therefore be an
   inferred latent variable in any real deployment, never a fixed
   hyperparameter. This is the single most actionable finding here.

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
three-factor signal at experiments 4–9 is the thread the search pulls on.
And the designer never once tests a single factor, where the historical
protocol spent 24 of its 61 experiments.

## Related work, and what is actually new here

Adaptive and Bayesian experimental design for biological perturbation
selection is **an established field**, and this repository does not claim to
introduce it:

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
   retrospective studies replay adaptive data collection over existing
   datasets and report predictive accuracy. This replays the path to a single
   known discovery under its own original assay, and compares against that
   discovery's *actual, published* experiment count.
2. **A reproducible benchmark.** 24 candidates, 10,626 hypotheses, an oracle
   pinned to published numbers, four strategies including an ablation, one
   dependency, under a minute on a laptop.
3. **Characterised structural failure modes**, in particular the silent
   zero-success failure under misspecified recipe cardinality, which to our
   knowledge is not reported elsewhere and has direct implications for
   deployment.

Claim 3 is the part most likely to be useful to others. Claims 1 and 2 are
infrastructure.

## Reproducing

```bash
git clone https://github.com/OWNER/yamanaka-replay
cd yamanaka-replay
pip install -r requirements.txt
./reproduce.sh
```

Runs in about two minutes on two CPU cores. **No GPU, no network access, no
external data.** Every number in this README is produced by those four
scripts and nothing else.

Reference environment: Linux 6.18, x86_64, 2 cores, 7.8 GB RAM, Python 3.13.16,
numpy 2.5.3.

## Repository layout

```
yamanaka_replay/
  oracle.py        colony-formation oracle, calibrated to the 2006 figures
  inference.py     exact posterior over 10,626 hypotheses; expected info gain
  strategies.py    historical protocol, adaptive design, ablation
experiments/
  01_main_comparison.py    the headline table
  02_prior_robustness.py   degrading the literature prior
  03_misspecification.py   where the method breaks
  04_trace_single_run.py   one run, experiment by experiment
results/          CSV outputs
```

An implementation note worth flagging: the colony-count distribution depends on
a hypothesis only through `|pool ∩ hypothesis|`, which takes five values. Both
the likelihood and the expected information gain therefore collapse to a
five-way mixture, so the cost of choosing an experiment is independent of how
many hypotheses remain alive. The design step scales to much larger candidate
sets than 24.

## Limitations

- **The oracle is a simulator.** Calibrated to published summary figures, not
  to raw data. It is a model of the assay, not the assay.
- **The hypothesis space assumes a fixed recipe size of four.** Experiment 3
  shows this is the method's most dangerous assumption; lifting it is the
  obvious next step.
- **The historical protocol is reconstructed from the paper's description**, not
  from lab notebooks. The 61-experiment count follows the documented procedure;
  the real campaign involved further work not captured here.
- **No real perturbation data is used.** The natural next validation is to
  repeat this against measured combinatorial perturbation datasets.

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
