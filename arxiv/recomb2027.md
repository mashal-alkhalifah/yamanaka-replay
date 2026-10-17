# RECOMB 2027 submission

## Deadlines

| | |
|---|---|
| Abstract registration (mandatory) | **28 October 2026, 23:59 AoE** |
| Full paper | 4 November 2026, 23:59 AoE |
| Notification | 22 December 2026 |
| Conference | 17–20 May 2027, Toronto |

Submission system: EasyChair,
https://easychair.org/conferences/?conf=recomb2027

**The registered abstract is binding.** RECOMB states that placeholder
registrations are not accepted, that the abstract should closely match the
final paper in title, author list and abstract, and that substantial changes
to either may lead to removal without consideration. So the abstract below is
written as the abstract of the ten-page paper, not as a reservation of a slot.
After 28 October no new submissions are accepted, but the paper itself can be
revised until 4 November.

## Title

Factor-combination discovery is stochastic threshold group testing: a
counterfactual replay of the 2006 reprogramming screen

## Author

Mashal Saleh Alkhalifah, Institute of Smart Cities, King Abdulaziz City for
Science and Technology (KACST), Riyadh, Saudi Arabia.
mkhalifah@kacst.gov.sa

Sole author, so the four-submissions-per-author cap is not in play.

## Abstract to paste

See `recomb_abstract.txt`. It is the arXiv abstract trimmed to the five
results, since the ten-page paper carries all five in the main text with the
supporting campaigns in the appendix.

## What RECOMB requires that is already in place

- **Public code with a documented README by the submission deadline**, linked
  from the paper: https://github.com/mashal-alkhalifah/yamanaka-replay —
  done, and the link is already in the manuscript.
- **Preprints are allowed without restriction**, and accepted authors must
  post the final version to arXiv or bioRxiv anyway, so posting now is not
  only permitted but ahead of schedule.
- **Single-blind review**, so no anonymisation is needed.

## Built

- `paper/recomb10.tex` -> `recomb10.pdf`: **7 pages of body** plus a separate
  title page and 2 pages of references, against a 10-page limit that excludes
  both. Zero overfull boxes, no unresolved references. Four figures in the
  main text: the main comparison, cardinality, the measured gap, and the
  latent gap, plus the retrospective benchmark that the headroom allowed back
  in.
- `paper/appendix.tex` -> `appendix.pdf`: 10 pages, submitted as the single
  supplementary PDF. Carries the campaign table, the five displaced figures
  with full captions, and the complete statistical report verbatim.

The manuscript is 20 pages; RECOMB allows 10 (Letter, 10pt minimum, 1-inch
margins, bibliography excluded, separate title page allowed) plus an appendix
of unlimited length that reviewers may or may not read.

Proposed split, keeping one thesis in the main text and the supporting
campaigns behind it:

**Main text (10 pages)**
1. The screen as stochastic threshold group testing, and the historical
   protocol's $\Theta(n)$ cost.
2. The replay: 10 experiments against 61, 300 paired replicates, with the
   random-pool ablation.
3. Assumed cardinality as a symmetric silent failure, removed by latent size.
4. The measured audit: $\rho = 0.77$, not 0.10, across 171 combinations.
5. The designer run on that structure: 2.5% recovery at 0.991 confidence.
6. The gap response as a latent variable: 99.3%, and better than a pilot.
7. Discussion: every fixed assay parameter is a silent failure.

**Appendix**
Experiments 5–9 (the large campaigns and their statistics), 11 and 15 (the
classical baseline and the proposal-generator ablation), 14 (the retrospective
Norman benchmark), the stopping-rule analysis, and the full statistical
report.

The scaling and proposal-generator material is the clearest candidate for the
appendix: it is a self-contained methods result that does not carry the
paper's thesis, and it is the natural basis for a separate paper later.
