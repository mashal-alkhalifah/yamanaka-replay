# arXiv submission details

Everything the submission form asks for. The upload itself is
`yamanaka-replay-arxiv.tar.gz` in the repository root, built from this folder.

## Title

Factor-combination discovery is stochastic threshold group testing: a
counterfactual replay of the 2006 reprogramming screen

## Authors

Mashal Saleh Alkhalifah

Affiliation: Institute of Smart Cities, King Abdulaziz City for Science and
Technology (KACST), Riyadh, Saudi Arabia

Sole author. No co-authors, no acknowledged collaborators beyond the tooling
note in the manuscript.

## Abstract

Paste `abstract.txt` verbatim. It is 1,909 characters, inside arXiv's 1,920
limit, and matches the manuscript's abstract apart from LaTeX markup.

## Categories

| | |
|---|---|
| Primary | `q-bio.QM` — Quantitative Methods |
| Cross-list | `stat.ME` — Methodology |
| Cross-list | `cs.LG` — Machine Learning |

The primary belongs in q-bio because the object of study is a biological
screen and the audit is against measured perturbation data; `stat.ME` carries
the Bayesian experimental-design content and `cs.LG` the active-learning
audience. `math.ST` would also be defensible but the paper proves nothing, so
it would attract the wrong referees.

**Endorsement**: confirmed required, and it is the one blocker left. A first
submission to `q-bio.QM` needs an endorsement from an arXiv author who has
submitted **2 papers** to any `q-bio.*` category between three months and five
years ago. Switching archive does not avoid it: `cs.LG` needs the same thing
and asks for **3 papers**, so `q-bio.QM` is both the better fit and the lower
bar. arXiv emails a unique endorsement code to the account address; that code
is account-specific and is deliberately not recorded in this public
repository.

To check whether a particular person can endorse: open one of their arXiv
papers and use the "Which of the authors of this article can endorse?" link at
the foot of the abstract page.

The submission itself is saved on arXiv and can be resumed once the
endorsement lands; nothing has to be re-entered.

## Comments field

20 pages, 10 figures. Code, per-run results and the full statistical report:
https://github.com/mashal-alkhalifah/yamanaka-replay

## License

Recommended: **CC BY 4.0**. It is the most permissive of arXiv's options, it
matches the repository's MIT licence in spirit, and it imposes no barrier to
later journal submission — arXiv preprints under CC BY are accepted by
Bioinformatics, PLOS Computational Biology, Nature Communications and Cell
Systems. The alternative, arXiv's own non-exclusive licence, is more
restrictive for readers and buys nothing here.

## Building the upload

```bash
cd arxiv && tar czf ../yamanaka-replay-arxiv.tar.gz main.tex main.bbl fig*.pdf
```

`main.bbl` is included deliberately: arXiv runs `pdflatex` but **not**
`bibtex`, so a submission without the `.bbl` loses every citation. The package
has been test-compiled in a clean directory with `pdflatex` alone, twice, and
produces the 20-page PDF with no overfull boxes and no unresolved references.

`\graphicspath` is set to `{./}` here rather than `{../figures/}` as in
`paper/main.tex`, because arXiv flattens the upload into one directory.

## After acceptance of the posting

The repository's `CITATION.cff` and the README's citation block both carry a
placeholder for the arXiv identifier. Fill both in once the identifier is
issued, and add the DOI arXiv mints.
