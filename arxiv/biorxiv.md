# bioRxiv submission

The preprint route that needs no endorser. arXiv requires an endorsement from
an author with prior submissions in the archive, which this account does not
have and cannot quickly obtain; bioRxiv has no such system. RECOMB's call
accepts either — accepted authors must post to "bioRxiv or arXiv" — so nothing
is given up by posting here instead.

It is also the better-targeted venue for this particular paper. The audit in
Section 3.4 is against Norman et al. and Joung et al., and the submission
target is a computational biology conference; that readership is on bioRxiv.
arXiv's advantage is the machine-learning and statistics audience, which is
this paper's secondary one.

## What bioRxiv requires (checked against their FAQ)

- **No endorsement and no sponsor.** No such system exists there.
- **No fee.**
- **Computational papers are explicitly welcome**, provided they are complete
  research articles with methodological detail and results — which this is.
  What they turn away is simple reanalysis of public data or a software
  announcement without methods.
- **A single PDF is accepted.** LaTeX must be converted first, which it is.
- **CC BY is among the licence options**, matching the choice already made.
- Screening is for plagiarism, scope and article type, not peer review, and
  articles are explicitly not endorsed by bioRxiv.
- **Registration is required before depositing**, with an affiliation
  supplied, so that there is an institution to contact if misconduct concerns
  are raised. Registration is free.

## Fields

| Field | Value |
|---|---|
| Title | Factor-combination discovery is stochastic threshold group testing: a counterfactual replay of the 2006 reprogramming screen |
| Author | Mashal Saleh Alkhalifah (sole author, corresponding) |
| Affiliation | Institute of Smart Cities, King Abdulaziz City for Science and Technology (KACST), Riyadh, Saudi Arabia |
| Email | mkhalifah@kacst.gov.sa |
| Abstract | `arxiv/abstract.txt`, verbatim |
| Subject area | Bioinformatics (alternative: Systems Biology) |
| Licence | CC BY 4.0 |
| Manuscript file | `paper/main.pdf` — the full 20-page version |
| Article category | **New Results** |
| Supplementary | none needed; the appendix belongs to the RECOMB submission |

Article category: bioRxiv asks for New Results, Confirmatory Results or
Contradictory Results. **New Results** is correct — the paper advances a
field rather than replicating or contradicting a specific prior study. The
measured audit refutes an assumption of our own earlier work, not of someone
else's, so Contradictory Results would misdescribe it.

Subject area: **Bioinformatics** fits the method and the audit better than
Systems Biology, which skews toward models of pathways and networks rather
than experimental design over them.

## Which PDF

Post the **full 20-page `main.pdf`**, not the 10-page RECOMB version. The
preprint has no page limit, the longer version carries all sixteen
experiments, and a reader arriving from the conference paper should find more
here rather than the same thing again.

## After posting

bioRxiv issues a DOI. Fill it into `CITATION.cff` (the `preferred-citation`
block is already there, commented, with placeholders) and into the README's
citation section, and add it to the RECOMB submission's comments field.

## Still worth doing, and free

## The order that costs the least

Three accounts are needed and each has to be created by hand. Doing them in
this order means creating **one** password rather than three:

1. **ORCID** (orcid.org) — a new account, about two minutes. An author
   identifier that gates nothing and keeps "Mashal Saleh Alkhalifah" attached
   to the right person across bioRxiv, Zenodo, RECOMB and every later paper.
2. **Zenodo** (zenodo.org) — offers *Sign in with ORCID*, so no new password.
   Then Settings -> GitHub, and switch on `yamanaka-replay`. The switch has to
   be on **before** a release is cut: Zenodo archives releases made after it,
   not before. `.zenodo.json` in the repository root already carries the
   metadata that record will use.
3. **bioRxiv** (biorxiv.org) — registration with an affiliation, free.
