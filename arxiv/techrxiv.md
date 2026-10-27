# The interim preprint venue

arXiv is where this paper belongs and is being pursued (`endorsement.md`),
but the endorsement takes as long as it takes. This file is about where to
post in the meantime, after bioRxiv declined the paper at screening
(`biorxiv.md`).

## TechRxiv — ruled out, 22 October 2026

IEEE's preprint server was the obvious interim choice: free, no endorser,
DOI on posting, engineering and computer-science scope, four-day
moderation, appealable decline.

**Its submissions are closed.** The front page carries: "We are currently
preparing a transition to a new platform. Submissions are temporarily
closed during this process." The most recent posting on the site is dated
6 March 2026, so the freeze has run for roughly seven months with no
announced end. It is not an option today and cannot be planned around.

Source: https://www.techrxiv.org/

## What the Zenodo record does and does not do

The existing DOI `10.5281/zenodo.23436007` secures priority, but it is
typed **Software** and its only file is
`mashal-alkhalifah/yamanaka-replay-v1.3.0.zip`. The manuscript is inside
that archive, which means the paper is deposited but not *readable* — no
landing page shows an abstract, and no indexer treats it as a paper.

**A second Zenodo record, typed Preprint, holding `paper/main.pdf`
directly, fixes that in about two minutes.** It is free, needs no
endorsement or screening, issues its own DOI, and is harvested by OpenAIRE
and picked up by Google Scholar. Link it to the software record with a
*supplements* relation so the two read as one body of work.

This is the cheapest useful step and it conflicts with nothing later:
posting the same paper to arXiv afterwards is routine.

## Preprints.org — possible, with one caveat

Scope explicitly covers **Computer Science and Mathematics**. No
endorsement; an MDPI account is all that is required. CC BY 4.0 by
default, which matches the licence already chosen. A DOI is issued and
screening takes about 24 hours. Posting after submission to a conference is
allowed where the venue permits it, and RECOMB permits preprints without
restriction.

Two things to know before choosing it:

- It wants **source files**, Word or LaTeX, not a PDF alone. For LaTeX
  they ask for a zip containing everything needed to rebuild the PDF,
  including the `.bib`. That archive is prepared at
  `arxiv/preprints-org-source.zip`.
- Its instructions **recommend against posting the same paper to more than
  one preprint server**, on the grounds that metrics fragment and readers
  are confused by duplicates. Since arXiv is the destination, posting here
  first means doing the thing they advise against. It is a recommendation,
  not a prohibition, and nothing stops it — but it is a reason to prefer
  the Zenodo preprint record, which is a repository deposit rather than a
  competing preprint server.

Source: https://www.preprints.org/instructions-for-authors

## Recommendation

1. Post the **Zenodo preprint record** now. Free, immediate, no screening
   risk, no duplicate-server problem.
2. Pursue the **arXiv endorsement** in parallel; that is the posting that
   matters, both for this paper's audience and for RECOMB's requirement
   that accepted authors post to arXiv or bioRxiv.
3. Hold Preprints.org in reserve, for the case where no endorser is found
   by early December and RECOMB's notification is approaching.

## Fields, for whichever route is taken

| Field | Value |
|---|---|
| Title | Factor-combination discovery is stochastic threshold group testing: a counterfactual replay of the 2006 reprogramming screen |
| Resource type | Preprint |
| Author | Mashal Saleh Alkhalifah (sole author, corresponding) |
| Affiliation | Institute of Smart Cities, King Abdulaziz City for Science and Technology (KACST), Riyadh, Saudi Arabia |
| Email | mkhalifah@kacst.gov.sa |
| ORCID | 0009-0007-7058-0345 |
| Abstract | `arxiv/abstract.txt`, verbatim |
| File | `paper/main.pdf` — the full 20-page version |
| Licence | CC BY 4.0 |
| Funding | none |
| Conflict of interest | none declared |

Keywords, the eight lodged with RECOMB:

```
adaptive experimental design
Bayesian optimal experimental design
threshold group testing
combinatorial search
cellular reprogramming
model misspecification
latent variables
active learning
```

Related identifiers:

- `10.5281/zenodo.23436007` — *is supplemented by* (the software record)
- `https://github.com/mashal-alkhalifah/yamanaka-replay` — *is supplemented by*
- `10.1016/j.cell.2006.07.024` — *is derived from*
- `10.1126/science.aax4438` — *references*

## After posting

Add the new DOI to `CITATION.cff` (the `preferred-citation` block is
there, commented, with placeholders), to the README's citation section,
and to the RECOMB submission's comments field.
