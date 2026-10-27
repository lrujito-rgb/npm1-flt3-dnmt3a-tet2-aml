# NPM1-mutated AML: does DNMT3A/TET2 co-mutation modify the FLT3-ITD vs FLT3-TKD survival difference?

Code, protocol and results for a multi-cohort public-data analysis (TCGA-LAML, Beat AML, AMLSG).
Manuscript in preparation; the article citation will be added here on publication.

Licence: code MIT, docs/results CC BY 4.0 (see `LICENSE`). An archived, citable release (Zenodo DOI) is linked from the repository page.

## Question and estimand
Among adults with NPM1-mutated AML and a single FLT3 alteration class (ITD-only or TKD-only),
is the overall-survival hazard ratio for FLT3-ITD vs FLT3-TKD different in DNMT3A/TET2-mutated vs
DNMT3A/TET2-wild-type disease? Primary estimand: the **ratio of hazard ratios** (ITD × DNMT3A/TET2
interaction) from a Cox model stratified by cohort and adjusted for age (per 10 y) and sex.
The analysis is observational and exploratory; non-significance is not interpreted as equivalence.

## Repository layout
```
run_all.py                  one-command reproduction (data -> provenance -> statistics -> figures)
requirements.txt            pinned Python packages (tested: Python 3.11.16, Windows 11)
environment.yml             conda alternative
src/build_data.py           download (cBioPortal REST API, AMLSG GitHub, UniProt) + harmonisation + FLT3 classifier
src/provenance.py           timestamps, API version, AMLSG commit, SHA-256 of every raw file
src/analyze.py              Tables 1-3, phenotype models, KM/log-rank, Cox models, S1-S9 sensitivity, bootstrap
src/figures.py              Figures 1-4 (style: src/publication.mplstyle)
docs/                       protocol, feasibility report, reproducibility audit, data dictionary, literature log
results/                    reference outputs (tables, results.json, Figure1-4 .png, data_provenance.json)
tests/check_reproduction.py compares a fresh run with tests/reference_results.json
data/README.md              data sources and terms (raw data are downloaded, not redistributed)
```

## Reproduce
```bash
git clone https://github.com/lrujito-rgb/npm1-flt3-dnmt3a-tet2-aml.git
cd npm1-flt3-dnmt3a-tet2-aml
python -m venv .venv && . .venv/bin/activate        # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python run_all.py
python tests/check_reproduction.py                  # prints REPRODUCED or the mismatching values
```
Internet access is required on the first run (cBioPortal API, raw.githubusercontent.com, rest.uniprot.org).
Runtime is about 6–8 min on a laptop, most of it the 2,000-replicate bootstrap (seed 20261026).
Outputs are written to `results/` (override with `AML_OUT`); raw downloads go to `data/raw/` (`AML_RAW`).

## Pipeline and documentation (in execution order)
| Stage | Document / code |
|---|---|
| 1. Feasibility audit (cohorts, FLT3 annotation, power for the interaction) | `docs/Deliverable1_Feasibility_Report.md` |
| 2. Final protocol (pre-specified before outcome analysis) incl. TKD-extended addendum | `docs/Deliverable2_Final_Protocol.md` |
| 3. Extraction and harmonisation; variable definitions | `src/build_data.py`, `docs/data_dictionary.csv`, `results/cohort_flow.csv` |
| 4–5. Statistical analysis | `src/analyze.py` |
| 6. Figures | `src/figures.py` |
| 7. Literature verification (PMID/DOI checked) and novelty assessment | `docs/references_verified.csv`, `docs/literature_novelty_report.md` |
| 8–9. Reproducibility audit (what was done, deviations, principal results) | `docs/Deliverable8_Reproducibility_Publication_Audit.md` |

## Cohorts
| Cohort | Source | Primary set n | Deaths | Median follow-up (reverse KM, months) |
|---|---|---|---|---|
| TCGA-LAML | cBioPortal `laml_tcga_pub` | 29 | 20 | 41.4 |
| Beat AML | cBioPortal `aml_ohsu_2022` | 76 | 48 | 30.9 |
| AMLSG 07-04 / HD98A / HD98B | `gerstung-lab/AML-multistage` @ `a3f96ba` | 205 | 128 | 79.7 |

Cohorts are never pooled naively: all pooled models are stratified by cohort; cohort-specific estimates
and inverse-variance pooling are reported alongside. Duplicate cBioPortal versions of TCGA and the Beat AML
2018 subset are not combined with these studies.

## Data provenance of this release
Source data were first downloaded on 26 October 2026 and re-downloaded for a clean end-to-end reproduction on
2026-10-27 (cBioPortal v7.2.0, database schema 3.1.0;
AMLSG commit `a3f96ba727bfa3f30e47d9419a9faf1caab44d51`, dated 2019-03-25). The clean run reproduced every output table byte-for-byte.
SHA-256 checksums of all 13 raw files are in `results/data_provenance.json`. The AMLSG commit is pinned by default;
override with the environment variable `AMLSG_REF` only to test a different data version.

## Reference results (this release)
Primary Cox model, n = 310, 196 deaths; events per parameter 39.2; minimum Schoenfeld p = 0.108.

| Contrast | HR (95% CI) | p |
|---|---|---|
| HR ITD vs TKD, DNMT3A/TET2 wt | 1.67 (0.89–3.11) | 0.110 |
| HR ITD vs TKD, DNMT3A/TET2 mut | 2.16 (1.30–3.61) | 0.003 |
| HR DNMT3A/TET2 mut vs wt, within TKD | 1.21 (0.59–2.48) | 0.610 |
| HR DNMT3A/TET2 mut vs wt, within ITD | 1.57 (1.11–2.20) | 0.010 |
| Ratio of HRs (interaction) | 1.30 (0.59–2.87) | 0.518 |

Bootstrap (B = 2000, 2000 successful) ratio of HRs: median 1.29, percentile 95% CI 0.52–2.80.
The interaction estimate is imprecise; the data neither establish nor exclude effect modification.

## Deviations from the protocol (also reported in the article)
1. Analysis language: the protocol specified R (`survival`, `coxphf`); R could not be loaded in the analysis environment, so all analyses were run in Python (lifelines/statsmodels).
2. Firth-penalised Cox was replaced by a 2,000-replicate cohort-stratified bootstrap as the small-cell stability check.

Note on FLT3 annotation (pre-specified, not a deviation): FLT3-ITD is taken from the clinical PCR/CE field in Beat AML and
from clinical annotation in AMLSG; in TCGA it is called from sequence by the classifier in `src/build_data.py`,
which was validated against clinical PCR in `aml_ohsu_2018` (Deliverable 1). Beat AML sequence ITD calls are not used
because cBioPortal stores them as an identical placeholder record.


## How to cite
Cite the article and this release (see `CITATION.cff`). Also cite the data sources:
TCGA (NEJM 2013; PMID 23634996), Beat AML (Cancer Cell 2022; PMID 35868306), AMLSG (NEJM 2016, PMID 27276561; Nat Genet 2017, PMID 28092685), and cBioPortal.
