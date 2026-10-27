# Deliverable 8: Reproducibility and Publication-Readiness Audit (2026-10-26)

## A. What was actually done (verified)

| Phase | Status | Evidence |
|---|---|---|
| 1 Feasibility audit | Done | Deliverable 1 (cBioPortal API inventory of 10 AML studies; FLT3 classifier validated against clinical PCR, sensitivity 0.82 / specificity 0.975) |
| 2 Protocol locked before outcome analysis | Done | Deliverable 2 (crude cell counts were known at locking; disclosed) |
| 3 Extraction and harmonisation | Done | `build_data.py` → `analytic_dataset.csv`, `cohort_flow.csv`, `data_dictionary.csv` |
| 4 Statistics | Done, all pre-specified analyses except Firth | `analyze.py` → Tables 1–3, phenotype, Schoenfeld, sensitivity S1–S9, `results.json` |
| 5 Code | Done, executed end to end | `aml_pipeline/` (build_data → analyze → figures → build_docx) |
| 6 Figures | Done | Figures 1–4 (render-checked) |
| 7 Literature | Done (PubMed, Europe PMC, Crossref); supplements not read | 72 verified references (65 cited); novelty report |
| 8 Manuscript | Done | Revised_Manuscript_NPM1_FLT3_DNMT3A_TET2.docx |

## B. Deviations from the master prompt and protocol (must be disclosed)

1. **Data source scope.** AMLSG (non-cBioPortal) was added, as decided by the investigator, because cBioPortal alone gave 22 NPM1/TKD patients, with only 4 deaths in the TKD/epi-wild-type cell. The title was changed accordingly.
2. **Software.** The analysis was run in Python (lifelines 0.30.3, statsmodels 0.15.0) instead of R.
   - The R `survival` package crashed in this sandbox (Mingw-w64 relocation error after conda install; CRAN binaries are blocked from writable paths).
   - **Side effect:** the shared base R environment in this workspace now fails at kernel start-up and cannot be rolled back from here (additive-only policy). Recreate or repair it via Customize → Environments before using R again.
3. **Firth-penalised Cox** (coxphf) was not run, because the package was unavailable. It was replaced by a 2,000-replicate bootstrap. Every FLT3 × epi cell had ≥13 events in the pooled model.
4. **TCGA "de novo" flag** was taken from the source publication (all de novo), not from a data field.
5. **Literature sub-agent:** supplementary tables of Papaemmanuil 2016, HARMONY 2026, Cocciardi 2025 and Ni 2025 were not read. Bezerra 2020 and Boddu 2017 were read as abstract or title only. The novelty verdict is therefore provisional.

## C. Principal results (copied from output tables)
- **Cohort:** n=310 (TCGA 29, Beat AML 76, AMLSG 205); 196 deaths; 244 ITD and 66 TKD.
- **ITD vs TKD hazard ratio:**
  - epi-wild-type: 1.67 (0.89–3.11);
  - epi-mutated: 2.16 (1.30–3.61).
- **Interaction (ratio of HRs):** 1.30 (0.59–2.87), P=0.52; bootstrap CI 0.52–2.80; I²=0%. Across sensitivity analyses the estimate ranged from 0.75 to 2.05, and every CI included 1.
- **Epi-mutation effect by FLT3 subtype:**
  - within ITD: HR 1.57 (1.11–2.20);
  - within TKD: HR 1.21 (0.59–2.48).
- **WBC:** epi-mutation was associated with higher WBC in TKD (GMR 1.77, 1.02–3.09) and in ITD (GMR 1.64, 1.15–2.35), with no interaction. This is exploratory and not multiplicity-adjusted.
- **Co-mutation frequency:** no clear ITD vs TKD difference in DNMT3A/TET2 frequency (MH OR 1.46, 0.84–2.54).

_Sections on submission planning are kept with the authors and are not part of this public record._
