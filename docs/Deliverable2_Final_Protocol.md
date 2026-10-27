# Deliverable 2 — Final Research Protocol (locked 2026-10-26, before any survival estimate or HR was computed)

Disclosure: cell sizes and crude death counts per subgroup were tabulated during the feasibility audit. No Kaplan–Meier estimate, log-rank test, hazard ratio or phenotype comparison had been computed when this protocol was fixed.

## 1. Title (working)
DNMT3A and TET2 co-mutations and the prognostic difference between FLT3-ITD and FLT3-TKD in NPM1-mutated adult AML: a multi-cohort analysis of public clinico-genomic data (TCGA and Beat AML via cBioPortal; AMLSG).
The word "pathways" is dropped, because the study measures mutation status rather than pathway activity.

## 2. Questions and hypotheses
- **Primary:** among adults with NPM1-mutated AML and an assay- or sequence-verified FLT3-ITD-only or FLT3-TKD-only genotype, does DNMT3A and/or TET2 mutation status modify the association between FLT3 subtype and OS? Estimand: the ratio of hazard ratios (FLT3-ITD vs FLT3-TKD in epi-mutated) / (FLT3-ITD vs FLT3-TKD in epi-wild-type), from a cohort-stratified Cox model.
- **Secondary:**
  - (S1) DNMT3A, TET2 and composite mutation frequencies in ITD vs TKD.
  - (S2) Within NPM1/FLT3-TKD, is epi-mutation associated with baseline WBC and BM blast %? All four groups are also described.
  - (S3) Gene-specific (DNMT3A) modification.
- All hypotheses are two-sided and may be unsupported. "Risk equalization" is not an endpoint. No equivalence margin is specified, so no equivalence claim will be made.

## 3. Data sources
| Cohort | Source | ID | Specimen | ITD method | TKD method | DNMT3A/TET2 method |
|---|---|---|---|---|---|---|
| TCGA | cBioPortal | `laml_tcga_pub` | diagnostic | WGS/WES sequence-derived (validated classifier) | sequence | sequence (tumour–normal) |
| Beat AML | cBioPortal | `aml_ohsu_2022` | initial-diagnosis specimen | clinical PCR/CE `FLT3_ITD` | WES | WES |
| AMLSG | GitHub `gerstung-lab/AML-multistage/data` | AMLSG 07-04, AML-HD98A/B | diagnostic | clinical `FLT3_ITD` | clinical `FLT3_TKD` | 111-gene targeted panel |

## 4. Eligibility
- **Inclusion:**
  - age ≥18 at diagnosis;
  - AML (de novo, secondary or therapy-related);
  - pretreatment/diagnostic specimen (OHSU: `DISEASE_STAGE_AT_SPECIMEN_COLLECTION = Initial Diagnosis`, sequenced);
  - NPM1 mutated;
  - FLT3-ITD-only or FLT3-TKD-only;
  - OS time and status available.
- **One record per patient:**
  - OHSU: earliest diagnostic specimen; ties broken bone marrow > peripheral blood > leukapheresis, then sample ID.
  - TCGA and AMLSG: one sample each.
- **Excluded from primary analysis:**
  - ITD+TKD (described separately; sensitivity S6);
  - FLT3 non-ITD/non-TKD variants only;
  - FLT3-wild-type.
- Missing WBC/blasts do **not** exclude patients from survival analyses.

## 5. Variable definitions
- **FLT3-TKD (primary):** D835X or I836 alterations. **TKD-extended (sensitivity):** adds D839, N841, Y842 (cBioPortal cohorts only; AMLSG uses the clinical TKD call).
- **FLT3-ITD in TCGA (primary):** verified duplications plus probable JM insertions plus JM splice-region indels. **Strict definition (sensitivity):** verified duplications only.
- **DNMT3A/TET2 mutated (primary):** any protein-altering variant (missense, nonsense, frameshift, splice-site, in-frame indel). In AMLSG this means any called variant (`Result` ONCOGENIC or POSSIBLE).
- **Restricted definition (sensitivity):**
  - truncating variants, or
  - DNMT3A missense in the MTase domain (aa 634–912), or
  - TET2 missense in catalytic regions (aa 1104–1481 or 1843–2002).
  - In AMLSG, this is `Result = ONCOGENIC`.
- **Wild-type** is assigned only to sequenced samples covering the gene. All included samples were sequenced with assays covering DNMT3A and TET2.
- **OS:** months from diagnosis to death from any cause; survivors censored at last follow-up. AMLSG days are divided by 30.4375. Times of 0 are set to 0.01 months.
- **Covariates (a priori):** age (per 10 years) and sex.
- **Sensitivity covariate:** normal karyotype (harmonized from karyotype strings / the AMLSG `NK` flag).
- **Not used as covariates:** WBC and blasts (potential mediators), and ELN categories (contain NPM1/FLT3 information).

## 6. Analyses
- **A1 (co-mutation).** Per-cohort counts with denominators and Fisher exact OR (95% CI). Pooled Mantel–Haenszel OR stratified by cohort. BH correction across the 3 tests (DNMT3A, TET2, composite).
- **A2 (phenotype).**
  - Median (IQR) per group per cohort; Kruskal–Wallis across 4 groups per cohort, with Dunn post-hoc tests (BH).
  - Pooled linear models of log(WBC) and of BM blast % on FLT3 × epi with cohort fixed effects and HC3 SE. The key contrast is TKD epi-mut vs TKD epi-wt (geometric-mean ratio / mean difference).
- **A3 (OS).**
  - Kaplan–Meier by 4 groups per cohort: median OS (95% CI), 36-month OS, numbers at risk.
  - Log-rank tests stratified by cohort for the pre-specified pairwise contrasts: (G3 vs G1) ITD vs TKD in epi-wt; (G4 vs G2) ITD vs TKD in epi-mut; (G4 vs G3) epi effect within TKD.
- **A4 (primary model).**
  - Cox model: OS ~ ITD + EPI + ITD×EPI + age/10 + sex, strata(cohort). Report the interaction ratio of HRs (95% CI, Wald p) and stratum-specific HRs derived from model contrasts.
  - Proportional hazards: Schoenfeld residual tests.
  - Stability: events per parameter, plus 2,000-replicate cohort-stratified bootstrap percentile CI for the interaction.
  - Cohort-specific models with inverse-variance pooling of the interaction log-HR (fixed and DerSimonian–Laird; I²).
- **Sensitivity analyses.**
  - S1: restricted variant definition.
  - S2: DNMT3A-only modifier.
  - S3: adjustment for normal karyotype.
  - S4: leave-one-cohort-out.
  - S5: TCGA strict ITD / TKD-extended.
  - S6: ITD+TKD counted as ITD.
  - S7: bootstrap (see A4).
  - S8: EFS in AMLSG.
  - S9: de novo AML only.
- **Infeasible (pre-declared):** TET2-only modification, double-mutant subgroup, R882 vs non-R882 within TKD. These are reported descriptively only.
- Significance level α=0.05, two-sided. Findings are interpreted by estimates and CIs; non-significance is not read as equivalence.
- **VAF:** descriptive only; no clonal-order inference.

## 7. Software
The pre-specified tool was R. R `survival` could not be loaded in this sandbox (conda binary caused a Mingw-w64 relocation crash; CRAN binaries were blocked from loading from a writable path). Analyses are therefore run in Python 3.11:
- lifelines (KM, Cox, Schoenfeld tests);
- statsmodels (Mantel–Haenszel, OLS-HC3);
- scipy (Fisher, Kruskal–Wallis).

Firth-penalized Cox (coxphf) was unavailable and is replaced by bootstrap CIs.


## Addendum (2026-10-26, documentation correction; no change to analyses)
The implemented TKD-extended sensitivity definition (S5) includes all FLT3 activation-loop alterations at R834, D835, I836, M837, H838, D839, S840, N841 and Y842. §5 listed only D839/N841/Y842 in addition to D835/I836. Results were computed with the implemented (broader) definition.
