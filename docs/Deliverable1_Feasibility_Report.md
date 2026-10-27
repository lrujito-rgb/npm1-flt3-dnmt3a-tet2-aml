# Deliverable 1 — Phase 1 Scientific and Data-Feasibility Audit

**Project:** Impact of DNMT3A/TET2 co-mutations on the prognostic divergence of FLT3-ITD vs FLT3-TKD in NPM1-mutated AML
**Data extraction date:** 2026-10-26 (cBioPortal public REST API, https://www.cbioportal.org/api; AMLSG data from github.com/gerstung-lab/AML-multistage, `data/`)
**Status labels used:** [VERIFIED] = checked directly in downloaded data; [LITERATURE] = from a published source; [NOT VERIFIED] = could not be established.

---

## 0. Source-document review

**Central question shared by both documents.** Whether DNMT3A and/or TET2 co-mutations modify the OS difference between FLT3-ITD and FLT3-TKD in NPM1-mutated adult AML (effect modification / interaction).

**Material that can be retained.**
- Rationale for distinguishing FLT3-ITD vs FLT3-TKD within NPM1-mutated AML.
- The four-group operational design (ITD/TKD × DNMT3A-TET2 wt/mut).
- Non-parametric phenotype comparisons (WBC, BM blasts). KM and log-rank as descriptive tools. Cox with an interaction term.
- Citation of primary dataset publications alongside cBioPortal.

**Assumptions that failed verification or require correction.**
1. "Combined TCGA + Beat AML give hundreds of profiles" → after eligibility only **22 NPM1/FLT3-TKD-only patients** with OS are available across both cBioPortal sources (Section 4). [VERIFIED]
2. TCGA PanCancer Atlas has WBC/blast data → **false**; WBC, BM/PB blasts and cytogenetic risk exist only in `laml_tcga_pub` (NEJM 2013). Both are the same 200 patients and must never be pooled. [VERIFIED]
3. Beat AML 2018 and Beat AML 2022 are distinct → **false**; `aml_ohsu_2022` carries `COHORT = Waves1+2` for the 2018 cohort (649 samples). OHSU 2022 is a superset. [VERIFIED]
4. "FLT3-ITD defined by in-frame insertions in exons 14–15 via genomic coordinates" → partially unsafe. In OHSU 2022, 137 ITD records are an identical placeholder (`E611_F612insLDFREYEYDLKWEFPRENLE`, all at chr13:28608279) that mirrors the clinical PCR call. They are not true sequenced alleles, so their VAF and length cannot be interpreted. [VERIFIED]
5. Beat AML "independent external validation" of TCGA, and the TARGET pediatric "validation" → inappropriate. TARGET is pediatric (out of scope). [VERIFIED]
6. Draft states ELN classifies NPM1/FLT3-TKD as favourable "because of TKD" and recommends up-staging → ELN 2022 classifies NPM1-mutated AML **without FLT3-ITD** as favourable (TKD does not modify risk). NPM1 with FLT3-ITD is intermediate regardless of allelic ratio. Adverse cytogenetics override both. [LITERATURE: Döhner et al., Blood 2022;140:1345, doi:10.1182/blood.2022016867]

**Claims in the draft that must not be retained unless supported by results:** "risk equalization"; "neutralizes the favourable prognosis"; "phenotypic convergence driven by DNMT3A/TET2"; STAT5-bypass, HOXA9/MEIS1 chromatin priming, and TKD "minimal STAT5" mechanisms (not measured); HMA/menin-inhibitor combinations as indicated therapy; recommendation to up-stratify ELN risk or counsel transplant; inference of clone order from VAF; GraphPad Prism 10 / R 4.3.1 (not used).

---

## 1. Candidate datasets [VERIFIED unless noted]

| Study ID (cBioPortal) | Official name | Source publication (PMID) | Samples / patients | Adult? | Decision |
|---|---|---|---|---|---|
| `laml_tcga_pub` | Acute Myeloid Leukemia (TCGA, NEJM 2013) | 23634996 | 200 / 200 | yes (18–88 y) | **Use** (TCGA source; clinical variables present) |
| `laml_tcga_pan_can_atlas_2018` | Acute Myeloid Leukemia (TCGA, PanCancer Atlas) | 29625048 et al. | 200 / 200 | yes | Exclude — same patients as above; no WBC/blast/cytogenetic risk |
| `laml_tcga`, `aml_tcga_gdc` | TCGA Firehose Legacy / GDC 2025 | — | 200 / 200 | yes | Exclude — same patients |
| `aml_ohsu_2022` | Acute Myeloid Leukemia (OHSU, Cancer Cell 2022) | 35868306 | 942 / 805 | mostly | **Use** (Beat AML source; superset of 2018) |
| `aml_ohsu_2018` | Acute Myeloid Leukemia (OHSU, Nature 2018) | 30333627 | 672 / 562 | mostly | Exclude from analysis (subset of 2022). **Used only to validate** the sequence-based ITD classifier |
| `heme_msk_impact_2022` | MSK-IMPACT Heme Tumors (MSK, 2022) | — | 2383 / 1937 | unknown (no age field) | Exclude — no age, OS from sequencing date, only 4 NPM1/TKD patients |
| `mnm_washu_2016` | AML or MDS (WashU, 2016) | 27959731 | 136 / 116 | yes | Exclude — decitabine trial, 2 FLT3 records |
| `aml_target_2018_pub`, `aml_target_gdc`, `aml_stjude_2024` | Pediatric AML | — | — | no | Exclude — pediatric |
| **Non-cBioPortal:** AMLSG (Papaemmanuil et al., NEJM 2016; Gerstung et al., Nat Genet 2017) | AMLSG 07-04, AML-HD98A, AML-HD98B trials | 27276561; 28092685 [LITERATURE — PMIDs confirmed in Phase 7 against PubMed/Crossref; see references_verified.csv] | 1540 / 1540 | yes (18–84 y) | **Candidate additional independent cohort** |

Population overlap: TCGA (Washington University, St Louis), Beat AML (OHSU-led US multi-centre), and AMLSG (German–Austrian trials) come from different institutions and periods. No overlap is expected, but this was not verifiable at the patient level. Within each source, all TCGA versions are identical patients, and OHSU 2018 is contained in OHSU 2022.

## 2. Variable availability (actual field names)

| Variable | `laml_tcga_pub` | `aml_ohsu_2022` | AMLSG |
|---|---|---|---|
| NPM1 | mutation file (WGS/WES; 53 frameshift ins) | clinical `NPM1` (PCR/CE consensus); 100% concordant with WES in eligible-stage samples (150/150 pos, 369/369 neg) | clinical `NPM1` |
| FLT3-ITD | **sequence-derived only** (no clinical ITD field) | clinical `FLT3_ITD` (PCR/CE consensus) + `ALLELIC_RATIO` | clinical `FLT3_ITD` + `AMLSG_FLT3ITD.txt` |
| FLT3-TKD | mutation file (D835/I836) | WES mutation file (D835/I836); `VARIANT_SUMMARY` (clinical text) supplementary | clinical `FLT3_TKD` |
| DNMT3A / TET2 | mutation file (tumour–normal, `Somatic`) | WES mutation file | `AMLSG_Genetic.txt` (111-gene panel; `Result` = ONCOGENIC/POSSIBLE) |
| VAF | alt/ref counts | alt/ref counts (FLT3-ITD placeholder VAF **not interpretable**) | `VAF` column |
| OS | `OS_MONTHS`, `OS_STATUS` (from diagnosis) | `OS_MONTHS`, `OS_STATUS` (from initial diagnosis; patient-level) | `OS` (days), `Status` |
| EFS/RFS | `DFS_*` (definition limited) | not available | `efs`/`EFSSTAT`, `rfs`/`rfsstat` |
| Age, sex | `AGE`, `SEX` | `AGE_AT_DIAGNOSIS` (sample), `SEX` | `AOD`, `gender` |
| WBC | `WBC` | generic assay `BLOOD_MEASUREMENTS: WBC_COUNT` | `wbc` |
| BM / PB blasts | `BM_BLAST_PERCENTAGE`, `PB_BLAST_PERCENTAGE` | `PERC_BLASTS_IN_BM`, `PERC_BLASTS_IN_PB` | `BM_Blasts`, `PB_Blasts` |
| Cytogenetics | `RISK_CYTO`, `CYTOGENETICS` | `KARYOTYPE` (ISCN string; risk must be derived) | `C_Risk` + karyotype indicators |
| Treatment | `INDUCTION`, `TRANSPLANT_TYPE` | `TYPE_INDUCTION_TX`, treatment types (incomplete) | trial arm, `TPL_o`, `TPL_type` |
| Specimen timing | primary diagnostic samples | `DISEASE_STAGE_AT_SPECIMEN_COLLECTION`, `TIME_OF_SAMPLE_COLLECTION_RELATIVE_TO_INCLUSION` | diagnostic samples |

Completeness among eligible NPM1/FLT3 (ITD-only or TKD-only) patients with OS:
- TCGA (n=29): WBC 29, BM blasts 29, PB blasts 29, cytogenetic risk 28.
- OHSU 2022 (n=76): WBC 61, BM blasts 64, PB blasts 53, karyotype 64, ITD allelic ratio 60, induction type 63.
- AMLSG (n=205): WBC 204, BM blasts 196, PB blasts 196, cytogenetic risk 186, HSCT 205.

## 3. Reliability of FLT3-ITD/TKD classification

**Pre-specified classifier for sequence-only data (TCGA).**
- ITD = insertion, duplication, or delins starting at codons 570–650 whose inserted residues contain at least 5 consecutive amino acids matching native FLT3 sequence (UniProt P36888, residues 560–650), or an explicit `dup` annotation ("verified").
- JM insertions without demonstrable duplication (e.g., inserts containing stop codons, consistent with intron-14-spanning ITDs), plus JM splice-region indels, are "probable" ITDs.
- TKD = D835X or I836 (del/sub/ins). Other activation-loop variants (D839G, N841X, Y842C) are a sensitivity-analysis category. All other FLT3 missense variants (e.g., A680V, N676K, V491L, S451F) are non-classifiable and are neither TKD nor ITD.

**Validation against clinical PCR/CE (Beat AML 2018; 622 sequenced samples with a clinical call; no placeholder present).**

| Sequence definition | TP | FN | FP | TN | Sensitivity | Specificity | PPV |
|---|---|---|---|---|---|---|---|
| verified only | 111 | 39 | 10 | 462 | 0.74 | 0.979 | 0.917 |
| verified + probable | 123 | 27 | 12 | 460 | 0.82 | 0.975 | 0.911 |

Interpretation: sequence-derived ITD calls in TCGA are acceptable but imperfect. About 18% of true ITDs may be missed, so some TCGA "TKD-only" patients could carry an undetected ITD (misclassification toward the null). This is a declared limitation, not fatal. OHSU 2022 and AMLSG use assay-based ITD calls and are the more reliable sources.

**Dual ITD+TKD:** OHSU 2022 = 1, TCGA = 0, AMLSG = 9. These are excluded from the primary contrast and described separately.

## 4. Eligible patients and deaths per primary subgroup

Primary epigenetic status = DNMT3A and/or TET2 protein-altering variant (AMLSG: `Result = ONCOGENIC`).

**Cohort flow.**
- OHSU 2022: 942 samples / 805 patients → AML 905/771 → initial-diagnosis specimen 574/556 → sequenced 551/538 → age ≥18 532/519 → one sample per patient 519 → NPM1-mutated 150 → FLT3 status: ITD-only 65, TKD-only 13, ITD+TKD 1, no ITD/TKD 71 → ITD-only or TKD-only with OS: **76**.
- TCGA: 200 → age ≥18 200 → NPM1-mutated 54 → ITD-only 19, TKD-only 10, none 25 → with OS: **29**.
- AMLSG: 1540 → NPM1-mutated 436 → ITD-only 161, TKD-only 44, ITD+TKD 9, none 222 → with OS: **205**.

| Group | TCGA n (deaths) | OHSU 2022 n (deaths) | cBioPortal total | AMLSG n (deaths) |
|---|---|---|---|---|
| G1 NPM1/ITD/epi-wt | 8 (5) | 25 (13) | 33 (18) | 61 (35) |
| G2 NPM1/ITD/epi-mut | 11 (8) | 39 (28) | 50 (36) | 100 (76) |
| G3 NPM1/TKD/epi-wt | 3 (2) | 5 (2) | **8 (4)** | 22 (9) |
| G4 NPM1/TKD/epi-mut | 7 (5) | 7 (5) | **14 (10)** | 22 (8) |

Gene-specific counts in the TKD stratum (cBioPortal): DNMT3A 12/22, TET2 3/22, both 1. Gene-specific analyses in the TKD stratum are therefore **infeasible**, and TET2-specific analysis in TKD patients is infeasible in all sources.

Sensitivity check (OHSU, any AML specimen rather than initial-diagnosis only): TKD-only 18 (12 deaths), ITD-only 83 (57 deaths). Even this relaxed definition (with left-truncation risk) does not resolve the sparsity.

**Approximate precision of the FLT3 × epigenetic interaction hazard ratio** (Wald SE ≈ √Σ1/dᵢ over the four cells):

| Data used | SE(log HR) | 95% CI spans × / ÷ | Min. detectable ratio of HRs (80% power, α=0.05) |
|---|---|---|---|
| TCGA only | 1.01 | 7.3 | 17.0 |
| OHSU 2022 only | 0.90 | 5.9 | 12.5 |
| cBioPortal pooled | 0.66 | 3.6 | 6.3 |
| AMLSG only | 0.53 | 2.8 | 4.4 |
| cBioPortal + AMLSG | 0.40 | 2.2 | 3.1 |

## 5. Main sources of missingness and bias
1. Sparse key cell (NPM1/TKD/epi-wt): 8 patients and 4 deaths in cBioPortal data.
2. ITD under-ascertainment in TCGA (sequence-only; sensitivity ≈ 0.82).
3. OHSU specimen timing: 35% of AML specimens were not taken at initial diagnosis. Restricting to diagnostic specimens avoids left truncation but reduces n.
4. Treatment heterogeneity: TCGA (2001–2010, intensive plus some decitabine/lenalidomide), OHSU (mixed intensive and non-intensive, era with midostaurin), AMLSG (intensive trials; no FLT3 inhibitor). FLT3-inhibitor exposure is not reliably captured.
5. Confounding by age: DNMT3A/TET2 mutations are associated with older age.
6. WBC/blasts are missing for about 20% of OHSU eligible patients. PB blasts are least complete.
7. No mutation-pathogenicity annotation (OncoKB) is available offline for cBioPortal sources; protein-altering variants are used, with a domain/truncating restriction as sensitivity analysis. AMLSG provides curated ONCOGENIC/POSSIBLE calls.

## 6. Recommendation: **PROCEED WITH MODIFICATIONS**

Criteria:
- (a) All essential variables (NPM1, assay-classifiable ITD vs TKD, DNMT3A, TET2, OS, age, sex) exist and are verifiable → the project is executable.
- (b) Using cBioPortal sources only, the primary interaction cannot be estimated with useful precision. Only a ≥6-fold difference between HRs would be detectable, and the reference cell has 4 deaths. A cBioPortal-only analysis must be framed as **descriptive and hypothesis-generating**.
- (c) Adding the independent, publicly released AMLSG cohort roughly halves the interaction SE (0.66 → 0.40). With all three cohorts it remains underpowered for moderate effects (minimum detectable ratio ≈ 3), so conclusions must stay estimation-focused.
- (d) Required design changes:
  - Never pool TCGA versions or OHSU 2018 with 2022.
  - Use cohort-stratified Cox models (baseline hazard by cohort) plus cohort-specific estimates.
  - Exclude ITD+TKD patients from the primary contrast.
  - Make the composite DNMT3A/TET2 exposure primary; DNMT3A-only analyses are secondary, and TET2-only and double-mutant analyses are infeasible.
  - Do not adjust the primary model for WBC/blasts.
  - Do not use ELN 2017/2022 categories as covariates (they encode FLT3-ITD/NPM1).
  - Retitle the study away from "pathways" and away from "cBioPortal" alone if AMLSG is used.

**Novelty note (preliminary; full check in Phase 7).** The NPM1/FLT3-ITD/DNMT3A triple-mutant has been studied extensively [LITERATURE: e.g., Bezerra et al., Blood 2020;135:870]. FLT3-TKD/NPM1 has been reported as a favourable group [LITERATURE: Boddu et al., Blood Adv 2017]. A recent HARMONY study of NPM1-mutated AML co-mutation patterns [LITERATURE: PMC12875867] must be checked for an overlapping TKD × DNMT3A analysis before any novelty claim. AMLSG data were used in the original Papaemmanuil/Gerstung gene–gene interaction models, so any AMLSG-based result is a secondary analysis of a previously modelled dataset.
