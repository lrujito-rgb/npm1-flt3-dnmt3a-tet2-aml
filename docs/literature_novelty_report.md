# Literature verification and novelty assessment
**Question:** In adults with NPM1-mutated AML carrying FLT3-ITD or FLT3-TKD, do DNMT3A and/or TET2 co-mutations modify the association between FLT3 subtype and overall survival?
**Date of searches:** 26 Oct 2026. **Sources:** PubMed E-utilities, Europe PMC / PMC full text, Crossref (DOI resolution for all 72 references), Unpaywall (ELN 2022 PDF), and web search (3 queries). OpenAlex was not used.
Reference keys refer to `references_verified.csv`. All 72 PMIDs exist and all DOIs resolve in Crossref, with titles concordant.

## Scope and deviations (read first)
- **Full text read (18):** ELN 2022, WHO 2022 (Khoury), Papaemmanuil 2016, Gerstung 2017, HARMONY (Hernández-Sánchez 2026), Yao 2024 (BCJ), Cocciardi 2025, Ni 2025, Arora 2025, Jahn 2022, Loghavi 2014, Oñate 2022, Wu 2024, Eisfeld 2020, Cocciardi 2019, Stone 2017, Smith 2012, Daver 2019.
- **Not accessible (abstract only, or less):**
  - **Boddu 2017:** full text unavailable. The abstract does not mention DNMT3A or TET2.
  - **Voso 2020:** abstract only.
  - **ELN 2017:** full text unavailable. Table 5 was read from the PMC table page via web search.
  - **ICC (Arber 2022):** abstract only. The 10% blast rule was confirmed from the ELN 2022 text.
  - **Bezerra 2020:** no PubMed abstract and full text blocked (HTTP 403). Only the title and a publisher-page snippet were seen, so how FLT3 subtype was handled is unknown.
  - **ADMIRAL (Perl 2019):** abstract only. TKD eligibility was **not** confirmed from retrieved text.
  - **Altman & Bland 1995 and 2003:** records only; there is no PubMed abstract.
- **Supplementary materials were not accessed for any paper.** This matters most for:
  - Papaemmanuil 2016 Table S10 (full interaction model) and Fig. S15 (NPM1 subgroups);
  - the HARMONY supplement (Figs S4–S19);
  - Cocciardi 2025 Tables S7–S8;
  - Ni 2025 Fig. 2D.
  A TKD × DNMT3A/TET2 result could sit in these supplements, so the novelty verdict below is **provisional**.
- **Search coverage:**
  - **Systematic PubMed searches:** three queries returning 35, 127 and 96 records, plus a TKD/NPM1 title-abstract query (2015–2026) returning 58. All titles were screened and NPM1-relevant abstracts read.
  - **Not searched:** Embase, conference abstracts (ASH/EHA), and non-English full texts.
- **Reference count:** 72 verified references, above the 30–45 target, so the manuscript can choose. A suggested core set is marked in §(e).

---

## (a) ELN and WHO/ICC rules relevant to this study (paraphrased)
**ELN 2022 [Dohner2022, Table 6; read from full text]**
- **Favorable:** mutated NPM1 *without* FLT3-ITD.
- **Intermediate:** mutated NPM1 *with* FLT3-ITD. Wild-type NPM1 with FLT3-ITD is also intermediate if there are no adverse lesions.
- **Allelic ratio dropped:** the FLT3-ITD allelic ratio is no longer used. FLT3-ITD AML is intermediate whatever the ratio or NPM1 status. Reasons given were assay standardisation problems, the effect of midostaurin, and the growing role of MRD.
- **FLT3-TKD:** not used for risk assignment at all. ELN 2022 says only that FLT3 screening should cover both ITD and TKD. **NPM1/FLT3-TKD is therefore favorable because it is "NPM1 without FLT3-ITD", not because of TKD.**
- **Adverse myelodysplasia-related genes:** ASXL1, BCOR, EZH2, RUNX1, SF3B1, SRSF2, STAG2, U2AF1, ZRSR2, plus TP53 (VAF ≥10%). **DNMT3A and TET2 are not on the list and do not change ELN risk.**
- A footnote says MR-gene mutations should not, for now, count as adverse when they co-occur with favorable-risk subtypes.
- NPM1-mutated AML with adverse-risk cytogenetics is adverse.
- Risk is mainly based on intensively treated patients and may change with MRD.

**ELN 2017 [Dohner2017, Table 5 via PMC table page]**
- **Favorable:** mutated NPM1 without FLT3-ITD, or with FLT3-ITD^low (allelic ratio <0.5).
- **Intermediate:** NPM1-mut with FLT3-ITD^high (≥0.5).
- **Adverse:** NPM1-wt with FLT3-ITD^high.
- FLT3-TKD, DNMT3A and TET2 are not risk-defining. The adverse gene mutations are RUNX1, ASXL1 and TP53.
- **Implication:** in TCGA and Beat AML, "ELN 2017" versus "ELN 2022" assignment differs for NPM1/ITD-low patients. State which version you use.

**WHO 2022 [Khoury2022, full text]:** AML with NPM1 mutation is a genetically defined entity and can be diagnosed at any blast count.

**ICC 2022 [Arber2022; rule confirmed via Dohner2022 text]:** AML-defining recurrent genetic abnormalities, including NPM1, need ≥10% blasts in marrow or blood. NPM1 cases with 10–19% blasts are classed as AML, not MDS/AML.

---

## (b) Prior studies relevant to the novelty question
| Study (key) | Cohort / n | What was analysed | Key result as reported | TKD × DNMT3A/TET2 tested? |
|---|---|---|---|---|
| Papaemmanuil 2016 (AMLSG; **same data as your AMLSG set**) | 1540 (NPM1-mut 436) | Multivariable Cox with gene–gene interaction terms | NPM1–FLT3-ITD–DNMT3A triple genotype in 93/1540 (6%); interaction HR 1.5 (1.2–1.9), P=0.0002, q=0.004. FLT3-ITD harm was "most clinically relevant" with both NPM1 and DNMT3A, and held regardless of ITD ratio. NPM1–FLT3-TKD interaction in 53 (3%), HR 0.7 (0.6–0.9), q=0.06. MLL-PTD–FLT3-TKD interaction poorer than expected (q=0.008). | **No** TKD×DNMT3A or TKD×TET2 term among the reported (q<0.1) effects. Full model (Table S10) not accessed. |
| Gerstung 2017 (same 1540) | 1540 | Knowledge-bank multistage models including interaction terms | No specific NPM1/FLT3/DNMT3A effect reported in main text. Notes most trials are underpowered for interactions involving genes mutated in <20% of patients. | No (main text) |
| Loghavi 2014 (MDACC) | 178 with DNMT3A/FLT3/NPM1 mutations; triple group n=35 (24 ITD, 9 TKD, 2 both) | Clinical features, OS/EFS by mutation group | DNMT3A independently associated with OS (p=0.026). **Within triple mutants, OS did not differ between FLT3-ITD and FLT3-TKD (p=0.822).** | **Closest direct precedent**: ITD vs TKD compared within DNMT3A-mut/NPM1-mut, but very small, descriptive, and no formal interaction test or TET2 |
| Boddu 2017 (MDACC) | TKD+NPM1 21; TKD-only 18; NPM1-only 117; ITD+NPM1 107 | NPM1/TKD vs NPM1-only | TKD/NPM1 had better RFS (median NR vs 18.3 mo, P=.03), trend for OS, and was independently favorable. | Not in abstract; full text not accessed |
| Perry 2018 | 126 (NPM1/TKD 12) | TKD with or without NPM1 | NPM1/TKD favorable | No (abstract) |
| Bacher 2008 (MLL) | 3082 (TKD 147) | TKD prognosis by context | No overall TKD effect; additional favorable EFS effect in NPM1- or CEBPA-mutated AML | No (pre-NGS era) |
| Voso 2020 (RATIFY TKD) | 163 TKD | Midostaurin effect; co-variables | EFS better with midostaurin (45.2% vs 30.1%, P=.044), OS similar. NPM1 independently favorable within TKD. | Not in abstract; full text not accessed |
| Jahn 2022 (RATIFY genomics) | 475 FLT3-mut (TKD 116) | Random survival forests; pairwise interactions | FLT3 type (ITD-high vs TKD) prognostic. **No relevant interaction for FLT3:DNMT3A:NPM1** with FLT3 pooled. | Not by FLT3 subtype |
| Arora 2025 (MDACC) | 124 TKD (ITD excluded) + 126 NPM1mut/TKDwt | NPM1 within TKD AML | NPM1 associated with better 3-yr OS (intensive therapy: 74% vs 40%, p=.03) | **No**: DNMT3A/TET2 given only as baseline frequencies |
| Cocciardi 2025 (AMLSG 09-09) | 568 NPM1-mut (TKD 25.9%; ITD 17.3%) | Multivariable EFS/OS in ELN-2022 favorable NPM1 (n=470, which includes NPM1/TKD) | DNMT3A R882 adverse for EFS and OS. FLT3-TKD non-significant covariate (OS HR 1.13, 0.77–1.64). | **No** (TKD as covariate only). Validation used the 435 Papaemmanuil NPM1 patients, **overlapping your AMLSG data** |
| HARMONY (Hernández-Sánchez 2026) | Training 1001; validation 762 + 585 | Hierarchical classification by FLT3-ITD, DNMT3A, IDH, TET2 | ITD+DNMT3A: 2-yr OS 29.1%. ITD-absent DNMT3A-mut IDH-wt: 77.2%. ITD-absent DNMT3A-wt TET2-mut: 56.9% vs TET2-wt 73.5%. | **No**: TKD (17%) not used. The "ITD-absent" branch mixes TKD and FLT3-wt. Supplement not accessed |
| Ni 2025 (BCJ) | 884 AML; DNMT3A-mut ≈143 | Co-mutations within DNMT3A-mut subgroups | In DNMT3Amut/NPM1mut/FLT3-ITDwt, TET2-mut had the worst outcome. FLT3-TKD was one stratifier among signalling genes. | Partly (TKD as a separate stratum); no interaction test; subgroup ≈47 (inferred; provisional) |
| Yao 2024 (BCJ) | 192 NPM1-mut (178 followed) | Co-mutation patterns | ITD and DNMT3A adverse; triple worst. TET2 and FLT3-TKD showed non-significant trends. | No |
| Wu 2024 | 91 NPM1mut/ITDwt (TKD 27) | Co-mutations in ITD-wt NPM1 | TET1/2 worse RFS (p=0.038), not OS. TKD and DNMT3A no significant effect. | No |
| Tian 2014 | 373 CN-AML | TET2 in NPM1 subgroups | TET2 shorter OS in NPM1+/FLT3-ITD− (9.5 vs 32.2 mo, p=0.013) | No (ITD-negative pools TKD and FLT3-wt) |
| Oñate 2022 (CETLAM) | 164 NPM1-mut | DNMT3A by ITD absent/low/high | DNMT3A did not change OS or ITD-ratio stratification; delayed MRD clearance | No (ITD only) |
| Duan 2024 | 234 NPM1-mut | Co-mutations and MRD | DNMT3A independent adverse (OS HR 2.2). ITD+DNMT3A poorest. | No (ITD only) |
| Gil 2026 (PETHEMA) | 1360 (688 intensive) | Co-mutation signatures | DNMT3A inferior RFS/EFS, not OS. "FLT3-OTHER" associated with improved OS. | No (abstract; FLT3-OTHER definition not checked) |
| Bezerra 2020 | 507 de novo AML | DNMT3A/NPM1/FLT3 co-occurrence | Triple co-occurrence adverse (title) | Unknown: content not accessible |

**Q3: Papaemmanuil/Gerstung and the triple interaction (confirmed from full text).**
- Papaemmanuil 2016 reported the three-way NPM1–DNMT3A–FLT3-ITD interaction in the same AMLSG cohort you plan to use.
- It was the most frequent three-gene co-occurrence (93/1540).
- Table 2 gives HR 1.5 (1.2–1.9), P=0.0002, q=0.004, with interaction components of 1.1–1.2.
- Their wording, in short: the deleterious effect of FLT3-ITD was "most clinically relevant in patients with concomitant NPM1 and DNMT3A mutations".
- Gerstung 2017 re-used the same 1540 patients. Its main text does not restate this interaction.
- **Consequence:** any DNMT3A × FLT3-ITD signal in the AMLSG arm is **not an independent replication** of Papaemmanuil. Report AMLSG as overlapping or derivation data, or show cohort-specific estimates.

---

## (c) Novelty verdict and suggested knowledge-gap wording
**Verdict (provisional, pending supplements):** None of the main texts or abstracts I read formally tests whether DNMT3A and/or TET2 **modify the OS difference between FLT3-ITD and FLT3-TKD** within NPM1-mutated AML. None reports DNMT3A/TET2 effects specifically within **NPM1/FLT3-TKD** patients as a primary analysis either.

The question is only **partly novel**, for four reasons:
1. **Loghavi 2014** already compared ITD and TKD within DNMT3A/FLT3/NPM1 triple mutants and found no OS difference (n=35, only 9 TKD-only).
2. In many datasets, the NPM1/ITD-absent subgroups analysed for TET2 (HARMONY, Ni 2025, Wu 2024, Tian 2014) or DNMT3A (Cocciardi 2025) **include TKD patients, but pooled with FLT3-wt**.
3. The AMLSG data have already been analysed for gene–gene interactions (Papaemmanuil 2016), and Cocciardi 2025 reused its 435 NPM1 patients.
4. With roughly 66 NPM1/TKD-only patients with OS across cohorts (per your memory note; not re-verified here), the analysis can estimate effects but cannot confirm or exclude moderate interactions [Altman1995; Gerstung2017].

**Recommended wording:**
> *What is known:* In NPM1-mutated AML, FLT3-ITD is associated with inferior outcome. This association is strongest when DNMT3A is co-mutated, as shown by a three-way NPM1–DNMT3A–FLT3-ITD interaction in AMLSG patients and by several independent cohorts [Papaemmanuil2016; Loghavi2014; Yao2024; Duan2024; HernandezSanchez2026], although not uniformly [Onate2022; Jahn2022]. FLT3-TKD is not part of ELN risk assignment [Dohner2017; Dohner2022]; observational data suggest NPM1/FLT3-TKD is favorable [Bacher2008; Boddu2017; Perry2018; Voso2020; Arora2025], with conflicting reports [Whitman2008; Li2023]. DNMT3A (particularly R882) and TET2 have been associated with poorer outcome in NPM1-mutated, FLT3-ITD-negative AML in some cohorts [Cocciardi2025; HernandezSanchez2026; Ni2025; Metzeler2011; Tian2014] but not others [Gaidzik2012; Gaidzik2013; Wu2024 (OS)].
> *What is not known:* Whether DNMT3A or TET2 co-mutation modifies the prognostic contrast between FLT3-ITD and FLT3-TKD in NPM1-mutated AML has not, to our knowledge, been formally tested with an interaction model. Prior analyses either pooled FLT3-TKD with FLT3-wild-type patients, modelled FLT3-TKD only as a covariate, or compared FLT3 subtypes within very small triple-mutant groups.

Avoid "first study to show" phrasing. Prefer "to our knowledge, not formally tested" and present the work as hypothesis-generating. Before submission, check Papaemmanuil Table S10/Fig S15, the HARMONY supplement, Cocciardi 2025 Tables S7–S8, and the Boddu 2017 full text for TKD-stratified co-mutation analyses.

**How the new analysis differs:**
1. It tests the FLT3-subtype × DNMT3A/TET2 interaction explicitly, ideally with stratum-specific and joint effects [Knol2012; Altman2003].
2. It keeps NPM1/FLT3-TKD separate from NPM1/FLT3-wt.
3. It uses public multi-cohort data (TCGA, Beat AML, AMLSG) with cohort stratification.
4. It reports OS specifically.

Limitations to state:
- sparse cells (Firth correction [Heinze2001]; EPV [Peduzzi1995; Vittinghoff2007]);
- heterogeneous treatment, mostly pre-midostaurin [Stone2017];
- AMLSG overlap with Papaemmanuil, Gerstung and Cocciardi 2025;
- no MRD data [Ivey2016];
- no inference of clonal order from VAF [Miles2020].

---

## (d) Draft-manuscript claims checked against the literature
| Draft claim | Verdict | Evidence / suggested correction |
|---|---|---|
| "NPM1/FLT3-TKD is favourable because TKD is favourable per ELN" | **Contradicted** | ELN 2017 and 2022 do not use FLT3-TKD for risk. NPM1/TKD is favorable only as "NPM1 without FLT3-ITD" [Dohner2017; Dohner2022]. TKD's favorable association in NPM1 comes from observational studies [Bacher2008; Boddu2017; Perry2018; Voso2020; Arora2025]. The meta-analysis shows no overall TKD effect on OS (HR 0.98) [Li2023], and one study found worse DFS [Whitman2008]. |
| "DNMT3A/TET2 place patients in the adverse ELN category" | **Contradicted** | Neither gene is in the ELN 2022 MR-gene list or the ELN 2017 adverse list [Dohner2022; Dohner2017]. Proposals to upgrade DNMT3A risk come from individual studies (e.g. intermediate-risk DNMT3A in Eisfeld2020; Ni2025; HARMONY), not guidelines. |
| "FLT3-TKD signals with minimal STAT5" | **Context-dependent; needs qualifying** | True in models without NPM1c [Choudhary2005; Grundler2005]. But with NPM1c, FLT3-TKD activates STAT5, mislocalises to the ER, and aberrant STAT5 activation appears in patient cells with both mutations [Rudorf2019]. In an NPM1-mutant manuscript, do not state "minimal STAT5" without this qualification. |
| "Type II FLT3 inhibitors are ineffective against TKD" | **Broadly supported; needs nuance** | Type II inhibitors (quizartinib, sorafenib) bind the inactive conformation and generally do not target TKD, whereas type I inhibitors (gilteritinib, crenolanib) act on both ITD and TKD [Daver2019]; acquired D835 mutations confer quizartinib resistance [Smith2012]. However, sensitivity varies by substitution (D835E/N relatively sensitive; D835Y/V/I/F resistant) [Smith2015]. QuANTUM-First enrolled FLT3-ITD only [Erba2023]. |
| "Midostaurin benefit applies only to FLT3-ITD" (if stated) | **Contradicted** | RATIFY stratified by TKD and observed longer OS with midostaurin in the TKD subgroup (not powered) [Stone2017]. TKD post hoc analysis: EFS improved, OS similar [Voso2020]. |
| "ELN 2022 uses FLT3-ITD allelic ratio" (if stated) | **Contradicted** | Allelic ratio was removed; NPM1/ITD is intermediate at any ratio [Dohner2022]. |
| "DNMT3A is consistently adverse in NPM1-mut AML" (if stated) | **Overstated** | Conflicting evidence: adverse [Loghavi2014; Duan2024; Cocciardi2025 (R882); Herold2020]; no overall effect [Gaidzik2013; Onate2022]; RFS/EFS but not OS [Gil2026]. Effects differ by R882 vs non-R882 and age [Gaidzik2013; Marcucci2012]. |
| "TET2 is adverse in NPM1-mut AML" (if stated) | **Mixed evidence** | Adverse: Metzeler2011 (ELN favorable CN-AML), Tian2014, HARMONY (intermediate step), Ni2025, Poonsombudlert2025 (EFS). No effect: Gaidzik2012. RFS only, not OS: Wu2024. |
| "DNMT3A/TET2 VAF shows they are founder events" (if stated) | **Methodologically weak** | Bulk VAF cannot definitively determine clonal co-occurrence or order [Miles2020; Morita2020]. Founder status is supported biologically by pre-leukaemic HSC and CH data [Shlush2014; CorcesZimmerman2014; Jaiswal2014; Genovese2014], not by your VAFs. |
| "AMLSG data independently validate the triple-mutant effect" (if stated) | **Incorrect framing** | Same patients as Papaemmanuil2016, Gerstung2017, and the Cocciardi2025 validation cohort. |

---

## (e) Suggested core reference set (about 45)
- **Guidelines and classification:** Dohner2017, Dohner2022, Khoury2022, Arber2022.
- **Data and tools:** TCGA2013, Tyner2018, Bottomly2022, Papaemmanuil2016, Gerstung2017, Cerami2012, Gao2013, deBruijn2023.
- **NPM1/FLT3:** Falini2005, Thiede2006, Schlenk2008, Gale2008.
- **TKD:** Mead2007, Bacher2008, Whitman2008, Boddu2017, Voso2020, Arora2025, Li2023, Jahn2022.
- **Co-mutation:** Loghavi2014, Bezerra2020 (cite with caution; content not verified), Onate2022, Yao2024, HernandezSanchez2026, Cocciardi2025, Ni2025, Wu2024.
- **DNMT3A/TET2:** Ley2010, Metzeler2011, Gaidzik2012, Gaidzik2013.
- **Clonal hematopoiesis and clonality:** Jaiswal2014, Shlush2014, Miles2020.
- **Biology and therapy:** Rudorf2019, Daver2019, Smith2015, Stone2017, Erba2023, Perl2019.
- **Statistics and reporting:** Peduzzi1995, Vittinghoff2007, Heinze2001, Knol2012, Altman2003, Altman1995, vonElm2007.
- **Not core:** the remaining entries in the CSV are verified alternatives.
