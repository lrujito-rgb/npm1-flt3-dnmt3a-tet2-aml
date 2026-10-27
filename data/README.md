# data/

`data/raw/` is created by `src/build_data.py` and is **not** version-controlled (see `.gitignore`).
All inputs are re-downloaded from their public sources:

| Cohort | Source | Identifier | Licence / terms |
|---|---|---|---|
| TCGA-LAML | cBioPortal REST API | `laml_tcga_pub` | TCGA open-access tier; cBioPortal terms of use |
| Beat AML | cBioPortal REST API | `aml_ohsu_2022` (+ `aml_ohsu_2018` for ITD-classifier validation only) | cBioPortal terms of use; cite Bottomly et al. 2022 |
| AMLSG | GitHub `gerstung-lab/AML-multistage/data` | commit `a3f96ba` (pinned default; override via `AMLSG_REF`) | see source repository; cite Papaemmanuil 2016 and Gerstung 2017 |
| FLT3 reference sequence | UniProt | P36888 | CC BY 4.0 |

Raw and patient-level derived data are deliberately not redistributed here; the scripts regenerate them.
Checksums of the exact files analysed are in `results/data_provenance.json`.
