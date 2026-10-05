# Data availability

> Pseudonymized ≠ anonymous, especially for small cohorts.

This project analyses human milk from **5 donors** (plus one excluded donor). With a cohort this small, even coded records can be linked back to people, so the input data are handled as restricted.

## 1. Public data in this repository

| What | Where | Notes |
|---|---|---|
| Analysis code | `analysis_non_targeted.py`, `analysis_non_targeted.ipynb` | MIT licence |
| Synthetic test data generator | `tests/conftest.py` | made-up donors and intensities in the real file format |

## 2. Restricted input (not distributed)

| File (expected name) | Content | Why it is restricted |
|---|---|---|
| `input_data/nontargeted table 5 000 threshold.csv` | untargeted LC-MS feature table: 752 features × 55 injections (40 samples, 10 pooled QC, 1 blank, 4 excluded) | per-sample lipid profiles of identifiable donors |
| `input_data/metadata_M1_M12 2 05 2026.xlsx` | sample, sequence and donor metadata, including age, BMI, delivery mode, parity, diet, supplementation and health status | health data about individual donors |

- **Source:** course material for the *Lipidomics and Glycomics* project (InfoBioChem MSc, Gdańsk University of Technology, 2026), from the research group that collected the samples.
- **Access:** only through the course instructors or that research group. The data are not in any public repository.
- **Minimisation:** the analysis reads only the `Sample_metadata` and `Sequence` sheets. The demographic and health sheets (`Donor_metadata`, `Timepoint_metadata`) are deliberately not loaded.
- **Pseudonymisation in outputs:** donor codes in the source files look like initials. The script replaces them at run time with neutral aliases (D1–D6, assigned in sorted order), so the original codes never appear in figures, result files or the code itself. `input_data/` and `processed_data/` are in `.gitignore`.

## 3. Derived aggregate results (public)

| File | Level of detail |
|---|---|
| `results/results_summary.json`, `filtering_summary.csv`, `metadata_audit_summary.json` | counts and summary statistics |
| `results/feature_filter_metrics.csv`, `final_feature_metrics.csv` | per **feature** (m/z, RT, presence, QC RSD), aggregated over samples |
| `results/univariate_statistics_m1_vs_m12.csv` | per feature: mean log2 fold change, paired tests, FDR |
| `results/pca_loadings.csv` | per-feature loadings (no sample scores) |
| `figures/*.png` | group-level plots. Figures 05 (PCA) and 07 (heatmap) show donor-level points or z-scores under aliases D1–D6. They do not reveal raw intensities or any metadata beyond timepoint. |
| `poster/`, `report/` | same content, donor codes replaced by the aliases |

Not published: per-sample matrices (`final_matrix_*.csv`), donor-aggregated matrices, PCA sample scores, and the inferred column-metadata table.

## 4. What can and cannot be reproduced without the restricted data

| Step | Reproducible from this repository? |
|---|---|
| Code paths: parsing, filtering, imputation, PCA, statistics, plotting | **Yes**, on synthetic data (`pytest`) |
| Blank / QC / presence filtering counts (752 → 575 → 160 → 94) | **No**, needs the feature table |
| QC RSD distribution, PCA, volcano plot, heatmap | **No**, needs the feature table and sample metadata |
| Paired M1 vs M12 statistics and FDR | **No**, needs the feature table and sample metadata |
| Interpretation and figures as published | Yes, from `results/` and `figures/` |

With access to the restricted files, the published results are reproduced exactly (checked byte-for-byte on 2026-10-05):

```bash
python analysis_non_targeted.py --input-dir input_data --output-dir .
```
