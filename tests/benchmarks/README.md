# Benchmark Datasets for DTA Meta-Analysis Pro Validation

These datasets are used to validate DTA Pro's bivariate GLMM implementation
against the R mada package (v0.5.12).

## Datasets

| File | Description | Studies (k) | Source |
|------|-------------|-------------|--------|
| dementia_mmse.csv | MMSE screening for dementia | 33 | mada::Dementia / Mitchell 2009 |
| scheidler_mri.csv | MRI for lymph node staging | 8 | Scheidler et al., mada package |
| cd64_sepsis.csv | CD64 biomarker for sepsis | 10 | Published DTA reviews |
| glas_fdgpet.csv | FDG-PET for colorectal recurrence | 9 | Glas et al. 2003, mada::Glas |

## CSV Format

All files use the same format:
- `study`: Study identifier (first author + year)
- `tp`: True positives
- `fp`: False positives
- `fn`: False negatives
- `tn`: True negatives

## Validation

Run `../validate_dta_pro.R` using R (>=4.3.0) with the mada package (>=0.5.12):

    Rscript tests/validate_dta_pro.R

Results are saved to `../r_validation_results.json`.

## License

These datasets are from published studies and the mada R package (GPL).
They are redistributed here for validation purposes only.
