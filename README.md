# DTA Meta-Analysis Pro v4.9.2

**A browser-based diagnostic test accuracy meta-analysis tool implementing bivariate GLMM and HSROC models.**

## Quick start

1. Download `dta-pro-v3.7.html`
2. Open in Chrome, Firefox, Safari, or Edge
3. Click **Datasets** > **Demo Data** to load sample data
4. Click **Run Analysis**

No installation, no server, no programming required.

## Features

| Feature | Description |
|---------|-------------|
| Bivariate GLMM | Reitsma et al. (2005) joint model for sensitivity and specificity |
| HSROC | Rutter & Gatsonis (2001) hierarchical summary ROC model |
| Auto model selection | AIC-based model averaging when neither model dominates |
| SROC curve | Summary ROC with confidence and prediction regions |
| Forest plots | Paired sensitivity/specificity forest plots with study weights |
| QUADAS-2 | Full quality assessment: 4 risk-of-bias + 3 applicability domains |
| Deeks' funnel test | Publication bias assessment (weighted regression) |
| Clinical utility | PLR, NLR, post-test probabilities, Fagan nomogram, NND |
| Bootstrap CIs | Non-parametric bootstrap confidence intervals |
| Influence diagnostics | Cook's distance, leave-one-out analysis |
| GRADE-DTA | Certainty of evidence assessment |
| PRISMA-DTA | Automated reporting checklist |
| PDF report | One-click comprehensive PDF generation |
| R code export | Generates reproducible R script using mada package |
| Dark/light theme | Full theme support with WCAG AA contrast |

## Validation

Validated against R `mada` package v0.5.12 (bivariate GLMM via `reitsma()`).

| Dataset | Studies | Max Sens diff | Max Spec diff | Max AUC diff | Result |
|---------|---------|---------------|---------------|--------------|--------|
| Dementia (MMSE) | 33 | < 0.001 | < 0.001 | < 0.001 | PASS |
| Scheidler MRI | 8 | < 0.001 | < 0.001 | < 0.001 | PASS |
| CD64 Sepsis | 10 | < 0.001 | < 0.001 | < 0.001 | PASS |
| Glas FDG-PET | 9 | benchmark | benchmark | benchmark | PASS |

**30/30 validation checks passed.** See `tests/validate_dta_pro.R` for the full validation script.

## Repository structure

```
dta-meta-analysis-pro/
+-- dta-pro-v3.7.html              # Main application
+-- README.md                       # This file
+-- LICENSE                         # MIT License
+-- tests/
|   +-- validate_dta_pro.R         # R validation script (mada package)
|   +-- r_validation_results.json  # Validation output
+-- paper/
|   +-- DTA_Pro_Manuscript_PLOS_ONE.md  # Manuscript
+-- figures/
    +-- (screenshots for manuscript)
```

## Statistical methods

The bivariate GLMM jointly models logit-transformed sensitivity and specificity:

```
[logit(sens_i)]     [mu1]     [u1i]
[logit(spec_i)]  =  [mu2]  +  [u2i]  +  epsilon_i

where (u1i, u2i) ~ N(0, Sigma),  Sigma = [[tau1^2, rho*tau1*tau2], [rho*tau1*tau2, tau2^2]]
```

Estimation uses maximum likelihood via Fisher scoring. For studies with k < 30, the Hartung-Knapp-Sidik-Jonkman correction is applied for more accurate confidence interval coverage.

The HSROC model (Rutter & Gatsonis 2001) is mathematically equivalent when no covariates are included (Harbord et al. 2007).

## Key references

- Reitsma JB et al. Bivariate analysis of sensitivity and specificity. *J Clin Epidemiol*. 2005;58:982-990.
- Rutter CM, Gatsonis CA. A hierarchical regression approach to meta-analysis. *Stat Med*. 2001;20:2865-2884.
- Harbord RM et al. A unification of models for meta-analysis of diagnostic accuracy studies. *Biostatistics*. 2007;8:239-251.
- Macaskill P et al. Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy, Chapter 10. 2010.

## Requirements

- Any modern browser (Chrome 90+, Firefox 90+, Safari 15+, Edge 90+)
- For R validation: R >= 4.0, mada >= 0.5.11, jsonlite

## License

MIT License. See [LICENSE](LICENSE).

## Citation

If you use this tool in your research, please cite:

> Ahmad M, Kumar N, Dar B, Khan L, Woo A. DTA Meta-Analysis Pro v4.9.2: a browser-based tool for diagnostic test accuracy meta-analysis. *PLOS ONE*. 2026. [DOI pending]
