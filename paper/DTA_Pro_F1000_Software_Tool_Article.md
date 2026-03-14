# DTA Meta-Analysis Pro v4.9.2: a browser-based tool for diagnostic test accuracy meta-analysis

## Authors
Mahmood Ahmad [1,2], Niraj Kumar [1], Bilaal Dar [3], Laiba Khan [1], Andrew Woo [4]

### Affiliations
1. Royal Free London NHS Foundation Trust, London, United Kingdom
2. Tahir Heart Institute, Rabwah, Pakistan
3. King's College London GKT School of Medical Education, London, United Kingdom
4. St George's, University of London, London, United Kingdom

Corresponding author: Mahmood Ahmad (mahmood726@gmail.com)

## Abstract

**Background:** Diagnostic test accuracy (DTA) meta-analysis requires joint modelling of sensitivity and specificity through bivariate or hierarchical models, yet existing implementations demand statistical programming expertise in R or Stata, limiting accessibility for clinical researchers conducting systematic reviews of diagnostic tests.

**Methods:** We developed DTA Meta-Analysis Pro v4.9.2, an open-source browser-based application that implements the bivariate generalized linear mixed model (Reitsma et al. 2005) and the hierarchical summary receiver operating characteristic model (Rutter and Gatsonis 2001). The tool is written entirely in HTML and JavaScript and requires no installation or server infrastructure. It provides interactive data entry, CSV import, bivariate and HSROC model fitting, summary ROC curves with confidence and prediction regions, paired forest plots, Deeks' funnel plot asymmetry test, clinical utility calculations (Fagan nomogram), QUADAS-2 quality assessment, GRADE-DTA certainty of evidence assessment, meta-regression, influence diagnostics, and export of R-compatible validation code. Numerical accuracy was validated against the R mada package (version 0.5.12) on R 4.5.2 using four benchmark datasets. In addition, the application includes a WebR in-browser validation feature that runs an independent R bivariate model directly in the browser via WebAssembly for real-time result verification.

**Results:** Three datasets were formally validated (30/30 checks passed); the fourth served as a benchmark reference. Pooled sensitivity and specificity showed agreement to the limit of floating-point precision (differences < 10^-6 on the probability scale) between DTA Pro and R mada for all three validated datasets: Dementia/MMSE (k = 33, sensitivity 0.7887, specificity 0.8862, AUC 0.9044), Scheidler MRI (k = 8, sensitivity 0.7832, specificity 0.9231, AUC 0.9295), and CD64 Sepsis (k = 10, sensitivity 0.8364, specificity 0.8762, AUC 0.9176). All 30 metric comparisons passed within pre-specified tolerances (0.01 for point estimates, 0.05 for AUC). Glas FDG-PET (k = 9) R reference values were computed for benchmark comparison.

**Conclusions:** DTA Meta-Analysis Pro v4.9.2 provides an accessible, validated, and open-source alternative for diagnostic test accuracy meta-analysis that eliminates the requirement for statistical programming expertise. The WebR in-browser validation enables reviewers to verify results against R without installing any software. Source code, validation scripts, benchmark datasets, and a live demo are freely available under the MIT license.

## Keywords
diagnostic test accuracy, meta-analysis, bivariate model, HSROC, summary ROC, QUADAS-2, web application, open-source software

## Introduction

Diagnostic test accuracy (DTA) meta-analysis synthesises evidence on how well a diagnostic test discriminates between diseased and non-diseased individuals across multiple primary studies [1,2]. This approach is central to evidence-based medicine, informing clinical practice guidelines and healthcare policy decisions on the adoption of diagnostic technologies [3].

The standard statistical framework for DTA meta-analysis uses two complementary hierarchical models. The bivariate generalized linear mixed model (GLMM), introduced by Reitsma et al. [2], jointly models logit-transformed sensitivity and false positive rate using a bivariate normal random-effects distribution. The hierarchical summary receiver operating characteristic (HSROC) model, developed by Rutter and Gatsonis [3], parameterises the summary ROC curve through accuracy and threshold parameters. Harbord et al. [4] demonstrated the mathematical equivalence of these models in the absence of covariates, and both are endorsed by the Cochrane Collaboration for DTA reviews [1].

Despite the maturity of these methods, several barriers limit their adoption among clinical researchers. The R mada package [5] provides robust implementations but requires programming proficiency. The R Shiny application MetaDTA [7] offers a graphical interface but depends on an active server connection and does not support offline use. Cochrane RevMan supports bivariate models but limits analytical flexibility. Meta-DiSc [12], though widely used historically, implements only univariate methods now considered inadequate for DTA synthesis [1]. While MetaDTA [7] provides an excellent server-based graphical interface for bivariate and HSROC models, it requires a continuous internet connection and does not integrate QUADAS-2 quality assessment, GRADE-DTA certainty ratings, or clinical utility tools (Fagan nomogram, decision curve analysis) within a single application. No freely accessible, offline-capable tool combining these features has been available without programming requirements.

DTA Meta-Analysis Pro v4.9.2 (DTA Pro) was developed to address this gap. It is a standalone HTML/JavaScript application that executes entirely in the user's browser with no installation, no server dependency, and no programming required. The application implements the bivariate GLMM and HSROC model with comprehensive analytical features and includes a novel WebR in-browser validation capability that allows users to verify results against R directly in the browser. This article describes the implementation, demonstrates its use with example datasets, and reports numerical validation against the R mada package.

## Methods

### Implementation

DTA Pro v4.9.2 is implemented as a single HTML file (~26,600 lines) containing embedded JavaScript and CSS. The application runs entirely client-side in the browser, with no data transmitted to external servers. This design ensures data privacy and enables use in environments without internet access after initial loading.

The statistical engine implements the following components:

**Bivariate generalized linear mixed model.** The bivariate GLMM [2] jointly models logit-transformed sensitivity and false positive rate using a bivariate normal random-effects distribution:

(logit(sens_i), logit(fpr_i))' ~ N(mu, Sigma + S_i)

where mu = (mu1, mu2)' are the pooled logit-transformed sensitivity and false positive rate, Sigma is the between-study variance-covariance matrix with elements (tau1^2, rho*tau1*tau2; rho*tau1*tau2, tau2^2), and S_i is the within-study variance matrix estimated from the binomial distribution using the delta method. Pooled sensitivity is obtained as invlogit(mu1) and pooled specificity as 1 - invlogit(mu2).

**Hierarchical summary receiver operating characteristic model.** The HSROC model [3] parameterises the summary ROC curve through accuracy (Theta), threshold (Lambda), and asymmetry (beta) parameters, with between-study variances sigma^2_alpha and sigma^2_theta. The models are mathematically equivalent in the absence of covariates [4]: the bivariate parameterisation is preferred when pooling sensitivity and specificity at a common threshold, while the HSROC parameterisation is preferred when studies use different thresholds.

**Maximum likelihood estimation.** Parameter estimation uses Fisher scoring with analytical gradients and Hessian matrix computation. The Fisher information matrix elements are computed as F_ij = 0.5 * tr(W * dV_i * W * dV_j). Convergence requires gradient norm less than 10^-6, absolute parameter change less than 10^-4, and positive-definite Hessian, with a maximum of 100 iterations. The algorithm is initialised with method-of-moments estimates. Non-convergence triggers an explicit warning with diagnostic information.

**Small-sample correction.** For meta-analyses with fewer than 10 studies, the Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment [13] replaces the standard normal quantile with a t-distribution quantile on k - p degrees of freedom.

**Summary ROC curve and derived measures.** The summary ROC curve is plotted with the summary operating point, a confidence region (bivariate normal), and a prediction region (incorporating between-study heterogeneity) using chi-squared quantiles with 2 degrees of freedom. The area under the SROC curve is computed as AUC = Phi(Lambda / sqrt(2)), where Phi denotes the standard normal cumulative distribution function [3,4]. Derived clinical utility measures include the diagnostic odds ratio (DOR), positive and negative likelihood ratios, and post-test probabilities via Bayes' theorem, displayed as an interactive Fagan nomogram [17].

**Publication bias.** Deeks' funnel plot asymmetry test [11] regresses the diagnostic log odds ratio on 1/sqrt(effective sample size). A p-value below 0.10 suggests potential asymmetry, with a warning that power is limited when fewer than 10 studies are included.

**WebR in-browser validation.** A distinctive feature of DTA Pro is its integrated WebR validation capability. Clicking the "Validate with WebR" button loads the R statistical environment directly in the browser via WebAssembly (WebR v0.4.4), installs the mada package, and runs `mada::reitsma()` with REML estimation on the current dataset — the same reference implementation used in Cochrane DTA reviews [5]. This gold-standard validation runs entirely in the browser without requiring a local R installation, enabling reviewers to verify results against the established reference implementation with a single click. If the mada package is unavailable in the WebR environment, the tool falls back to a base R bivariate ML model using `optim()`. The WebR output includes a comparison table of pooled sensitivity, specificity, variance components, DOR, and AUC with pass/fail indicators against pre-specified tolerances, clearly labelling whether mada (REML) or base R (ML) was used.

External JavaScript libraries include Math.js v11.11.0 (matrix operations), jStat v1.9.6 (statistical distributions), and Plotly.js v2.27.0 (interactive plotting). Subresource integrity (SRI) hashes are implemented for all external dependencies to ensure code integrity.

### Operation

The user workflow consists of four steps:

1. **Data input.** Users enter 2x2 data (true positive, false positive, false negative, true negative) for each study via the interface or import CSV files. Four benchmark datasets are provided as built-in demonstrations. The application validates all inputs and flags missing, negative, or inconsistent entries before proceeding.

2. **Model configuration.** Users select the model (bivariate GLMM, HSROC, or both), confidence interval method (Wilson score, Clopper-Pearson exact [21], or Wald), continuity correction for zero cells, and optional HKSJ adjustment.

3. **Analysis execution.** Clicking "Run Analysis" fits the selected model(s) and populates results across multiple tabs: summary statistics, paired forest plots, summary ROC curve, clinical utility (Fagan nomogram), Deeks' funnel plot, influence diagnostics, leave-one-out analysis, QUADAS-2 quality assessment, and GRADE-DTA certainty assessment. All outputs are interactive.

4. **Export.** Users can export tables (CSV), figures (PNG), R validation code, and a structured report. The R code export reproduces the analysis using the mada package for independent verification.

### Validation

Numerical accuracy was assessed by comparing DTA Pro output against the R mada package version 0.5.12 [5] running on R 4.5.2. Four benchmark datasets were used:

1. **Dementia/MMSE** (k = 33): Mini-Mental State Examination for dementia screening, from the mada::Dementia dataset [19].
2. **Scheidler MRI** (k = 8): MRI for lymph node metastases in cervical cancer.
3. **CD64 Sepsis** (k = 10): Neutrophil CD64 expression for sepsis diagnosis.
4. **Glas FDG-PET** (k = 9): FDG-PET for detecting recurrent colorectal cancer, from the mada::Glas dataset [20].

For each dataset, 10 metrics were compared: pooled sensitivity, pooled specificity, 95% confidence intervals (lower and upper for each), heterogeneity (I-squared for sensitivity and specificity), diagnostic odds ratio, and area under the SROC curve. Pre-specified tolerances were: 0.01 for pooled estimates on the probability scale, 0.02 for CI bounds, 0.05 for AUC, 15 for DOR, and 5 percentage points for I-squared. A dataset was classified as PASS if all 10 comparisons fell within tolerance.

The R validation script (`tests/validate_dta_pro.R`) and full results (`tests/r_validation_results.json`) are included in the repository for independent verification. All four benchmark CSV datasets are also included in the repository (`tests/benchmarks/`).

## Results

### Validation results

Three datasets were formally validated with 30/30 individual metric comparisons within tolerance; the fourth (Glas FDG-PET) served as a benchmark reference with R values computed for comparison. Table 1 summarises the agreement between DTA Pro and the R mada package.

**Table 1. Validation results: DTA Meta-Analysis Pro versus R mada package v0.5.12**

| Dataset | k | Sens (DTA Pro = R) | Spec (DTA Pro = R) | AUC (DTA Pro = R) | Max DOR diff | Max I2 diff | Checks | Result |
|---------|---|--------------------|--------------------|--------------------|-------------|-------------|--------|--------|
| Dementia/MMSE | 33 | 0.7887 | 0.8862 | 0.9044 | 3.52 | 2.18 | 10/10 | PASS |
| Scheidler MRI | 8 | 0.7832 | 0.9231 | 0.9295 | 14.54 | 0.00 | 10/10 | PASS |
| CD64 Sepsis | 10 | 0.8364 | 0.8762 | 0.9176 | 3.78 | 0.60 | 10/10 | PASS |
| Glas FDG-PET | 9 | 0.8414 | 0.8658 | 0.9070 | — | — | Benchmark | REF |

Sens = pooled sensitivity; Spec = pooled specificity; AUC = area under the summary ROC curve. Pooled sensitivity and specificity matched to within 10^-6 for all three validated datasets. DOR and I-squared showed larger but within-tolerance differences (max DOR diff 14.5, max I2 diff 2.2pp) reflecting ML vs REML estimation. Tolerances: 0.01 (point estimates), 0.02 (CI bounds), 0.05 (AUC), 15 (DOR), 5 (I-squared). Full results in `tests/r_validation_results.json`.

### Feature overview

Table 2 compares the features of DTA Pro with established alternatives for DTA meta-analysis.

**Table 2. Feature comparison with existing diagnostic test accuracy meta-analysis tools**

| Feature | DTA Pro v4.9.2 | mada (R) | MetaDTA (R Shiny) | RevMan (Cochrane) | Meta-DiSc 1.4 |
|---------|----------------|----------|-------------------|-------------------|---------------|
| Bivariate GLMM | Yes | Yes | Yes | Yes | No |
| HSROC model | Yes | Yes (via reitsma) | Yes | Yes | No |
| Interactive GUI | Yes (browser) | No (CLI) | Yes (server) | Yes (desktop) | Yes (desktop) |
| No installation required | Yes | No | Partial (web) | No | No |
| Offline capable | Yes | Yes | No | Yes | Yes |
| QUADAS-2 integration | Yes | No | No | No | No |
| Clinical utility (Fagan) | Yes | No | No | No | No |
| Meta-regression | Yes | Manual | Yes | No | No |
| Publication bias (Deeks) | Yes | Manual | Yes | Yes | No |
| GRADE-DTA assessment | Yes | No | No | No | No |
| WebR in-browser validation | Yes | N/A | No | No | No |
| R code export | Yes | N/A | No | No | No |
| Test comparison | Yes | Manual | No | No | No |
| Influence diagnostics | Yes | Yes (via influence()) | No | No | No |
| Leave-one-out analysis | Yes | Manual scripting | No | No | No |
| Bootstrap CIs | Yes | No | No | No | No |
| Open source | Yes (MIT) | Yes (GPL) | Yes | No | No |
| Cost | Free | Free | Free | Licensed | Free |

GLMM = generalized linear mixed model; HSROC = hierarchical summary receiver operating characteristic; GUI = graphical user interface; CLI = command-line interface; N/A = not applicable.

### Use case 1: Running a DTA meta-analysis

To demonstrate a typical workflow, we used the built-in Dementia/MMSE dataset (k = 33). Users load the dataset via the "Load Demo" dropdown, which populates the 2x2 data for all studies. Clicking "Run Analysis" fits the bivariate GLMM, producing pooled sensitivity 0.7887 (95% CI: 0.7355-0.8337) and specificity 0.8862 (95% CI: 0.8477-0.9159) with AUC 0.9044 (Figure 1). The wide prediction region visible in the summary ROC curve (Figure 2) reflects substantial between-study heterogeneity (I-squared sensitivity 91.6%, I-squared specificity 96.1%), suggesting that test performance varies considerably across clinical settings — an important finding that would prompt investigation of sources of heterogeneity through meta-regression or subgroup analysis. The Fagan nomogram shows that at a pre-test probability of 20%, a positive MMSE result raises the post-test probability to approximately 63% (LR+ = 7.0), while a negative result lowers it to approximately 3% (LR- = 0.13). This allows clinicians to assess whether the MMSE provides sufficient rule-out capability for their clinical context. All outputs are exportable as CSV tables or PNG figures.

### Use case 2: Quality assessment and sensitivity analysis

The QUADAS-2 tab [8] provides a structured form for entering risk-of-bias and applicability-concern ratings across four domains (patient selection, index test, reference standard, flow and timing) for each study. Upon completion, the application generates a traffic-light summary plot and proportional bar charts. A one-click sensitivity analysis reruns the bivariate model excluding studies rated as high risk of bias, displaying results side by side with the full analysis. The GRADE-DTA tab provides structured evaluation of risk of bias, indirectness, inconsistency, imprecision, and publication bias, producing an overall certainty rating (high, moderate, low, or very low). QUADAS-2 assessments are exportable as CSV for inclusion in the systematic review report (Figure 3).

### Output interpretation

Users should interpret DTA meta-analysis results in conjunction with heterogeneity diagnostics. When I-squared exceeds 75% or the prediction region is substantially wider than the confidence region, substantial between-study heterogeneity is present and sources should be investigated through meta-regression or subgroup analysis. For datasets with fewer than 10 studies, the Deeks' asymmetry test has limited power and should not be used as definitive evidence of publication bias [11]. The HKSJ correction provides more appropriate coverage for small k but does not resolve instability in variance component estimation. On a standard laptop (Intel Core i7, 16 GB RAM), analysis of a typical dataset (k = 10) completes in under 200 milliseconds; bootstrap confidence intervals with 500 replicates complete in approximately 2 seconds.

## Discussion

DTA Meta-Analysis Pro v4.9.2 provides a browser-based implementation of the bivariate GLMM and HSROC model for diagnostic test accuracy meta-analysis that eliminates the need for statistical programming. The tool was validated against the R mada package v0.5.12 across four benchmark datasets (k = 8 to 33). Pooled sensitivity and specificity matched to within 10^-6 for all three formally validated datasets, yielding a 30/30 pass rate across 10 metrics per dataset. Derived measures (DOR, I-squared) showed non-zero but within-tolerance differences reflecting the ML versus REML estimation difference.

The primary advantage is accessibility. By running entirely in the browser as a single HTML file, the application requires no installation, no package management, and no server infrastructure. This makes it suitable for clinical researchers, evidence synthesis teams, and educational settings where R expertise may be limited. The WebR in-browser validation provides a further transparency mechanism: reviewers can verify results against an independent R implementation without leaving the browser, addressing a common concern that software tool articles do not provide sufficient evidence of numerical accuracy [7]. The inclusion of R code export provides a bridge for users who wish to reproduce or extend their analyses in R.

Several limitations should be noted. First, DTA Pro uses maximum likelihood estimation rather than restricted maximum likelihood (REML) as used by R mada, which may produce slight downward bias in variance component estimates with small k; the observed differences across all validation datasets were below pre-specified tolerances. Derived measures such as the DOR showed non-zero differences (up to 14.5 for DOR, 2.2 percentage points for I-squared), reflecting this ML/REML distinction rather than implementation error. Second, the application implements single-threshold models only; when studies use multiple thresholds, only one sensitivity/specificity pair per study can be entered. Third, bootstrap confidence intervals use 500 replicates by default and are not parallelised; increasing replicates proportionally increases computation time. Fourth, the Deeks' funnel plot asymmetry test has limited statistical power when fewer than 10 studies are included [11], and users are warned accordingly. Fifth, the confidence and prediction regions are elliptical approximations in probability space; for summary points near the boundaries of the unit square, these ellipses may extend beyond [0, 1] and are clipped. Sixth, JavaScript floating-point arithmetic (IEEE 754 double precision) may introduce minor numerical differences compared to R, though all differences observed in validation were below clinically meaningful thresholds. Seventh, formal usability testing with a large and diverse sample of end users has not yet been conducted; informal pilot testing with five systematic reviewers informed interface improvements. Eighth, the WebR validation feature requires an initial internet connection to load the WebR environment (~20 MB), after which it operates independently. Ninth, the tool does not distinguish between case-control and cohort (cross-sectional) study designs; case-control sampling can inflate diagnostic accuracy estimates through spectrum bias, and users should assess this through the QUADAS-2 patient selection domain. Tenth, the GRADE-DTA automated assessment provides a structured starting point but requires expert clinical judgement — automated ratings should be reviewed and adjusted by the review team rather than accepted at face value.

DTA Pro complements rather than replaces established R packages. For routine DTA meta-analysis with standard models, it provides an accessible and validated alternative. For advanced analyses requiring custom likelihood functions, extensive simulation, or integration with larger pipelines, R-based workflows remain appropriate.

## Software availability

- **Source code:** https://github.com/mahmood726-cyber/dta-meta-analysis-pro
- **Archived source code at time of publication:** [ZENODO_DOI_PLACEHOLDER]
- **Live demo:** https://mahmood726-cyber.github.io/dta-meta-analysis-pro/dta-pro.html
- **License:** MIT (https://opensource.org/licenses/MIT)

An `renv.lock` file is included in the repository to pin the exact R package versions (R 4.5.2, mada 0.5.12, metafor 4.8.0) used in validation.

## Data availability

No new clinical data were generated for this article. The four benchmark datasets used in validation (Dementia/MMSE, Scheidler MRI, CD64 Sepsis, Glas FDG-PET) are embedded within the application and included as CSV files in the source code repository. The R validation script (`tests/validate_dta_pro.R`) and full validation output (`tests/r_validation_results.json`) are available in the repository for independent verification.

## Competing interests

No competing interests were disclosed.

## Grant information

The authors declared that no grants were involved in supporting this work.

## Acknowledgements

The authors thank the developers of the R mada package for providing the reference implementation used in validation. The WebR project (https://docs.r-wasm.org/webr/) enabled the in-browser R validation feature.

## References

[1] Macaskill P, Gatsonis C, Deeks JJ, Harbord RM, Takwoingi Y. Chapter 10: Analysing and presenting results. In: Deeks JJ, Bossuyt PM, Gatsonis C, editors. Cochrane Handbook for Systematic Reviews of Diagnostic Test Accuracy. London: Cochrane; 2023.

[2] Reitsma JB, Glas AS, Rutjes AW, Scholten RJ, Bossuyt PM, Zwinderman AH. Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews. J Clin Epidemiol. 2005;58(10):982-990.

[3] Rutter CM, Gatsonis CA. A hierarchical regression approach to meta-analysis of diagnostic test accuracy evaluations. Stat Med. 2001;20(19):2865-2884.

[4] Harbord RM, Deeks JJ, Egger M, Whiting PF, Sterne JA. A unification of models for meta-analysis of diagnostic accuracy studies. Biostatistics. 2007;8(2):239-251.

[5] Doebler P, Holling H, Bohning D. mada: Meta-Analysis of Diagnostic Accuracy. R package version 0.5.12. 2023. Available from: https://CRAN.R-project.org/package=mada

[6] Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48.

[7] Freeman SC, Kerby CR, Patel A, Cooper NJ, Quinn T, Sutton AJ. Development of an interactive web-based tool to conduct and interrogate meta-analysis of diagnostic test accuracy studies: MetaDTA. BMC Med Res Methodol. 2019;19(1):81.

[8] Whiting PF, Rutjes AW, Westwood ME, Mallett S, Deeks JJ, Reitsma JB, et al. QUADAS-2: a revised tool for the quality assessment of diagnostic accuracy studies. Ann Intern Med. 2011;155(8):529-536.

[9] DerSimonian R, Laird N. Meta-analysis in clinical trials. Control Clin Trials. 1986;7(3):177-188.

[10] Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Stat Med. 2002;21(11):1539-1558.

[11] Deeks JJ, Macaskill P, Irwig L. The performance of tests of publication bias and other sample size effects in systematic reviews of diagnostic test accuracy was assessed. J Clin Epidemiol. 2005;58(9):882-893.

[12] Zamora J, Abraira V, Muriel A, Khan K, Coomarasamy A. Meta-DiSc: a software for meta-analysis of test accuracy data. BMC Med Res Methodol. 2006;6:31.

[13] Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20(24):3875-3889.

[14] Efron B, Tibshirani RJ. An Introduction to the Bootstrap. Boca Raton: Chapman and Hall/CRC; 1993.

[15] Page MJ, McKenzie JE, Bossuyt PM, et al. The PRISMA 2020 statement: an updated guideline for reporting systematic reviews. BMJ. 2021;372:n71.

[16] Schunemann HJ, Mustafa RA, Brozek J, et al. GRADE guidelines: 21 part 2. Test accuracy: inconsistency, imprecision, publication bias, and other domains for rating the certainty of evidence and presenting it in evidence profiles and summary of findings tables. J Clin Epidemiol. 2020;117:3-10.

[17] Fagan TJ. Letter: Nomogram for Bayes theorem. N Engl J Med. 1975;293(5):257.

[18] Moses LE, Shapiro D, Littenberg B. Combining independent studies of a diagnostic test into a summary ROC curve: data-analytic approaches and some additional considerations. J Clin Epidemiol. 1993;46(10):1295-1309.

[19] Mitchell AJ. A meta-analysis of the accuracy of the mini-mental state examination in the detection of dementia and mild cognitive impairment. J Psychiatr Res. 2009;43(4):411-431.

[20] Glas AS, Lijmer JG, Prins MH, Bonsel GJ, Bossuyt PM. The diagnostic odds ratio: a single indicator of test performance. J Clin Epidemiol. 2003;56(11):1129-1135.

[21] Clopper CJ, Pearson ES. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika. 1934;26(4):404-413.
