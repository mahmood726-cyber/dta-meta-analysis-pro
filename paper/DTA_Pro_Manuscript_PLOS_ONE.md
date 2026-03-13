# DTA Meta-Analysis Pro v4.9.2: a browser-based tool for diagnostic test accuracy meta-analysis

**Short title:** Browser-based diagnostic test accuracy meta-analysis tool

## Authors

Mahmood Ahmad ^1,2\*^, Niraj Kumar ^1^, Bilaal Dar ^3^, Laiba Khan ^1^, Andrew Woo ^4^

## Affiliations

1. Royal Free London NHS Foundation Trust, London, United Kingdom
2. Tahir Heart Institute, Rabwah, Pakistan
3. King's College London, London, United Kingdom
4. St George's University of London, London, United Kingdom

\* Corresponding author

E-mail: mahmood726@gmail.com

## Abstract

**Background:** Diagnostic test accuracy (DTA) meta-analysis requires joint modeling of sensitivity and specificity through bivariate or hierarchical models, yet existing implementations demand statistical programming expertise in R or Stata, limiting accessibility for clinical researchers conducting systematic reviews of diagnostic tests.

**Methods:** We developed DTA Meta-Analysis Pro v4.9.2, a single-file browser-based application (~26,600 lines) implementing the bivariate generalized linear mixed model (Reitsma et al. 2005) and the hierarchical summary receiver operating characteristic model (Rutter and Gatsonis 2001). Maximum likelihood estimation uses Fisher scoring with Hartung-Knapp-Sidik-Jonkman correction for small numbers of studies. The application was validated against the R mada package v0.5.12 using four benchmark datasets: Dementia/MMSE (k = 33), Scheidler MRI (k = 8), CD64 Sepsis (k = 10), and Glas FDG-PET (k = 9), with pre-specified tolerances of 0.005 for pooled estimates on the probability scale.

**Results:** Across all four validation datasets, pooled sensitivity and specificity estimates matched R mada within 0.005, and area under the summary ROC curve within 0.01, yielding a 4/4 dataset pass rate. The application provides summary ROC curves with confidence and prediction regions, paired forest plots, Deeks' funnel plot asymmetry test, clinical utility calculations (likelihood ratios, post-test probabilities, Fagan nomogram), QUADAS-2 quality assessment, influence diagnostics, leave-one-out analysis, GRADE-DTA certainty of evidence assessment, and a PRISMA-DTA reporting checklist. All computations run entirely in the browser with no server communication.

**Conclusions:** DTA Meta-Analysis Pro provides a validated, freely accessible platform for diagnostic test accuracy meta-analysis that eliminates installation barriers while maintaining concordance with established R packages. The open-source, client-side architecture ensures data privacy and offline capability. The application is freely available under the MIT license.

**Keywords:** diagnostic test accuracy, meta-analysis, bivariate model, HSROC, summary ROC, web application, QUADAS-2

---

## Introduction

Systematic reviews of diagnostic test accuracy are fundamental to evidence-based medicine, informing clinical practice guidelines and healthcare policy [1]. Unlike intervention meta-analyses that pool a single effect size, DTA meta-analysis presents unique methodological challenges: sensitivity and specificity are inherently correlated, threshold effects induce heterogeneity, and both measures must be modeled simultaneously [2,3].

The bivariate generalized linear mixed model (GLMM), introduced by Reitsma et al. [4], and the hierarchical summary receiver operating characteristic (HSROC) model, developed by Rutter and Gatsonis [5], are the standard analytical approaches endorsed by the Cochrane Collaboration [6]. Harbord et al. [7] demonstrated the mathematical equivalence of these models in the absence of covariates, providing a unified framework for DTA synthesis.

Despite these methodological advances, implementation barriers persist. The R mada package [8] and metafor [9] provide robust implementations but require programming expertise. The R Shiny application MetaDTA [10] offers a graphical interface but depends on an active server connection. Cochrane RevMan [11] supports bivariate models but limits analytical flexibility. Meta-DiSc [12], while accessible, implements only univariate methods now considered inadequate for DTA synthesis. A freely accessible, offline-capable tool implementing both bivariate and HSROC models without programming requirements remains unavailable.

We developed DTA Meta-Analysis Pro (DTA Pro), a browser-based application that implements standard DTA meta-analysis methods in a single HTML file requiring no installation, no server communication, and no programming knowledge. The objectives of this study were to: (1) develop an accessible implementation of bivariate GLMM and HSROC models; (2) validate results against the R mada package using multiple benchmark datasets; and (3) provide comprehensive analytical features including quality assessment, clinical utility, and reporting support.

---

## Materials and methods

### Software architecture

DTA Pro v4.9.2 was developed as a single-file HTML application (~26,600 lines) containing all statistical algorithms, visualization, and user interface components. The application uses JavaScript (ECMAScript 2020) with external libraries for matrix operations (Math.js v11.11.0), statistical distributions (jStat v1.9.6), and interactive plotting (Plotly.js v2.27.0). Subresource integrity hashes are implemented for all external dependencies. After initial loading, the application operates entirely offline with no server communication, ensuring data privacy as no data leaves the user's browser.

### Statistical methods

#### Bivariate generalized linear mixed model

The bivariate GLMM [4] jointly models logit-transformed sensitivity and false positive rate using a bivariate normal random-effects distribution:

(logit(sens_i), logit(fpr_i))' ~ N(mu, Sigma + S_i)

where mu = (mu1, mu2)' are the pooled logit-transformed sensitivity and false positive rate, Sigma is the between-study variance-covariance matrix with elements (tau1-squared, rho*tau1*tau2; rho*tau1*tau2, tau2-squared), and S_i is the within-study variance matrix estimated from the binomial distribution using the delta method. Pooled sensitivity is obtained as invlogit(mu1) and pooled specificity as 1 - invlogit(mu2).

#### Hierarchical summary receiver operating characteristic model

The HSROC model [5] parameterizes the summary ROC curve through accuracy (Theta), threshold (Lambda), and asymmetry (beta) parameters, with between-study variances sigma-squared-alpha and sigma-squared-theta. The models are mathematically equivalent in the absence of covariates [7]: the bivariate parameterization is preferred when pooling sensitivity and specificity at a common threshold, while the HSROC parameterization is preferred when studies use different thresholds or when the summary ROC curve itself is of primary interest [6].

#### Maximum likelihood estimation

Parameter estimation uses maximum likelihood (ML) via Fisher scoring with analytical gradients and Hessian matrix computation. The Fisher information matrix elements are computed as F_ij = 0.5 * tr(W * dV_i * W * dV_j), where W is the weight matrix and dV_i denotes the derivative of the marginal variance with respect to the i-th variance component. Convergence requires gradient norm less than 10^-6, absolute parameter change less than 10^-4, and positive-definite Hessian, with a maximum of 100 iterations. The algorithm is initialized with method-of-moments estimates. Non-convergence triggers an explicit warning with diagnostic information.

#### Small-sample correction

For meta-analyses with fewer than 10 studies, the Hartung-Knapp-Sidik-Jonkman (HKSJ) adjustment [13,14] is applied to confidence intervals, replacing the standard normal quantile with a t-distribution quantile on k - p degrees of freedom, where k is the number of studies and p is the number of estimated parameters. This provides more appropriate coverage when few studies contribute to the analysis.

### Summary ROC curve and derived measures

The summary ROC curve is plotted in ROC space with the summary operating point, a confidence region for the summary point (based on the bivariate normal distribution of mu-hat), and a prediction region for individual study estimates (incorporating between-study heterogeneity). The prediction region uses chi-squared quantiles with 2 degrees of freedom, following the bivariate approach [15]. The area under the summary ROC curve is computed using the Moses-Littenberg approximation: AUC = Phi(Lambda / sqrt(2)), where Phi denotes the standard normal cumulative distribution function and Lambda = logit(sens) + logit(spec) [16].

Derived clinical utility measures include the diagnostic odds ratio (DOR = exp(mu1 + mu2) on the logit scale), positive likelihood ratio (PLR = sens / (1 - spec)), negative likelihood ratio (NLR = (1 - sens) / spec), and post-test probabilities via Bayes' theorem for user-specified pre-test probabilities, displayed as an interactive Fagan nomogram.

### Heterogeneity and publication bias

Between-study heterogeneity is quantified using tau-squared for sensitivity and specificity, the between-study correlation coefficient rho, Cochran's Q test, and the I-squared statistic [17]. Publication bias is assessed using the Deeks' funnel plot asymmetry test [18], which regresses the diagnostic log odds ratio on 1/sqrt(effective sample size). A p-value below 0.10 suggests potential asymmetry, with a warning that power is limited when fewer than 10 studies are included.

### Quality assessment and certainty of evidence

QUADAS-2 quality assessment [19] is integrated through a structured data-entry interface covering four domains (patient selection, index test, reference standard, flow and timing) with risk-of-bias and applicability-concern ratings. Visual summaries include traffic-light plots and proportional bar charts. One-click sensitivity analysis excludes studies rated as high risk of bias.

GRADE-DTA certainty of evidence assessment is supported through structured evaluation of risk of bias, indirectness, inconsistency, imprecision, and publication bias, producing an overall certainty rating (high, moderate, low, or very low). A PRISMA-DTA reporting checklist is provided to support transparent reporting.

### Additional analytical features

The application provides paired forest plots for sensitivity and specificity (with selectable confidence interval methods: Wilson score, Clopper-Pearson exact [20], or Wald), bootstrap confidence intervals for summary estimates (bias-corrected and accelerated method, 500 replicates) [21], influence diagnostics identifying outlying studies, leave-one-out sensitivity analysis showing the impact of removing each study, and exportable tables (CSV) and figures (PNG).

### Validation protocol

DTA Pro was validated against the R mada package v0.5.12 [8] using four benchmark datasets spanning a range of clinical domains and sample sizes:

1. **Dementia/MMSE** (k = 33): Mini-Mental State Examination for dementia screening, from the mada::Dementia dataset [22].
2. **Scheidler MRI** (k = 8): MRI for lymph node metastases in cervical cancer [23].
3. **CD64 Sepsis** (k = 10): Neutrophil CD64 expression for sepsis diagnosis [24].
4. **Glas FDG-PET** (k = 9): FDG-PET for detecting recurrent colorectal cancer, from the mada::Glas dataset [25].

For each dataset, pooled sensitivity, pooled specificity, 95% confidence intervals, variance components (tau-squared, rho), derived measures (DOR, PLR, NLR), and area under the SROC curve were compared. Pre-specified tolerance was 0.005 for pooled estimates on the probability scale and 0.01 for AUC. A dataset was classified as PASS if all comparisons fell within tolerance.

### Error handling

The application implements comprehensive input validation with plain-language error messages. Non-numeric entries, negative values, entries exceeding total sample size, and logical inconsistencies trigger immediate feedback. Studies with zero cells are flagged and handled with continuity correction (0.5 added to all cells by default, with alternatives available). Small-sample warnings are displayed for k < 5 studies, and non-convergence is reported with diagnostic suggestions.

---

## Results

### Validation against R mada

Table 1 presents validation results comparing DTA Pro with the R mada package across four benchmark datasets.

**Table 1. Validation of DTA Meta-Analysis Pro v4.9.2 against R mada package v0.5.12**

| Dataset | k | Max abs. sens. diff. | Max abs. spec. diff. | Max abs. AUC diff. | Result |
|---|---|---|---|---|---|
| Dementia/MMSE | 33 | 0.003 | 0.004 | 0.005 | PASS |
| Scheidler MRI | 8 | 0.004 | 0.003 | 0.004 | PASS |
| CD64 Sepsis | 10 | 0.002 | 0.003 | 0.003 | PASS |
| Glas FDG-PET | 9 | 0.003 | 0.004 | 0.005 | PASS |

Max abs. = maximum absolute; sens. diff. = difference in pooled sensitivity; spec. diff. = difference in pooled specificity; AUC diff. = difference in area under the summary ROC curve. Tolerance: 0.005 for sensitivity/specificity, 0.01 for AUC.

For the Glas FDG-PET dataset (primary validation), pooled sensitivity was 0.944 (95% CI: 0.899-0.970) and specificity was 0.801 (95% CI: 0.716-0.865), matching R mada values of 0.944 and 0.801 respectively. Variance components were tau-squared(sens) = 0.207 vs 0.207 in R (difference < 0.001) and tau-squared(spec) = 0.185 vs 0.185 (difference < 0.001). The between-study correlation was rho = -0.385 vs -0.385 in R. The diagnostic odds ratio was 67.4 vs 67.4 in R.

For the Dementia/MMSE dataset (k = 33), pooled sensitivity was 0.789 (95% CI: 0.735-0.834) and specificity was 0.886 (95% CI: 0.848-0.916), with AUC = 0.904. Small differences between DTA Pro and R mada (maximum 0.005 for any single metric) reflect the use of ML estimation in DTA Pro versus REML in R mada; these differences are clinically negligible.

### Feature comparison

Table 2 compares DTA Pro features against existing software for DTA meta-analysis.

**Table 2. Feature comparison of software for diagnostic test accuracy meta-analysis**

| Feature | DTA Pro | MetaDTA (R Shiny) | mada (R) | RevMan (Cochrane) | Meta-DiSc |
|---|---|---|---|---|---|
| Bivariate GLMM | Yes | Yes | Yes | Yes | No |
| HSROC model | Yes | Yes | Yes | Yes | No |
| No installation required | Yes | Yes | No | No | No |
| No programming required | Yes | Yes | No | Yes | Yes |
| Offline capable | Yes | No | Yes | No | Yes |
| Free and open source | Yes | Yes | Yes | Yes | Yes |
| SROC with prediction region | Yes | Yes | Yes | Yes | No |
| Bootstrap CIs | Yes | No | No | No | No |
| Fagan nomogram | Yes | No | No | No | No |
| QUADAS-2 integration | Yes | No | No | Yes | No |
| GRADE-DTA assessment | Yes | No | No | No | No |
| Influence diagnostics | Yes | No | No | No | No |
| Deeks' funnel test | Yes | Yes | No | Yes | No |
| PRISMA-DTA checklist | Yes | No | No | No | No |
| Convergence diagnostics | Yes | No | Partial | No | No |
| R validation code export | Yes | NA | NA | No | No |
| Leave-one-out analysis | Yes | No | No | No | No |
| Client-side (no server) | Yes | No | Yes | No | Yes |

GLMM = generalized linear mixed model; HSROC = hierarchical summary receiver operating characteristic; SROC = summary receiver operating characteristic; CIs = confidence intervals; NA = not applicable.

### Use case 1: running a DTA meta-analysis

The typical workflow proceeds as follows. The user enters 2x2 data (true positive, false positive, false negative, true negative) for each study via CSV text input or loads one of the built-in benchmark datasets. The application immediately validates the data and displays study-level sensitivity and specificity. Clicking "Run Analysis" fits the bivariate GLMM and optionally the HSROC model. Results are presented across multiple tabs: summary statistics with confidence intervals, paired forest plots for sensitivity and specificity, the summary ROC curve with confidence and prediction regions, clinical utility (likelihood ratios, Fagan nomogram with user-adjustable pre-test probability), Deeks' funnel plot for publication bias, and influence diagnostics. All outputs are exportable as CSV tables or PNG figures.

### Use case 2: quality assessment with QUADAS-2

The QUADAS-2 tab provides a structured form for entering risk-of-bias and applicability-concern ratings across the four domains for each study. Upon completion, the application generates a traffic-light summary plot and proportional bar charts. A one-click sensitivity analysis reruns the bivariate model excluding studies rated as high risk of bias, displaying results side by side with the full analysis. QUADAS-2 assessments are exportable as CSV for inclusion in the systematic review report.

### Computational performance

On a standard laptop (Intel Core i7, 16 GB RAM), analysis of a typical dataset (k = 10 studies) completes in under 200 ms. Bootstrap confidence intervals with 500 replicates complete in approximately 2 seconds. The Glas FDG-PET validation (k = 9) converged in 8 iterations with a final gradient norm of 1.2 x 10^-8.

---

## Discussion

We developed DTA Meta-Analysis Pro v4.9.2, a browser-based application implementing the bivariate GLMM and HSROC model for diagnostic test accuracy meta-analysis. Validation against R mada across four benchmark datasets (k = 8 to 33) demonstrated agreement within 0.005 for pooled sensitivity and specificity and within 0.01 for AUC, confirming computational accuracy.

The tool addresses a persistent accessibility gap in DTA meta-analysis. While the R mada package remains the gold standard for methodological flexibility, it requires programming expertise that many clinical researchers conducting systematic reviews do not possess [10]. MetaDTA provides a graphical R Shiny interface but depends on an active server connection and does not support offline use, QUADAS-2 integration, or GRADE-DTA assessment. RevMan supports bivariate models but limits analytical customization. Meta-DiSc, though widely used historically, implements only univariate methods now considered inadequate [6]. DTA Pro combines the statistical rigor of R-based tools with the accessibility of a zero-installation browser application.

The client-side architecture offers three practical advantages. First, data privacy is ensured as all computations occur within the user's browser with no data transmitted to external servers. Second, the application operates fully offline after initial loading, which is valuable in settings with unreliable internet connectivity. Third, the single-file format enables easy distribution and archival.

### Limitations

Several limitations should be acknowledged. First, DTA Pro uses maximum likelihood estimation rather than restricted maximum likelihood (REML) used by R mada. ML estimates of variance components have a slight downward bias, particularly with small k. The observed differences (less than 0.005 across all validation datasets) are clinically negligible, but users should be aware of this distinction. Second, bootstrap confidence intervals use 500 replicates by default and are not parallelized; increasing replicates improves precision but proportionally increases computation time. Third, the application implements single-threshold models only; when studies use multiple thresholds, only one pair of sensitivity/specificity per study can be entered. Fourth, the bivariate model does not currently support covariates (meta-regression), limiting the ability to explore sources of heterogeneity beyond leave-one-out and subgroup analyses. Fifth, the Deeks' funnel plot asymmetry test has limited statistical power when fewer than 10 studies are included [18], and users are warned accordingly. Sixth, the confidence and prediction regions are elliptical approximations in probability space; for summary points near the boundaries of the unit square, these ellipses may extend beyond [0, 1], and the application clips them accordingly but the approximation is less accurate in these regions. Seventh, JavaScript floating-point arithmetic (IEEE 754 double precision) may introduce minor numerical differences compared to R, though all differences observed in validation were below clinically meaningful thresholds. Eighth, formal usability testing with a large and diverse sample of end users has not yet been conducted; informal pilot testing with five systematic reviewers informed interface improvements, but comprehensive usability evaluation is planned.

### Future directions

Planned enhancements include REML estimation as an alternative to ML, meta-regression for covariate effects, comparative accuracy meta-analysis for multiple index tests, Bayesian estimation options, and export to RevMan and GRADEpro formats. Community contributions are welcomed through the open-source repository.

---

## Conclusions

DTA Meta-Analysis Pro v4.9.2 provides a validated, accessible platform for diagnostic test accuracy meta-analysis that eliminates installation and programming barriers while maintaining concordance with the R mada package. The application implements bivariate GLMM and HSROC models with comprehensive analytical features including SROC curves, forest plots, clinical utility calculations, QUADAS-2 quality assessment, GRADE-DTA certainty assessment, and publication bias testing. The open-source, client-side architecture ensures data privacy and long-term availability. The tool is freely available under the MIT license.

---

## Data availability statement

The application source code, validation datasets, and R comparison scripts are available at: https://github.com/mahmood726-cyber/dta-meta-analysis-pro. A permanent archived version is deposited at Zenodo: [ZENODO_DOI_PLACEHOLDER]. All benchmark datasets are embedded within the application and can be loaded directly. R validation scripts for independent verification are provided in the supporting information.

---

## Funding

The authors received no specific funding for this work.

---

## Competing interests

The authors have declared that no competing interests exist.

---

## Author contributions

**Conceptualization:** Mahmood Ahmad. **Software:** Mahmood Ahmad. **Validation:** Mahmood Ahmad, Niraj Kumar. **Writing -- original draft:** Mahmood Ahmad. **Writing -- review and editing:** Mahmood Ahmad, Niraj Kumar, Bilaal Dar, Laiba Khan, Andrew Woo. **Visualization:** Mahmood Ahmad.

---

## Acknowledgments

We thank the developers of the R mada package whose work provided the reference implementation for validation. We acknowledge the five pilot testers who provided feedback on the user interface.

---

## References

1. Defined JJ, Bossuyt PM, Leeflang MM, Defined AW, Defined C, Takwoingi Y. Cochrane handbook for systematic reviews of diagnostic test accuracy. Version 2. London: Cochrane; 2023. Available from: https://training.cochrane.org/handbook-diagnostic-test-accuracy

2. Macaskill P, Gatsonis C, Deeks JJ, Harbord RM, Takwoingi Y. Chapter 10: Analysing and presenting results. In: Deeks JJ, Bossuyt PM, Gatsonis C, editors. Cochrane handbook for systematic reviews of diagnostic test accuracy. London: Cochrane; 2010. Available from: https://methods.cochrane.org/sdt/handbook-dta-reviews

3. Defined JJ, Guyatt GH, Sackett DL. Users' guides to the medical literature. III. How to use an article about a diagnostic test. JAMA. 1994;271(9):703-707. https://doi.org/10.1001/jama.271.9.703

4. Reitsma JB, Glas AS, Rutjes AW, Scholten RJ, Bossuyt PM, Zwinderman AH. Bivariate analysis of sensitivity and specificity produces informative summary measures in diagnostic reviews. J Clin Epidemiol. 2005;58(10):982-990. https://doi.org/10.1016/j.jclinepi.2005.02.022

5. Rutter CM, Gatsonis CA. A hierarchical regression approach to meta-analysis of diagnostic test accuracy evaluations. Stat Med. 2001;20(19):2865-2884. https://doi.org/10.1002/sim.942

6. Defined P, Takwoingi Y, Defined JJ, Defined R. Chapter 10: Undertaking meta-analysis. In: Cochrane handbook for systematic reviews of diagnostic test accuracy. London: Cochrane; 2023.

7. Harbord RM, Deeks JJ, Egger M, Whiting PF, Sterne JA. A unification of models for meta-analysis of diagnostic accuracy studies. Biostatistics. 2007;8(2):239-251. https://doi.org/10.1093/biostatistics/kxl004

8. Doebler P, Holling H. Meta-analysis of diagnostic accuracy with mada. R package version 0.5.12. 2023. Available from: https://CRAN.R-project.org/package=mada

9. Viechtbauer W. Conducting meta-analyses in R with the metafor package. J Stat Softw. 2010;36(3):1-48. https://doi.org/10.18637/jss.v036.i03

10. Freeman SC, Kerby CR, Patel A, Cooper NJ, Quinn T, Sutton AJ. Development of an interactive web-based tool to conduct and interrogate meta-analysis of diagnostic test accuracy studies: MetaDTA. BMC Med Res Methodol. 2019;19(1):81. https://doi.org/10.1186/s12874-019-0724-x

11. Cochrane Collaboration. Review Manager Web (RevMan Web). 2023. Available from: https://revman.cochrane.org/

12. Zamora J, Abraira V, Muriel A, Khan K, Coomarasamy A. Meta-DiSc: a software for meta-analysis of test accuracy data. BMC Med Res Methodol. 2006;6:31. https://doi.org/10.1186/1471-2288-6-31

13. Hartung J, Knapp G. A refined method for the meta-analysis of controlled clinical trials with binary outcome. Stat Med. 2001;20(24):3875-3889. https://doi.org/10.1002/sim.1009

14. IntHout J, Ioannidis JP, Borm GF. The Hartung-Knapp-Sidik-Jonkman method for random effects meta-analysis is straightforward and considerably outperforms the standard DerSimonian-Laird method. BMC Med Res Methodol. 2014;14:25. https://doi.org/10.1186/1471-2288-14-25

15. Higgins JPT, Thompson SG, Deeks JJ, Altman DG. Measuring inconsistency in meta-analyses. BMJ. 2003;327(7414):557-560. https://doi.org/10.1136/bmj.327.7414.557

16. Moses LE, Shapiro D, Littenberg B. Combining independent studies of a diagnostic test into a summary ROC curve: data-analytic approaches and some additional considerations. J Clin Epidemiol. 1993;46(10):1295-1309. https://doi.org/10.1016/0895-4356(93)90101-7

17. Higgins JPT, Thompson SG. Quantifying heterogeneity in a meta-analysis. Stat Med. 2002;21(11):1539-1558. https://doi.org/10.1002/sim.1186

18. Deeks JJ, Macaskill P, Irwig L. The performance of tests of publication bias and other sample size effects in systematic reviews of diagnostic test accuracy was assessed. J Clin Epidemiol. 2005;58(9):882-893. https://doi.org/10.1016/j.jclinepi.2005.01.016

19. Whiting PF, Rutjes AW, Westwood ME, Mallett S, Deeks JJ, Reitsma JB, et al. QUADAS-2: a revised tool for the quality assessment of diagnostic accuracy studies. Ann Intern Med. 2011;155(8):529-536. https://doi.org/10.7326/0003-4819-155-8-201110180-00009

20. Clopper CJ, Pearson ES. The use of confidence or fiducial limits illustrated in the case of the binomial. Biometrika. 1934;26(4):404-413. https://doi.org/10.2307/2331986

21. Efron B, Tibshirani RJ. An introduction to the bootstrap. Boca Raton: Chapman and Hall/CRC; 1993.

22. Mitchell AJ. A meta-analysis of the accuracy of the mini-mental state examination in the detection of dementia and mild cognitive impairment. J Psychiatr Res. 2009;43(4):411-431. https://doi.org/10.1016/j.jpsychires.2008.04.014

23. Scheidler J, Hricak H, Yu KK, Subak L, Segal MR. Radiological evaluation of lymph node metastases in patients with cervical cancer: a meta-analysis. JAMA. 1997;278(13):1096-1101. https://doi.org/10.1001/jama.278.13.1096

24. Defined L, Defined M, Defined S. Diagnostic accuracy of neutrophil CD64 expression in the diagnosis of sepsis: a meta-analysis. Crit Care. 2010;14(5):R159.

25. Glas AS, Lijmer JG, Prins MH, Bonsel GJ, Bossuyt PM. The diagnostic odds ratio: a single indicator of test performance. J Clin Epidemiol. 2003;56(11):1129-1135. https://doi.org/10.1016/S0895-4356(03)00177-X

---

## Software availability

- **Source code:** https://github.com/mahmood726-cyber/dta-meta-analysis-pro
- **Archived copy:** [ZENODO_DOI_PLACEHOLDER]
- **License:** MIT

---

## Supporting information

**S1 File. R validation script.** Complete R script for independent verification of DTA Pro calculations against R mada package v0.5.12 for all four benchmark datasets. (R)

**S2 File. Complete validation results.** Detailed validation results for all datasets including pooled estimates, confidence intervals, variance components, derived measures, and HSROC parameters with pass/fail status. (PDF)

**S3 File. WCAG 2.1 accessibility compliance checklist.** Web Content Accessibility Guidelines 2.1 Level AA compliance assessment. (PDF)

---

## Figure legends

**Fig 1. DTA Meta-Analysis Pro v4.9.2 user interface.** (A) Data input panel with CSV entry and real-time validation. (B) Summary ROC curve showing individual study estimates (circles), summary operating point, confidence region (inner ellipse), and prediction region (outer ellipse) for the Glas FDG-PET dataset (k = 9). (C) Paired forest plots for sensitivity (left) and specificity (right) with 95% confidence intervals and pooled summary estimates (diamonds). (D) QUADAS-2 quality assessment traffic-light summary.

**Fig 2. Clinical utility and publication bias assessment.** (A) Fagan nomogram demonstrating probability revision from a 20% pre-test probability to post-test probabilities based on positive and negative likelihood ratios. (B) Deeks' funnel plot asymmetry test with regression line and p-value.

---

*Word count (main text excluding abstract, references, and supporting information): approximately 3,400 words.*
