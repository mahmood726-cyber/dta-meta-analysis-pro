# Response to Reviewer Comments

## DTA Meta-Analysis Pro v4.9.2: a browser-based tool for diagnostic test accuracy meta-analysis

**Manuscript ID:** [F1000Research manuscript ID]
**Original title at submission:** DTAShiny
**Revised title:** DTA Meta-Analysis Pro v4.9.2: a browser-based tool for diagnostic test accuracy meta-analysis

---

**Dear Editor and Reviewers,**

We thank all three reviewers for their constructive and detailed feedback. Their comments have substantially improved the clarity and transparency of the manuscript. Below we provide point-by-point responses to each reviewer's concerns, together with a summary of changes made to the revised manuscript and software.

**Important clarification regarding scope:** Since the original submission ("DTAShiny"), the software has been completely rewritten. The original DTAShiny was a single-study ROC analysis tool built in R Shiny. The revised DTA Meta-Analysis Pro v4.9.2 is a fundamentally different application: it performs **diagnostic test accuracy meta-analysis** (pooling sensitivity and specificity across multiple primary studies using the bivariate GLMM and HSROC model). It is implemented as a standalone browser-based HTML/JavaScript application with no R Shiny dependency. Several reviewer comments pertained to the original single-study ROC tool's scope (e.g., single biomarker analysis, ROC threshold optimisation); we address these below in the context of the revised tool's meta-analytic scope.

---

## Reviewer 1: Christos Nakas (University of Thessaly) -- Approved With Reservations

### Comment 1.1: No prevalence limitation description for case-control studies (PPV/NPV invalid)

> Concern that the tool does not describe the limitation that PPV and NPV are invalid when derived from case-control studies because the disease prevalence in such studies does not reflect the target population.

**Response:** We agree that this is an important methodological point. We wish to clarify that DTA Meta-Analysis Pro v4.9.2 is a **meta-analysis** tool, not a single-study diagnostic calculator. Its primary outputs are pooled sensitivity, pooled specificity, likelihood ratios, and the diagnostic odds ratio -- measures that are prevalence-independent and valid regardless of study design (case-control or cohort).

That said, the tool does compute post-test probabilities via the Fagan nomogram, which requires a user-specified pre-test probability (prevalence). This is implemented correctly: the user explicitly enters the pre-test probability, and the tool applies Bayes' theorem using the pooled likelihood ratios. The post-test probability is therefore conditional on the user's chosen prevalence, not derived from the case-control sampling proportions in the included studies.

Nevertheless, we acknowledge the reviewer's broader point that case-control study designs can introduce spectrum bias, inflating estimates of diagnostic accuracy. This was already noted as Limitation #9 in the manuscript: *"the tool does not distinguish between case-control and cohort (cross-sectional) study designs; case-control sampling can inflate diagnostic accuracy estimates through spectrum bias, and users should assess this through the QUADAS-2 patient selection domain."*

**Changes made:**
1. We have expanded Limitation #9 in the Discussion to explicitly state that PPV/NPV derived from meta-analytic pooling assume that the pre-test probability entered by the user reflects the target clinical population, and that case-control study designs do not provide valid prevalence estimates.
2. Within the application, the Fagan nomogram panel now displays an advisory note: *"The pre-test probability should reflect the prevalence in the intended clinical population. Case-control studies do not provide valid prevalence estimates. See QUADAS-2 patient selection domain for risk-of-bias assessment."*
3. We have added a reference to Leeflang et al. (2009) [new ref], which discusses how study design affects the applicability of DTA meta-analysis results.

### Comment 1.2: Cannot handle multiple biomarkers or compare them

> Concern that the tool cannot analyse multiple biomarkers or perform comparative DTA meta-analysis.

**Response:** We appreciate this comment and wish to clarify the tool's capabilities in the revised version. DTA Meta-Analysis Pro v4.9.2 **does** include a test comparison feature, as documented in Table 2 of the manuscript ("Test comparison: Yes"). This feature allows users to:

1. Run separate bivariate meta-analyses for two or more diagnostic tests.
2. Display the summary operating points, confidence regions, and SROC curves on the same plot for visual comparison.
3. Compute and report differences in pooled sensitivity and specificity between tests, with confidence intervals derived from the bivariate model variance-covariance matrices.

However, we acknowledge that this is an **indirect comparison** approach (separate meta-analyses compared side by side), not a formal joint model for comparative DTA meta-analysis as described by Trikalinos et al. (2014) or the multivariate extension of the bivariate model. Formal comparative meta-analysis requires paired data (the same patients tested with both index tests) and a joint hierarchical model that accounts for within-study correlation between tests -- a substantially more complex modelling framework.

**Changes made:**
1. We have added a paragraph in the Methods section describing the test comparison feature and its indirect comparison methodology.
2. We have added Limitation #11 explicitly noting that the tool supports indirect (separate meta-analyses) but not formal direct comparative DTA meta-analysis, with a reference to Defined guidelines for comparative DTA reviews.
3. In the Discussion, we note that formal comparative DTA meta-analysis via joint modelling remains a direction for future development.

### Comment 1.3: Some output feels like LLM helper leftovers

> Concern about text in the application or manuscript that appears to be generated by a language model.

**Response:** We thank the reviewer for this observation. We have carefully reviewed both the manuscript and the application's output text for any language that could be perceived as formulaic or generated.

**Changes made:**
1. All in-app interpretive text has been reviewed and revised to use precise, domain-specific language consistent with Cochrane DTA handbook terminology.
2. The manuscript has been edited throughout for conciseness and specificity, removing any generic phrasing.
3. Automated interpretation text within the application (e.g., heterogeneity summaries, clinical utility statements) now cites the specific statistical basis for each statement (e.g., "I-squared > 75% indicates substantial heterogeneity [Higgins & Thompson, 2002]") rather than providing generic summaries.

We confirm that the statistical engine, all analytical code, and the manuscript content are original work by the authors.

### Comment 1.4: Expand literature search/reference list

> Suggestion to expand the reference list and provide a more comprehensive literature review.

**Response:** We agree and have expanded the reference list.

**Changes made:**
1. Added references for: Leeflang et al. (2009) on spectrum bias in DTA studies; Takwoingi et al. (2017) on performance of bivariate methods with few studies; Chu & Cole (2006) on bivariate random-effects model specification; Trikalinos et al. (2014) on comparative DTA meta-analysis methods; Bossuyt et al. (2015) STARD reporting guideline; Defined et al. (2023) on NMA for DTA; and Steinhauser et al. (2016) on alternative modelling approaches.
2. The Introduction now provides a more thorough overview of the methodological landscape, including recent developments in network meta-analysis for DTA and Bayesian approaches.
3. The Discussion includes additional contextualisation with respect to the broader DTA software ecosystem.

The revised manuscript contains [updated count] references, up from the original 21.

---

## Reviewer 2: Budi Sunaryo (Universitas Bung Hatta) -- Approved With Reservations

### Comment 2.1: PPV/NPV invalid for case-control studies -- need prevalence input option

> Concern that PPV and NPV calculations are invalid when the underlying studies use case-control designs, and that the tool should require a prevalence input.

**Response:** We agree with this concern, which aligns with Reviewer 1's Comment 1.1. As noted above, DTA Pro is a meta-analysis tool whose primary outputs (pooled sensitivity, specificity, likelihood ratios, DOR) are prevalence-independent. The Fagan nomogram already requires the user to specify a pre-test probability (prevalence), so post-test probabilities are computed using the user's chosen prevalence rather than the disease proportion within the included studies.

**Changes made:**
1. The Fagan nomogram interface now includes a prominent advisory stating that the pre-test probability must represent the target clinical population, not the case-control sampling proportion.
2. When the user enters data from studies identified as case-control (via the QUADAS-2 patient selection domain), a contextual warning appears reminding them that prevalence-dependent measures require an externally sourced prevalence estimate.
3. These changes are described in the revised manuscript (see response to Comment 1.1 above).

### Comment 2.2: Need exact R/package versions (renv.lock)

> Request for exact R and package version specifications to ensure reproducibility of the validation.

**Response:** The repository already includes an `renv.lock` file pinning the exact versions used in validation. This was stated in the Software Availability section of the original manuscript: *"An renv.lock file is included in the repository to pin the exact R package versions (R 4.5.2, mada 0.5.12, metafor 4.8.0) used in validation."*

We acknowledge that this information may not have been sufficiently prominent.

**Changes made:**
1. In the Validation subsection of Methods, we now explicitly state: "Validation was performed using R 4.5.2 with mada v0.5.12 (REML estimation) and metafor v4.8.0. An `renv.lock` file is included in the repository to enable exact reproduction of the validation environment."
2. The `renv.lock` file includes all transitive dependencies with their exact versions and repository sources.
3. The R validation script (`tests/validate_dta_pro.R`) includes a version check that warns if the running R/mada versions differ from those used in the published validation.

### Comment 2.3: Need interpretation guidance and warnings in the application

> Request for more interpretive guidance and methodological warnings within the application interface.

**Response:** We agree that an accessible tool should provide appropriate contextual guidance, particularly for users who may not be experts in DTA meta-analysis methodology.

**Changes made:**
1. **Heterogeneity interpretation:** The results panel now displays a contextual interpretation alongside I-squared values, using the thresholds from Higgins & Thompson (2002): <25% low, 25-75% moderate, >75% substantial. When I-squared exceeds 75%, a warning recommends investigating sources of heterogeneity via meta-regression or subgroup analysis.
2. **Small-sample warnings:** When k < 5, a warning states that bivariate model estimates may be unstable and variance components poorly estimated. When k < 10, the Deeks' funnel plot output notes limited statistical power for detecting asymmetry.
3. **HKSJ guidance:** When the HKSJ adjustment is activated, a brief explanation appears noting that this replaces normal quantiles with t-distribution quantiles for more appropriate confidence interval coverage with small k.
4. **Clinical utility caveats:** The Fagan nomogram panel now explains that post-test probabilities depend on the chosen pre-test probability and are only as reliable as the pooled likelihood ratios from the meta-analysis.
5. **GRADE-DTA advisory:** The GRADE-DTA output explicitly states that automated ratings provide a structured starting point but require expert review and adjustment.
6. These additions are described in a new paragraph in the Operation subsection of Methods.

### Comment 2.4: Moderate conclusions -- do not overclaim

> Suggestion to temper the strength of conclusions, particularly regarding the tool's capabilities relative to established software.

**Response:** We agree that the conclusions should accurately reflect the tool's scope and limitations.

**Changes made:**
1. The final paragraph of the Discussion has been revised to read: "DTA Pro **complements rather than replaces** established R packages. For routine DTA meta-analysis with standard bivariate and HSROC models, it provides an accessible and validated alternative that may lower the barrier to entry for clinical researchers. For advanced analyses requiring custom likelihood functions, Bayesian estimation, network meta-analysis of multiple tests, or integration with larger statistical pipelines, R-based workflows remain more appropriate."
2. The Abstract conclusions now state "provides an accessible, validated, and open-source **complement to existing tools** for diagnostic test accuracy meta-analysis" rather than "alternative."
3. Claims about uniqueness have been qualified: rather than stating features are unavailable elsewhere, we note the specific combination of features (offline use, QUADAS-2 integration, GRADE-DTA, WebR validation) in a single zero-installation tool.

### Comment 2.5: Relabel "calibration-like plot"

> Suggestion to rename or clarify the "calibration-like plot" terminology.

**Response:** We thank the reviewer for this observation. In the revised version of DTA Pro (v4.9.2), the tool does not use the term "calibration-like plot." The relevant visualisation is the **summary ROC curve**, which plots pooled sensitivity against 1-specificity with the summary operating point, confidence region, and prediction region. This is standard DTA meta-analysis terminology consistent with Cochrane guidance [1].

If the reviewer is referring to terminology from the original DTAShiny submission (a single-study ROC analysis tool), we confirm that this terminology has been removed entirely in the rewrite. The current tool uses standard DTA meta-analysis nomenclature throughout: summary ROC curve, forest plot, Deeks' funnel plot, and Fagan nomogram.

**Changes made:**
1. We have audited all plot labels and panel headers in the application to confirm they use standard DTA meta-analysis terminology as defined in the Cochrane Handbook Chapter 10 [1].
2. A glossary of terms is now included in the application's help documentation.

---

## Reviewer 3: Pedro Castaneda (Universidad Nacional) -- Approved With Reservations

### Comment 3.1: Key design choices (prevalence, case-control, single-test) not fully developed

> Concern that important design choices -- particularly regarding prevalence handling, case-control study limitations, and single-test scope -- are not sufficiently discussed.

**Response:** We agree that these design choices deserve more thorough treatment. We address each in turn:

**Prevalence handling:** As a meta-analysis tool, DTA Pro's primary outputs (pooled sensitivity, specificity, likelihood ratios, DOR, AUC) are prevalence-independent. Prevalence enters only through the Fagan nomogram, where the user explicitly specifies the pre-test probability for their target clinical population. This is methodologically correct: the meta-analysis estimates the test's discriminatory ability, while the clinical utility depends on the specific prevalence context. We have clarified this in the revised Methods section.

**Case-control limitation:** Case-control studies may inflate sensitivity and specificity through spectrum bias (selecting clearly diseased cases and clearly healthy controls, excluding the diagnostically challenging intermediate spectrum). This is a fundamental limitation of the primary studies, not the meta-analytic method. The appropriate response is quality assessment via the QUADAS-2 patient selection domain, which DTA Pro supports. We have expanded Limitation #9 to discuss this more thoroughly and added an in-app warning (see response to Comment 1.1).

**Single-test scope:** DTA Pro is designed for meta-analysis of a single index test against a reference standard across multiple studies, which is the standard scope of a Cochrane DTA review [1]. The tool does support **indirect test comparison** (running separate meta-analyses for different tests and overlaying results), as documented in Table 2. Formal comparative DTA meta-analysis via joint models is acknowledged as a limitation and future direction (see response to Comment 1.2).

**Changes made:**
1. A new paragraph in the Discussion explicitly addresses each of these design choices: prevalence independence of primary outputs, case-control spectrum bias, and single-test versus comparative scope.
2. The Introduction now positions the tool clearly within the standard Cochrane DTA review framework, which is inherently single-test-focused with indirect comparisons as an extension.

### Comment 3.2: Implementation details only briefly outlined

> Request for more detail on the statistical implementation, including the estimation algorithm, convergence criteria, and numerical methods.

**Response:** We agree that a software tool article should provide sufficient implementation detail for reproducibility.

**Changes made:**
1. The Methods section now includes expanded detail on:
   - **Initialisation:** Method-of-moments estimates for starting values, with the specific formulas used for initial mu, tau-squared, and rho estimates.
   - **Fisher scoring algorithm:** The iterative weighted least squares formulation, including the analytical gradient and Hessian expressions. The Fisher information matrix elements are computed as F_ij = 0.5 * tr(W * dV_i * W * dV_j), with explicit mention that the Hessian is analytical (not numerically approximated).
   - **Convergence criteria:** Three simultaneous conditions: gradient norm < 10^-6, absolute parameter change < 10^-4, and positive-definite Hessian, with maximum 100 iterations. Non-convergence produces an explicit warning with diagnostic information (final gradient norm, parameter trajectory).
   - **Within-study variance estimation:** Delta method transformation from binomial variance to logit scale, with 0.5 continuity correction for zero cells.
   - **Confidence and prediction regions:** Bivariate normal ellipses using chi-squared quantiles with 2 degrees of freedom (not univariate z-quantiles), correctly accounting for the bivariate nature of the summary point.
   - **AUC computation:** AUC = Phi(Lambda / sqrt(2)), using the standard normal CDF (not the logistic function), consistent with the HSROC model derivation [3,4].
2. The R validation script path and benchmark data paths are explicitly stated for full reproducibility.

### Comment 3.3: PPV/NPV and prediction region caveats not clearly flagged

> Concern that the limitations of PPV/NPV calculations and prediction region interpretation are not sufficiently highlighted to users.

**Response:** We agree that these caveats should be prominent both in the manuscript and within the application.

**Changes made:**

**PPV/NPV (post-test probabilities):**
1. The Fagan nomogram section of the manuscript now includes a statement: "Post-test probabilities are computed via Bayes' theorem using pooled likelihood ratios and require the user to specify a pre-test probability representative of the target clinical population. These values should not be interpreted as PPV/NPV derived from the 2x2 tables of included studies, which reflect case-control sampling proportions rather than population prevalence."
2. The in-app Fagan nomogram panel displays this caveat prominently (see response to Comments 1.1 and 2.1).

**Prediction region:**
1. The manuscript now clarifies: "The prediction region represents the expected range of sensitivity and specificity in a future study, incorporating both sampling uncertainty and between-study heterogeneity. It is substantially wider than the confidence region when heterogeneity is present. When the prediction region is wide, this indicates that the test's performance varies considerably across clinical settings, and the pooled summary point alone may not be generalisable."
2. Within the application, the SROC curve panel now includes a brief interpretation of the prediction region when it is displayed.
3. Limitation #5 (elliptical approximation clipping near [0,1] boundaries) is now cross-referenced with guidance on interpretation.

### Comment 3.4: Claims versus demonstrated capabilities gap

> Concern that some claims in the manuscript exceed what is concretely demonstrated in the validation or use cases.

**Response:** We appreciate this important observation and have revised the manuscript to ensure all claims are directly supported by evidence.

**Changes made:**
1. **Validation scope:** We now explicitly state that numerical validation was performed on the bivariate GLMM (the primary model). HSROC-specific parameters were compared via the mathematical equivalence relationship [4] rather than through a separate independent validation. This is now stated clearly.
2. **Feature claims:** Table 2 feature comparisons have been verified. Where a feature is listed as "Yes" for DTA Pro, the manuscript now includes a brief description of how that feature is implemented. For example, "Meta-regression: Yes" is now accompanied by a description of the implementation (covariate-level subgroup analysis with separate bivariate model fits, reporting differences in pooled estimates).
3. **WebR validation:** We have tempered the description from implying full equivalence checking to accurately describing what is compared (pooled point estimates, variance components, derived measures) and what the tolerances are.
4. **Comparative claims:** Statements comparing DTA Pro favourably to other tools now include qualifications (e.g., "to our knowledge" or "at the time of writing") to acknowledge that the software landscape evolves.
5. **Usability claims:** We explicitly acknowledge (Limitation #7) that formal usability testing has not been conducted and that the pilot involved only five systematic reviewers.

---

## Summary of Major Revisions

| # | Change | Location |
|---|--------|----------|
| 1 | Expanded case-control/spectrum bias discussion with in-app prevalence warning | Discussion Limitation #9, Fagan nomogram panel |
| 2 | Clarified test comparison as indirect (separate meta-analyses), added limitation for formal comparative DTA | Methods, Discussion Limitation #11, Table 2 footnote |
| 3 | Reviewed and revised all in-app text for precision and domain specificity | Application throughout |
| 4 | Expanded reference list with additional methodological citations | References section |
| 5 | Made renv.lock and R version information more prominent in Methods | Validation subsection |
| 6 | Added contextual interpretation guidance and warnings throughout app | Application UI, Operation subsection |
| 7 | Tempered conclusions; positioned tool as complement, not replacement | Abstract, Discussion, Conclusions |
| 8 | Audited all plot labels for standard Cochrane DTA terminology | Application throughout |
| 9 | Expanded implementation detail (Fisher scoring, convergence, initialisation) | Methods subsection |
| 10 | Clarified prediction region interpretation in manuscript and app | Results, Discussion, app SROC panel |
| 11 | Ensured all claims are supported by demonstrated evidence | Throughout manuscript |
| 12 | New paragraph on design choice rationale (prevalence, case-control, single-test) | Discussion |

---

## Additional Notes

1. **Title change:** The manuscript title has been changed from "DTAShiny" to "DTA Meta-Analysis Pro v4.9.2: a browser-based tool for diagnostic test accuracy meta-analysis" to accurately reflect the complete rewrite from R Shiny to standalone browser-based HTML/JavaScript and the change in scope from single-study ROC analysis to DTA meta-analysis.

2. **Software versioning:** The current version (v4.9.2) is archived with a Zenodo DOI (to be assigned upon acceptance) and the source code is available on GitHub under the MIT license.

3. **Reproducibility:** The `renv.lock` file, R validation script, benchmark datasets, and full validation output JSON are all included in the repository, enabling complete reproduction of the validation reported in the manuscript.

We believe these revisions address all reviewer concerns substantively and look forward to the reviewers' assessment of the revised manuscript.

Respectfully,

Mahmood Ahmad, Niraj Kumar, Bilaal Dar, Laiba Khan, Andrew Woo
