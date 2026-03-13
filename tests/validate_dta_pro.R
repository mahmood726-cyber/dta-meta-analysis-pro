#!/usr/bin/env Rscript
# =============================================================================
# DTA Meta-Analysis Pro -- R Validation Script
# =============================================================================
#
# Purpose:
#   Validates the bivariate GLMM (Reitsma model) calculations in DTA Pro
#   against the R `mada` package (reitsma() function). This script is part
#   of the reproducibility materials for a PLOS ONE submission.
#
# Methodology:
#   The bivariate generalized linear mixed model (Reitsma et al. 2005) jointly
#   models sensitivity and specificity on the logit scale using a bivariate
#   normal random-effects distribution. This is the recommended approach for
#   DTA meta-analysis per Cochrane Handbook Chapter 10.
#
# Dependencies:
#   - R >= 4.0.0
#   - mada (>= 0.5.11) for reitsma() bivariate GLMM
#   - jsonlite for JSON output
#
# Usage:
#   Rscript validate_dta_pro.R
#
# Output:
#   - Console: PASS/FAIL for each metric per dataset
#   - File: tests/r_validation_results.json
#
# References:
#   Reitsma JB, Glas AS, Rutjes AWS, Scholten RJPM, Bossuyt PM, Zwinderman AH.
#   Bivariate analysis of sensitivity and specificity produces informative
#   summary measures in diagnostic reviews. J Clin Epidemiol. 2005;58(10):982-990.
#
#   Doebler P (2020). mada: Meta-Analysis of Diagnostic Accuracy.
#   R package version 0.5.11. https://CRAN.R-project.org/package=mada
#
# Author: DTA Pro validation suite
# Date: 2026-03-13
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Package installation and loading
# ---------------------------------------------------------------------------

install_if_missing <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(sprintf("Installing package '%s' ...", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
  }
}

install_if_missing("mada")
install_if_missing("jsonlite")

suppressPackageStartupMessages(library(mada))
suppressPackageStartupMessages(library(jsonlite))

cat("\n")
cat("===========================================================================\n")
cat("  DTA Meta-Analysis Pro -- R Validation (mada::reitsma)\n")
cat("===========================================================================\n")
cat(sprintf("  R version:    %s\n", R.version.string))
cat(sprintf("  mada version: %s\n", packageVersion("mada")))
cat(sprintf("  Date:         %s\n", Sys.time()))
cat("===========================================================================\n\n")

# ---------------------------------------------------------------------------
# 2. Tolerance thresholds
# ---------------------------------------------------------------------------
#
# Point estimates (sensitivity, specificity): tight tolerance because the
# bivariate GLMM should converge to the same logit-scale estimates regardless
# of implementation.
#
# CI bounds: slightly wider because they depend on Wald vs profile likelihood
# and optimizer convergence.
#
# I-squared: generous tolerance because DTA heterogeneity can be computed via
# multiple methods (univariate Cochran's Q, Zhou generalized I2, Holling).
# We use the univariate fixed-effect approach for separate sens/spec I2.
#
# AUC: depends on SROC parameterization (Rutter-Gatsonis vs regression-based).
#
# DOR: exponential scale amplifies small logit differences.

TOL_POINT <- 0.01    # pooled sens/spec (probability scale)
TOL_CI    <- 0.02    # confidence interval bounds
TOL_I2    <- 5.0     # I-squared (percentage points)
TOL_AUC   <- 0.05    # area under SROC curve
TOL_DOR   <- 15.0    # diagnostic odds ratio (absolute; wide because DOR = exp(logDOR)
                      #   amplifies small logit-scale differences, and different methods
                      #   [summary point vs HSROC] yield different DOR estimates)
TOL_TAU2  <- 0.10    # between-study variance (logit scale)
TOL_RHO   <- 0.15    # correlation parameter

# ---------------------------------------------------------------------------
# 3. Helper functions
# ---------------------------------------------------------------------------

#' Compare a computed value against an expected value within tolerance
#'
#' @param label   Character label for the metric
#' @param computed Numeric value from R computation
#' @param expected Numeric value from the DTA Pro app
#' @param tol     Numeric tolerance (absolute difference)
#' @return List with label, computed, expected, diff, tol, pass, reason
check_metric <- function(label, computed, expected, tol) {
  if (is.null(computed) || is.na(computed) || is.null(expected) || is.na(expected)) {
    cat(sprintf("  [FAIL] %-22s  computed=%-10s  expected=%-10s  (NA value)\n",
                label, as.character(computed), as.character(expected)))
    return(list(
      label = label, computed = NA, expected = NA,
      diff = NA, tol = tol, pass = FALSE, reason = "NA value"
    ))
  }

  diff_val <- abs(computed - expected)
  pass <- diff_val <= tol
  status <- ifelse(pass, "PASS", "FAIL")

  cat(sprintf("  [%s] %-22s  computed=%.6f  expected=%.6f  diff=%.6f  tol=%.4f\n",
              status, label, computed, expected, diff_val, tol))

  list(
    label    = label,
    computed = round(computed, 6),
    expected = round(expected, 6),
    diff     = round(diff_val, 6),
    tol      = tol,
    pass     = pass,
    reason   = ifelse(pass, "within tolerance",
                       sprintf("diff %.6f > tol %.4f", diff_val, tol))
  )
}


#' Compute univariate I-squared for sensitivity or specificity
#'
#' Uses the logit-transformed proportion with inverse-variance weights from
#' a fixed-effect model (Cochran's Q). Continuity correction of 0.5 is added
#' to all cells (standard in DTA meta-analysis to handle zero cells).
#'
#' @param TP,FP,FN,TN Integer vectors of 2x2 cell counts
#' @param measure Character: "sens" or "spec"
#' @return Named list with Q, I2 (percentage), p_value
compute_I2_univariate <- function(TP, FP, FN, TN, measure = "sens") {
  k <- length(TP)

  if (measure == "sens") {
    num <- TP + 0.5
    den <- FN + 0.5
  } else {
    num <- TN + 0.5
    den <- FP + 0.5
  }

  logit_est <- log(num / den)
  var_est   <- 1 / num + 1 / den
  w         <- 1 / var_est

  theta_hat <- sum(w * logit_est) / sum(w)
  Q         <- sum(w * (logit_est - theta_hat)^2)
  df        <- k - 1
  I2        <- max(0, (Q - df) / Q * 100)
  p_value   <- 1 - pchisq(Q, df)

  list(Q = Q, I2 = I2, p_value = p_value)
}


# ---------------------------------------------------------------------------
# 4. Dataset definitions
# ---------------------------------------------------------------------------
# Each dataset is a list with:
#   - name: descriptive label
#   - data: data.frame with columns TP, FP, FN, TN (uppercase for mada)
#   - expected: named list of expected values from the DTA Pro app
#   - source: citation or provenance

datasets <- list()

# --- Dataset 1: Dementia (MMSE) -- 33 studies ---
# Source: mada::Dementia (R package); Mitchell AJ 2009
# These are the exact cell counts embedded in the DTA Pro app's validation
# dataset, matching the mada::Dementia data.
datasets$dementia <- list(
  name = "Dementia (MMSE) - 33 studies",
  source = "mada::Dementia; Mitchell 2009",
  data = data.frame(
    TP = c(66,118,48,134,24,68,64,282,14,262,144,184,22,112,152,30,32,10,
           708,182,60,74,28,40,318,388,118,44,124,26,74,38,78),
    FP = c(240,10,64,28,44,48,0,20,44,30,30,34,152,590,126,26,4,12,
           1438,18,24,16,26,76,174,16,2,34,98,4,2,0,46),
    FN = c(4,12,20,8,6,16,18,64,2,20,18,34,0,0,82,26,6,4,
           88,108,30,24,12,6,52,116,66,8,46,44,32,46,34),
    TN = c(870,110,990,152,292,154,72,286,286,178,124,52,140,2092,1010,236,248,334,
           10448,184,74,144,210,528,578,54,44,396,310,172,226,440,376)
  ),
  expected = list(
    pooled_sens    = 0.788737,
    pooled_sens_lo = 0.735469,
    pooled_sens_hi = 0.833704,
    pooled_spec    = 0.886182,
    pooled_spec_lo = 0.847650,
    pooled_spec_hi = 0.915934,
    auc            = 0.904416,
    dor            = 32.59,
    sens_I2        = 93.1,
    spec_I2        = 98.3
  )
)

# --- Dataset 2: Scheidler MRI -- 8 studies ---
# Source: Scheidler et al. Radiology 1997
# MRI for lymph node metastases in cervical cancer.
# Contains zero cells (FP=0, FN=0 in some studies).
datasets$scheidler <- list(
  name = "Scheidler MRI - 8 studies",
  source = "Scheidler et al. Radiology 1997",
  data = data.frame(
    TP = c(9, 4, 3, 11, 8, 5, 22, 7),
    FP = c(3, 0, 0, 0, 2, 1, 2, 1),
    FN = c(1, 2, 1, 0, 2, 0, 3, 4),
    TN = c(25, 18, 10, 15, 35, 10, 43, 14)
  ),
  expected = list(
    pooled_sens    = 0.783185,
    pooled_sens_lo = 0.675581,
    pooled_sens_hi = 0.862369,
    pooled_spec    = 0.923134,
    pooled_spec_lo = 0.871856,
    pooled_spec_hi = 0.954952,
    auc            = 0.929533,
    dor            = 57.92,
    sens_I2        = 1.4,
    spec_I2        = 0.0
  )
)

# --- Dataset 3: CD64 Sepsis -- 10 studies ---
# Source: Meta-analysis of neutrophil CD64 expression for sepsis diagnosis
datasets$cd64 <- list(
  name = "CD64 Sepsis - 10 studies",
  source = "CD64 sepsis meta-analysis 2010",
  data = data.frame(
    TP = c(28, 19, 18, 25, 29, 22, 45, 33, 18, 41),
    FP = c(5, 8, 4, 7, 3, 6, 8, 5, 7, 4),
    FN = c(4, 6, 3, 8, 5, 4, 7, 6, 4, 6),
    TN = c(38, 42, 35, 40, 48, 33, 55, 44, 36, 52)
  ),
  expected = list(
    pooled_sens    = 0.836365,
    pooled_sens_lo = 0.791813,
    pooled_sens_hi = 0.872913,
    pooled_spec    = 0.876155,
    pooled_spec_lo = 0.842680,
    pooled_spec_hi = 0.903324,
    auc            = 0.917607,
    dor            = 32.38,
    sens_I2        = 0.0,
    spec_I2        = 0.0
  )
)

# --- Dataset 4: Glas 2003 (FDG-PET) -- 9 studies ---
# Source: Glas AS, et al. J Clin Epidemiol 2003;56:1129-1135
# FDG-PET for recurrent colorectal cancer detection.
# Contains zero cells (FP=0, FN=0 in some studies).
# This 9-study dataset serves as a benchmark-only entry: R computes the
# reference values, which can then be verified against the app interactively.
datasets$glas <- list(
  name = "Glas 2003 (FDG-PET) - 9 studies",
  source = "Glas et al. J Clin Epidemiol 2003",
  data = data.frame(
    TP = c(19, 12, 18, 3, 26, 25, 31, 8, 12),
    FP = c(0, 3, 1, 6, 7, 3, 0, 1, 1),
    FN = c(1, 3, 0, 0, 6, 3, 0, 1, 4),
    TN = c(20, 11, 8, 35, 39, 14, 30, 51, 19)
  ),
  expected = list(
    pooled_sens    = NULL,
    pooled_spec    = NULL,
    auc            = NULL,
    dor            = NULL,
    sens_I2        = NULL,
    spec_I2        = NULL
  )
)

# ---------------------------------------------------------------------------
# 5. Run validation for each dataset
# ---------------------------------------------------------------------------

all_results <- list()
total_pass  <- 0
total_fail  <- 0
total_skip  <- 0

for (ds_key in names(datasets)) {
  ds <- datasets[[ds_key]]

  cat(sprintf("\n--- %s ---\n", ds$name))
  cat(sprintf("  Source: %s\n", ds$source))
  cat(sprintf("  k = %d studies\n\n", nrow(ds$data)))

  # ---- 5a. Fit the bivariate GLMM (Reitsma model) ----
  # suppressWarnings: zero cells trigger continuity correction warnings from
  # the internal logit transformation. mada adds 0.5 to all cells by default.
  fit <- suppressWarnings(
    reitsma(ds$data, method = "reml")
  )

  s <- summary(fit)

  # ---- 5b. Extract pooled estimates ----
  # mada::reitsma parameterizes as (tsens, tfpr) on the logit scale.
  # The summary coefficients table includes both logit-scale and
  # back-transformed (probability-scale) rows:
  #   Row "sensitivity"     = back-transformed sensitivity (with CI)
  #   Row "false pos. rate" = back-transformed FPR (with CI)

  pooled_sens <- s$coefficients["sensitivity", "Estimate"]
  sens_lo     <- s$coefficients["sensitivity", "95%ci.lb"]
  sens_hi     <- s$coefficients["sensitivity", "95%ci.ub"]

  pooled_fpr  <- s$coefficients["false pos. rate", "Estimate"]
  fpr_lo      <- s$coefficients["false pos. rate", "95%ci.lb"]
  fpr_hi      <- s$coefficients["false pos. rate", "95%ci.ub"]

  # Specificity = 1 - FPR (CI bounds invert)
  pooled_spec <- 1 - pooled_fpr
  spec_lo     <- 1 - fpr_hi
  spec_hi     <- 1 - fpr_lo

  # ---- 5c. Variance components (between-study) ----
  Sigma     <- fit$Psi  # 2x2 variance-covariance matrix on logit scale
  tau2_sens <- Sigma[1, 1]
  tau2_fpr  <- Sigma[2, 2]

  # Correlation: guard against NaN when tau2 is ~0
  if (tau2_sens > 1e-10 && tau2_fpr > 1e-10) {
    rho <- Sigma[1, 2] / sqrt(tau2_sens * tau2_fpr)
  } else {
    rho <- NA_real_
  }

  # ---- 5d. I-squared (univariate, logit-scale, Cochran's Q) ----
  # Separate I2 for sensitivity and specificity, computed via fixed-effect
  # inverse-variance weighting on the logit scale. This is the standard
  # approach in DTA software (e.g., Review Manager, MetaDTA).
  i2_sens <- compute_I2_univariate(ds$data$TP, ds$data$FP, ds$data$FN, ds$data$TN, "sens")
  i2_spec <- compute_I2_univariate(ds$data$TP, ds$data$FP, ds$data$FN, ds$data$TN, "spec")

  sens_I2 <- i2_sens$I2
  spec_I2 <- i2_spec$I2

  # ---- 5e. DOR from pooled estimates ----
  # DOR = (sens / (1-sens)) * (spec / (1-spec))
  #     = exp(logit_sens) * exp(logit_spec)   where logit_spec = -logit_fpr
  #     = exp(logit_sens - logit_fpr)
  # This is the standard summary DOR from the bivariate model summary point.
  dor <- (pooled_sens / (1 - pooled_sens)) * (pooled_spec / (1 - pooled_spec))

  # ---- 5f. AUC from mada::AUC ----
  # Uses the Rutter-Gatsonis SROC parameterization (default in mada)
  auc_obj <- tryCatch(AUC(fit), error = function(e) NULL)
  auc     <- if (!is.null(auc_obj) && !is.null(auc_obj$AUC)) auc_obj$AUC else NA_real_

  # ---- 5g. Print all computed values ----
  cat("  Computed R values:\n")
  cat(sprintf("    Pooled Sensitivity: %.6f [%.6f, %.6f]\n", pooled_sens, sens_lo, sens_hi))
  cat(sprintf("    Pooled Specificity: %.6f [%.6f, %.6f]\n", pooled_spec, spec_lo, spec_hi))
  cat(sprintf("    tau2(sens):  %.6f\n", tau2_sens))
  cat(sprintf("    tau2(fpr):   %.6f\n", tau2_fpr))
  cat(sprintf("    rho:         %s\n", ifelse(is.na(rho), "NA (tau2 ~ 0)", sprintf("%.4f", rho))))
  cat(sprintf("    I2(sens):    %.1f%%  (Q=%.2f, p=%.4f)\n", sens_I2, i2_sens$Q, i2_sens$p_value))
  cat(sprintf("    I2(spec):    %.1f%%  (Q=%.2f, p=%.4f)\n", spec_I2, i2_spec$Q, i2_spec$p_value))
  cat(sprintf("    DOR:         %.2f\n", dor))
  cat(sprintf("    AUC:         %.6f\n", auc))
  cat("\n")

  # ---- 5h. Compare against expected values ----
  checks <- list()
  exp <- ds$expected

  if (!is.null(exp$pooled_sens)) {
    cat("  Validation checks:\n")
    checks$pooled_sens <- check_metric("Pooled Sensitivity",  pooled_sens, exp$pooled_sens,    TOL_POINT)
    checks$pooled_spec <- check_metric("Pooled Specificity",  pooled_spec, exp$pooled_spec,    TOL_POINT)
    checks$sens_ci_lo  <- check_metric("Sens CI Lower",       sens_lo,     exp$pooled_sens_lo, TOL_CI)
    checks$sens_ci_hi  <- check_metric("Sens CI Upper",       sens_hi,     exp$pooled_sens_hi, TOL_CI)
    checks$spec_ci_lo  <- check_metric("Spec CI Lower",       spec_lo,     exp$pooled_spec_lo, TOL_CI)
    checks$spec_ci_hi  <- check_metric("Spec CI Upper",       spec_hi,     exp$pooled_spec_hi, TOL_CI)
    checks$auc         <- check_metric("AUC",                 auc,         exp$auc,            TOL_AUC)
    checks$dor         <- check_metric("DOR",                 dor,         exp$dor,            TOL_DOR)
    checks$sens_I2     <- check_metric("I2 (Sensitivity)",    sens_I2,     exp$sens_I2,        TOL_I2)
    checks$spec_I2     <- check_metric("I2 (Specificity)",    spec_I2,     exp$spec_I2,        TOL_I2)

    for (ch in checks) {
      if (ch$pass) total_pass <- total_pass + 1 else total_fail <- total_fail + 1
    }
  } else {
    cat("  [INFO] No expected values provided -- reporting R-computed values only.\n")
    cat("         These values can be used as reference for app validation.\n")
    total_skip <- total_skip + 1
  }

  # ---- 5i. Store results ----
  all_results[[ds_key]] <- list(
    dataset  = ds$name,
    source   = ds$source,
    k        = nrow(ds$data),
    method   = "reitsma() bivariate GLMM, REML",
    computed = list(
      pooled_sens    = round(pooled_sens, 6),
      pooled_sens_ci = c(round(sens_lo, 6), round(sens_hi, 6)),
      pooled_spec    = round(pooled_spec, 6),
      pooled_spec_ci = c(round(spec_lo, 6), round(spec_hi, 6)),
      tau2_sens      = round(tau2_sens, 6),
      tau2_fpr       = round(tau2_fpr, 6),
      rho            = if (is.na(rho)) "NA" else round(rho, 4),
      sens_I2        = round(sens_I2, 1),
      spec_I2        = round(spec_I2, 1),
      dor            = round(dor, 2),
      auc            = round(auc, 6)
    ),
    checks = checks
  )
}

# ---------------------------------------------------------------------------
# 6. Overall summary
# ---------------------------------------------------------------------------

cat("\n")
cat("===========================================================================\n")
cat("  VALIDATION SUMMARY\n")
cat("===========================================================================\n")
cat(sprintf("  Datasets validated:  %d (%d with expected values, %d benchmark-only)\n",
            length(datasets), length(datasets) - total_skip, total_skip))
cat(sprintf("  Total checks:        %d\n", total_pass + total_fail))
cat(sprintf("  PASSED:              %d\n", total_pass))
cat(sprintf("  FAILED:              %d\n", total_fail))
cat("===========================================================================\n")

if (total_fail == 0 && total_pass > 0) {
  cat("  RESULT: ALL CHECKS PASSED\n")
} else if (total_fail > 0) {
  cat(sprintf("  RESULT: %d CHECK(S) FAILED -- review details above\n", total_fail))
} else {
  cat("  RESULT: NO CHECKS PERFORMED (no expected values provided)\n")
}

cat("===========================================================================\n\n")

# ---------------------------------------------------------------------------
# 7. Save results to JSON
# ---------------------------------------------------------------------------

output <- list(
  metadata = list(
    description  = "DTA Meta-Analysis Pro -- R validation results",
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S"),
    r_version    = R.version.string,
    mada_version = as.character(packageVersion("mada")),
    method       = "reitsma() bivariate GLMM, REML",
    reference    = paste0(
      "Reitsma JB, Glas AS, Rutjes AWS, Scholten RJPM, Bossuyt PM, ",
      "Zwinderman AH. Bivariate analysis of sensitivity and specificity ",
      "produces informative summary measures in diagnostic reviews. ",
      "J Clin Epidemiol. 2005;58(10):982-990."
    ),
    tolerances = list(
      point_estimate = TOL_POINT,
      ci_bound       = TOL_CI,
      I2             = TOL_I2,
      auc            = TOL_AUC,
      dor            = TOL_DOR,
      tau2           = TOL_TAU2,
      rho            = TOL_RHO
    )
  ),
  summary = list(
    total_datasets    = length(datasets),
    datasets_checked  = length(datasets) - total_skip,
    datasets_benchmark = total_skip,
    total_checks      = total_pass + total_fail,
    passed            = total_pass,
    failed            = total_fail,
    overall           = ifelse(total_fail == 0 && total_pass > 0, "PASS", "FAIL")
  ),
  datasets = all_results
)

# Determine output path: same directory as this script
script_dir <- tryCatch({
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    dirname(normalizePath(sub("--file=", "", file_arg[1])))
  } else {
    getwd()
  }
}, error = function(e) getwd())

json_path <- file.path(script_dir, "r_validation_results.json")

writeLines(
  toJSON(output, pretty = TRUE, auto_unbox = TRUE, digits = 6),
  json_path
)

cat(sprintf("Results saved to: %s\n\n", json_path))

# Return non-zero exit code if any checks failed (useful for CI/CD)
if (total_fail > 0) {
  quit(status = 1, save = "no")
}
