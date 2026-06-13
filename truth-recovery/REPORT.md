# Truth-recovery yardstick — dta-meta-analysis-pro

**Verdict: REAL but MODEST improvement (measured) + an honest labelling correction.**

The app is `mada::reitsma`-validated on its point outputs. That proves the
formulas; it does not show whether the pooled Se/Sp **confidence intervals
recover the true sensitivity/specificity** under bivariate between-study
heterogeneity. This harness supplies that evidence.

## Method
- DGP (`dgp-dta.mjs`): bivariate logit-normal random effects
  `(uA,uB) ~ N(0, Σ)`, `sens_i = expit(logit(Se)+uA)`, `spec_i = expit(logit(Sp)+uB)`,
  `TP ~ Bin(nDis, sens_i)`, `TN ~ Bin(nHeal, spec_i)`. Seeded → reproducible.
- Estimand: a logit RE pool targets `expit(muA)=Se`, `expit(muB)=Sp`.
- Engine (`engine.mjs`): `improvedBivariatePool` + helpers copied **verbatim**
  from `dta-pro-v3.7.html`. The SAME math the app ships is measured.
- 400–1000 reps/cell, `Se=0.85, Sp=0.80`, `k∈{5,10,20}`, four heterogeneity
  presets (`het_low/mod/high` + a negative-Se/Sp-correlation threshold preset).

## Headline (1000 reps/cell, mean over all 12 cells)

| method        | mean cov(Se) | mean cov(Sp) |
|---------------|-------------:|-------------:|
| shipped       |       0.9445 |       0.9452 |
| genuine HKSJ  |       0.9518 |       0.9546 |

Per-cell, the gap is largest under high heterogeneity:

| scenario  | k | ship cov(Sp) | HKSJ cov(Sp) |
|-----------|---|-------------:|-------------:|
| het_high  | 10|      0.93    |    0.95      |
| het_high  | 20|      0.90    |    0.93      |
| het_mod   | 20|      0.9375  |    0.945     |

## Findings (all measured)
1. **Labelling correction (the substantive finding):** the shipped interval
   reports `ciMethod: 'HKSJ-t'`, but the code only swaps the z critical value for
   `t_{k-2}` — it uses the *ordinary* RE standard error `1/√Σw*`. It does **not**
   apply the Hartung–Knapp–Sidik–Jonkman variance adjustment (the `Q/(k-1)`
   inflation, floored at 1). So the label overstates what the method does.
2. **Consequence:** with only the t-swap, the interval is slightly too narrow
   under heterogeneity — mean coverage ≈ 0.945 (below nominal), dropping to 0.90
   under high heterogeneity with k=20.
3. **Fix:** adding the genuine HKSJ adjustment (same point estimate, inflated SE)
   lifts mean coverage to ≈ 0.95 (on nominal) and **never lowers it** in any cell;
   gains are +2–3pp exactly where the shipped interval is worst (high het).
   → **Recommendation:** either rename the method honestly to `t_{k-2}` *or*
   implement the genuine HKSJ inflation (provided in `harness.mjs::reAxis`). The
   latter is the better fix and is statistically preferred for k<30.
4. The improvement is **modest** — the t-crit already does most of the small-k
   work, which is why the shipped interval is only ~0.5pp under nominal on
   average, not catastrophically so. Honest framing: a refinement, not a rescue.

## What did NOT transfer / scope notes
- NPE/conformal/SBC/PartialID are estimator-of-μ machinery; this is a moment-based
  bivariate pool, so only the known-truth harness + bivariate DGP transferred.
- The harness measures the moment-based pool the app ships; it does not test the
  full joint REML/HSROC/copula paths (those use heavier numerics). The Se/Sp
  marginal-coverage finding applies to the default `improvedBivariatePool` output,
  which feeds the app's headline pooled Se/Sp + CIs.

## Reproduce
```
node truth-recovery/harness.mjs --reps 400
node --test truth-recovery/test-truth-recovery.mjs
```
