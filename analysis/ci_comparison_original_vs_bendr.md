# Confidence intervals and NtC: original script vs. BendR

Comparison of the confidence interval (CI) and NtC calculation in the reference script of
Stadnicka-Michalak et al. (`inst/reference/StadnickaNtC.R`) and in BendR 0.3.0
(https://github.com/alxbetz/bendr, commit 0f2d407). Status of each issue in the fork
https://github.com/schuerc/bendr and in this package is given at the end.

Model (both): `survival = 100 / (1 + 10^((logEC50 - log10(conc)) * slope))`, where survival is
100% - effect, so the curve decreases and the slope is negative.

## 1. Shared approach

- Delta method: `CI(x) = t * sqrt(J(x) Cov(theta) J(x)')`, a pointwise band on the mean curve.
- Two parameters; top and bottom fixed at 100% and 0%.
- One pooled residual variance (equal scatter assumed at all concentrations).
- The band is not bounded to 0-100%.
- NtC rule: upper-CI criterion, lower-CI (> 90% survival) criterion, narrow-CI (10x ratio)
  correction, override by the measured values.

## 2. Differences

| Aspect | StadnickaNtC.R | BendR 0.3.0 |
|---|---|---|
| Control row | Always drops the first data row, whatever its concentration | Rows with concentration 0 removed only in `cleanAndPivot` |
| Fitted data / df | Concentration means; df = `nrow(file) - 3` (= number of means - 2 when one row is dropped) | All replicate values; df = number of values - 2 |
| Optimiser | `nlxb`, then `nls2(brute-force)` evaluated at that point | `nlsLM`, fallback `nls2(grid-search)` |
| Range for the CI criteria | Fixed: from log 0 (if the lowest used log concentration is > 0) or 3 x the lowest log concentration, to 2 x the highest log concentration | +/-20% of the tested range; `extend_ci` down to 99.9% survival only if the lower band never exceeds 90% |
| Upper criterion | upper band >= 99.999999999 | upper band > 100 |
| Measured values | Index arithmetic; fails with a single replicate column | `dplyr::filter`, robust |
| Extras | - | Bootstrap option, EC50/EC10 CIs |

## 3. Issues found

1. **Range too short (both).** A criterion whose crossing lies outside the range is either "not
   found" (script: `max(which())` on an empty set, which fails) or, for the upper criterion in
   BendR, treated as effect 0, which triggers the narrow-CI correction and gives a too low NtC.
   Quantified in `ntc_decomposition.md`.
2. **Unguarded `max(which(...))`** for the lower criterion (both); 0/0 in the narrow-CI ratio when
   both effects are 0.
3. **Bootstrap row overwrite (BendR):** `parameters[1,] <-` inside the loop overwrote the
   original-fit row on every iteration. Also, parameters were stored by position while `coef()`
   order follows the start vector, so slope and logEC50 were swapped between rows.
4. **Bootstrap fits not optimised (BendR):** `nls2(grid-search)` returns the best grid point, so the
   estimates are quantised; the grid (0.5-2 x start) collapses when logEC50 is near 0; no
   `tryCatch`, so one failed fit aborted the run; case resampling could drop concentration levels.
5. **EC10 CI (BendR):** the reparametrised EC10 formula was defined but not used; the refit used the
   original formula with `brute-force` from one start, so the "EC10 CI" was meaningless.
6. **Bootstrap and `extend_ci` (BendR):** `extend_ci` always used the delta method and the default
   confidence level, mixing methods; bootstrap bands never exceed 100%, so the upper criterion can
   never be met with them.
7. **Flat curves (BendR):** the default extension to 99.9% survival lies extremely far below the
   data for very flat curves; the grid could not be allocated.

## 4. Conceptual issues (both)

- The band is unbounded; "upper band reaches 100%" depends on that.
- Unequal scatter: replicate scatter is typically largest in the steep part of the curve and small
  near 0% and 100% survival; the pooled variance is wrong at both ends.
- The degrees of freedom depend on the replication unit. Wells on one plate are pseudo-replicates:
  treating them as independent makes the band too narrow.
- Pointwise, not simultaneous, bands; top and bottom fixed.

## 5. Status

| Issue | bendr fork (>= 0.4.0) | NtC package |
|---|---|---|
| 1 Range too short | Fixed: range extended until both criteria are found | Updated algorithm uses it; original reproduces the script and reports the failure |
| 2 Guards | Fixed: NA with a message | Message in the app |
| 3 Bootstrap overwrite / column order | Fixed | - |
| 4 Bootstrap fits | Fixed: `nlsLM` from the original estimates, failures counted, resampling within concentrations | Bootstrap option for EC50/EC10 |
| 5 EC10 CI | Fixed: fit of the reparametrised model | EC10 with CI shown |
| 6 Bootstrap and NtC | Fixed: NtC always from the delta-method band | Same |
| 7 Flat curves | Fixed (0.4.1): default extension limited | Message in the app |
| Replication unit | `replicates = "same_plate"` / `"independent"` | Checkbox, default same plate |
| Unbounded band, unequal scatter | Open | Unequal scatter flagged in the fit diagnostics |

## 6. Possible further improvements (not implemented)

1. Bounded band: compute the band on the linear predictor `(logEC50 - log conc) * slope` and
   transform. This changes the upper criterion (a bounded band never exceeds 100%), so the NtC rule
   would need a new threshold.
2. Model the unequal scatter: weighted least squares or `nlme::gnls` with a variance function.
3. Coverage simulation: simulate from fitted curves and compare the actual coverage of the bands and
   the bias of the NtC for the different settings.
