# Updated algorithm, based on the bendr fork (https://github.com/schuerc/bendr) -----------

#' Updated NtC algorithm
#'
#' Fit with Levenberg-Marquardt (bendr), delta-method confidence band, concentration grid
#' extended until both CI criteria are found, guards for all failure cases.
#'
#' @inheritParams ntc_original
#' @param replicates `"same_plate"` (default): replicate columns are wells on one plate; the
#'   curve is fitted to the concentration means, df = n_concentrations - 2.
#'   `"independent"`: replicate columns are independent experiments; the curve is fitted to all
#'   values, df = n_values - 2.
#' @param level confidence level of the band (default 0.95)
#' @param ec_ci method for the EC50/EC10/slope confidence intervals: `"delta"` (default) or
#'   `"bootstrap"` (1000 resamples within concentrations; independent replicates only). The NtC
#'   always uses the delta-method band.
#'
#' @return list with the NtC, intermediate values, the fit and the confidence band
#' @export
ntc_updated <- function(data, measured_criterion = TRUE, narrow_ci_correction = TRUE,
                        replicates = c("same_plate", "independent"), level = 0.95,
                        ec_ci = c("delta", "bootstrap")) {
  replicates <- match.arg(replicates)
  ec_ci <- match.arg(ec_ci)
  fail <- function(msg) ntc_abort(paste("Updated algorithm:", msg), "ntc_algorithm_error")
  if (!is.numeric(level) || level <= 0.5 || level >= 1) fail("the confidence level must be between 0.5 and 1.")
  if (ec_ci == "bootstrap" && replicates != "independent") {
    fail("bootstrap confidence intervals need independent replicates.")
  }

  eff <- data$effects
  df <- data.frame(conc = data$conc, eff, check.names = FALSE)
  names(df) <- c("conc", paste0("rep", seq_len(ncol(eff))))
  warns <- character()
  collect <- function(expr) {
    withCallingHandlers(expr,
      warning = function(w) { warns <<- c(warns, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) { warns <<- c(warns, trimws(conditionMessage(m))); invokeRestart("muffleMessage") })
  }

  fo <- tryCatch(collect(bendr::fitdr_replicates(drc_formula(), df, 2:ncol(df), conc, level = level,
                                                 ci_method = ec_ci, replicates = replicates)),
                 error = function(e) fail(paste("the curve fit failed:", conditionMessage(e))))
  if (all(is.na(fo$plot.data$ci.values.lower))) {
    fail("the confidence band could not be computed (the covariance matrix of the fit is not available).")
  }
  cf <- coef(fo$curve.fit)
  if (cf[["slope"]] >= 0) {
    fail("the fitted curve does not decrease: survival increases with concentration. Check that the values are survival in % (100 = no effect), not effect in %.")
  }
  d <- tryCatch(collect(bendr::findNTC(fo, measured.criterion = measured_criterion,
                                        narrow.ci.correction = narrow_ci_correction, details = TRUE)),
                error = function(e) fail(conditionMessage(e)))
  if (is.na(d$NtC)) fail(paste0(d$message, "."))

  resid <- unname(residuals(fo$curve.fit))
  conv <- tryCatch(isTRUE(fo$curve.fit$convInfo$isConv), error = function(e) NA)
  band <- data.frame(log_conc = d$plot.data$log.concentration, fit = d$plot.data$curve.predict,
                     lower = d$plot.data$ci.values.lower, upper = d$plot.data$ci.values.upper)
  band <- band[order(band$log_conc), ]

  list(NtC = d$NtC, NtC_model = d$NtC_model, Effect_NtC = d$Effect_NtC, rule = d$rule,
       NtC_lower = d$NtC_lower, NtC_upper = d$NtC_upper,
       Effect_NtC_lower = d$Effect_NtC_lower, Effect_NtC_upper = d$Effect_NtC_upper,
       measured_applied = d$measured.applied, NtC_measured = d$NtC_measured,
       logEC50 = unname(cf["logEC50"]), slope = unname(cf["slope"]), EC50 = unname(fo$ec50),
       EC50_ci = unname(fo$ec50.ci), EC10 = unname(fo$ec10), EC10_ci = unname(fo$ec10.ci),
       slope_ci = unname(fo$slope.ci), ec_ci = ec_ci,
       se = sqrt(diag(vcov(fo$curve.fit)))[c("logEC50", "slope")],
       df = nrow(fo$data) - 2, t = qt(1 - (1 - level) / 2, nrow(fo$data) - 2), level = level,
       fit_data = data.frame(log_conc = fo$data$logconc, effect = fo$data$effect),
       fit_to = if (replicates == "same_plate") "concentration means" else "all replicate values",
       n_fit = nrow(fo$data), residuals = resid, converged = conv,
       band = band, range = range(band$log_conc), n_extensions = d$n.extensions,
       warnings = unique(warns[nzchar(warns)]))
}
