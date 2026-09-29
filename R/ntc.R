#' Determine the non-toxic concentration (NtC)
#'
#' @param data an `ntc_data` object ([ntc_prepare()]), a data frame (first column
#'   concentrations, other columns replicates) or a path to a CSV file
#' @param algorithm `"updated"` (default) or `"original"` (Stadnicka-Michalak et al.)
#' @param control is the first row a control? `"auto"`: yes if its concentration is 0.
#'   Only used when `data` is not already an `ntc_data` object.
#' @param measured_criterion apply the measured value criterion (default TRUE)
#' @param narrow_ci_correction apply the correction for narrow confidence intervals (default TRUE)
#' @param replicates,level,ec_ci settings of the updated algorithm, see [ntc_updated()];
#'   ignored for the original algorithm
#'
#' @return an object of class `ntc_result`. If the calculation fails, `status` is `"failed"`
#'   and `message` explains why.
#' @export
#' @examples
#' f <- system.file("extdata", "demo_semicolon.csv", package = "NtC")
#' r <- ntc(f)
#' r
ntc <- function(data, algorithm = c("updated", "original"), control = "auto",
                measured_criterion = TRUE, narrow_ci_correction = TRUE,
                replicates = c("same_plate", "independent"), level = 0.95,
                ec_ci = c("delta", "bootstrap")) {
  algorithm <- match.arg(algorithm)
  replicates <- match.arg(replicates)
  ec_ci <- match.arg(ec_ci)
  if (is.character(data)) data <- read_ntc_file(data)
  if (!inherits(data, "ntc_data")) data <- ntc_prepare(data, control = control)

  settings <- list(algorithm = algorithm, control = data$control,
                   measured_criterion = measured_criterion, narrow_ci_correction = narrow_ci_correction)
  if (algorithm == "updated") {
    settings <- c(settings, list(replicates = replicates, level = level, ec_ci = ec_ci))
  } else {
    settings <- c(settings, list(replicates = "fit to concentration means", level = 0.95, ec_ci = "none"))
  }

  res <- tryCatch({
    fit <- if (algorithm == "original") {
      ntc_original(data, measured_criterion, narrow_ci_correction)
    } else {
      ntc_updated(data, measured_criterion, narrow_ci_correction, replicates, level, ec_ci)
    }
    list(status = "ok", message = NA_character_, fit = fit)
  }, ntc_error = function(e) list(status = "failed", message = conditionMessage(e), fit = NULL),
  error = function(e) list(status = "failed", message = paste("Unexpected error:", conditionMessage(e)), fit = NULL))

  out <- structure(c(res, list(data = data, settings = settings, versions = ntc_versions(),
                               time = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"))),
                   class = "ntc_result")
  out$diagnostics <- ntc_diagnostics(out)
  out
}

#' Fit diagnostics and warnings
#'
#' @param result an `ntc_result`
#' @return data frame with columns `level` ("info" or "warning") and `message`
#' @export
ntc_diagnostics <- function(result) {
  d <- data.frame(level = character(), message = character())
  add <- function(level, ...) d <<- rbind(d, data.frame(level = level, message = sprintf(...)))
  dat <- result$data
  f <- result$fit
  means <- rowMeans(dat$effects, na.rm = TRUE)

  add("info", "%d concentrations, %d replicate column(s), %d values used.", length(dat$conc),
      ncol(dat$effects), sum(!is.na(dat$effects)))
  for (n in dat$notes) add("info", "%s", n)
  if (ncol(dat$effects) == 1) add("warning", "Only one replicate column: the measured value criterion rests on single values.")
  if (means[1] < 90) {
    add("warning", "Mean survival at the lowest tested concentration is %.1f%% (< 90%%): no concentration without effect was tested, so the NtC is extrapolated below the tested range.", means[1])
  }
  if (min(means) > 50) add("warning", "Survival never falls below 50%%: the EC50 lies above the tested range (partial curve).")
  if (max(means) < 50) add("warning", "Survival is below 50%% at all concentrations: the EC50 lies below the tested range (partial curve).")
  if (ncol(dat$effects) >= 2) {
    sds <- apply(dat$effects, 1, sd, na.rm = TRUE)
    ok <- is.finite(sds)
    if (sum(ok) >= 2 && max(sds[ok]) > 5 * max(min(sds[ok]), 0.5)) {
      add("warning", "The scatter of the replicates differs strongly between concentrations (SD %.1f to %.1f). The confidence band assumes equal scatter.",
          min(sds[ok]), max(sds[ok]))
    }
  }
  if (is.null(f)) return(d)

  add("info", "Fit to %s: %d points, %d parameters, df = %d (t = %.2f), %s%% confidence band.",
      f$fit_to, f$n_fit, 2L, f$df, f$t, format(100 * f$level))
  if (isFALSE(f$converged)) add("warning", "The fit did not report convergence.")
  tss <- sum((f$fit_data$effect - mean(f$fit_data$effect))^2)
  r2 <- 1 - sum(f$residuals^2) / tss
  add(if (is.finite(r2) && r2 < 0.9) "warning" else "info", "R\u00b2 = %.3f, residual SD = %.2f percentage points.",
      r2, sqrt(sum(f$residuals^2) / f$df))
  if (f$slope >= 0) add("warning", "The fitted slope is not negative: survival does not decrease with concentration.")
  if (f$n_extensions > 0) {
    add("info", "The concentration range was extended by %d log unit(s) below the tested range to find both confidence criteria.", f$n_extensions)
  }
  if (is.finite(f$NtC) && f$NtC < min(dat$conc)) {
    add("warning", "The NtC (%s) is below the lowest tested concentration (%s): it is an extrapolation.", fmt_num(f$NtC), fmt_num(min(dat$conc)))
  }
  if (is.finite(f$NtC) && f$NtC > max(dat$conc)) {
    add("warning", "The NtC (%s) is above the highest tested concentration (%s): no toxic effect was reached, the value is an extrapolation.", fmt_num(f$NtC), fmt_num(max(dat$conc)))
  }
  if (f$rule == "narrow_ci_correction") add("info", "The NtC was set by the narrow confidence interval correction.")
  if (isTRUE(f$measured_applied)) add("info", "The NtC was set by the measured value criterion.")
  for (w in f$warnings) add("warning", "%s", w)
  d
}

#' Results as a table
#'
#' @param result an `ntc_result`
#' @return data frame with one row per quantity
#' @export
ntc_table <- function(result) {
  s <- result$settings
  f <- result$fit
  row <- function(q, v, note = "") data.frame(quantity = q, value = v, note = note)
  rule_txt <- c(lower = "lower-CI criterion", upper = "upper-CI criterion",
                narrow_ci_correction = "narrow-CI correction")
  base <- rbind(
    row("status", result$status, if (is.na(result$message)) "" else result$message),
    row("algorithm", s$algorithm, if (s$algorithm == "original") "Stadnicka-Michalak et al., as published" else "based on bendr (schuerc fork)"))
  if (is.null(f)) return(base)
  ci <- function(x) if (length(x) == 2 && all(is.finite(x))) sprintf("%s to %s", fmt_num(x[1]), fmt_num(x[2])) else ""
  rbind(base,
    row("NtC", fmt_num(f$NtC), if (isTRUE(f$measured_applied)) "set by measured value criterion" else paste("from", rule_txt[[f$rule]])),
    row("NtC_lowerCI", fmt_num(f$NtC_lower), sprintf("effect %s%%", fmt_num(f$Effect_NtC_lower, 3))),
    row("NtC_upperCI", fmt_num(f$NtC_upper), sprintf("effect %s%%", fmt_num(f$Effect_NtC_upper, 3))),
    row("NtC_measured", fmt_num(f$NtC_measured), "lowest tested concentration with >= 10% effect"),
    row("Effect at NtC (model)", fmt_num(f$Effect_NtC, 3), "%"),
    row("EC50", fmt_num(f$EC50), ci(f$EC50_ci)),
    row("EC10", if (is.null(f$EC10)) "" else fmt_num(f$EC10), ci(f$EC10_ci)),
    row("slope", fmt_num(f$slope), ci(f$slope_ci)),
    row("degrees of freedom", as.character(f$df), f$fit_to),
    row("confidence level", format(f$level), if (!is.null(f$ec_ci) && f$ec_ci == "bootstrap") "EC CIs: bootstrap" else ""))
}

#' @export
print.ntc_result <- function(x, ...) {
  cat(sprintf("NtC (%s algorithm), NtC %s, bendr %s\n", x$settings$algorithm, x$versions[["NtC"]], x$versions[["bendr"]]))
  if (x$status != "ok") cat("FAILED:", x$message, "\n")
  print(ntc_table(x)[-(1:2), ], row.names = FALSE)
  w <- x$diagnostics[x$diagnostics$level == "warning", "message"]
  if (length(w)) cat(paste0("Warning: ", w, collapse = "\n"), "\n")
  invisible(x)
}
