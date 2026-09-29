# Original algorithm (Stadnicka-Michalak et al.) ------------------------------
#
# Port of inst/reference/StadnickaNtC.R. The calculation follows the script line by line,
# including its quirks. Differences:
# - control row: handled by ntc_prepare() (the script always drops the first row);
# - degrees of freedom: number of fitted concentrations - 2 (the script's nrow(file) - 3 is
#   the same when one row is dropped as control);
# - where the script would fail or return a meaningless value, an ntc_error explains why;
# - the two criteria can be switched off (as in the published app).

#' Original NtC algorithm (Stadnicka-Michalak et al.)
#'
#' @param data an `ntc_data` object from [ntc_prepare()]
#' @param measured_criterion apply the measured value criterion
#' @param narrow_ci_correction apply the correction for narrow confidence intervals
#'
#' @return list with the NtC, intermediate values, the fit and the confidence band
#' @export
ntc_original <- function(data, measured_criterion = TRUE, narrow_ci_correction = TRUE) {
  fail <- function(msg) ntc_abort(paste("Original algorithm:", msg), "ntc_algorithm_error")

  conc <- data$conc
  logconc <- log10(conc)
  effect_rep <- as.data.frame(data$effects)
  if (ncol(effect_rep) == 1) effect_rep <- effect_rep[[1]]   # as data[2:n, 2:2] in the script
  effect_mean <- if (is.data.frame(effect_rep)) rowMeans(effect_rep, na.rm = TRUE) else effect_rep

  data_for_fit <- data.frame(concentration = logconc, effect = effect_mean)
  drc_formula <- effect ~ 100/(1+10^((logEC50-concentration)*slope))

  curvefit <- tryCatch({
    fit1 <- nlmrt::nlxb(drc_formula,
                        start = c(logEC50 = log10(conc[ceiling(length(conc)/2)]), slope = -2),
                        trace = FALSE, data = data_for_fit)
    nls2::nls2(drc_formula, data = data_for_fit, start = fit1$coefficients, algorithm = "brute-force")
  }, error = function(e) fail(paste("the curve fit failed:", conditionMessage(e))))

  coefficient <- coef(curvefit)
  logEC50_value <- unname(coefficient[1])
  slope_value <- unname(coefficient[2])
  covmatrix <- tryCatch(vcov(curvefit), error = function(e) NULL)
  if (is.null(covmatrix) || any(!is.finite(covmatrix))) {
    fail("the covariance matrix of the fit could not be computed (the curve is not identifiable from these data).")
  }

  # concentration range: fixed rule of the original method
  if (logconc[1] > 0) {
    x_start <- 0
  } else {
    x_start <- logconc[1] + 2*logconc[1]
  }
  x_end <- logconc[length(logconc)]*2
  if (x_end <= x_start) fail("the fixed concentration range of the original method is empty for these concentrations.")
  x_values <- seq(x_start, x_end, 0.001)

  curve.predict <- predict(curvefit, newdata = data.frame(concentration = c(x_values)))
  pp <- c(logEC50 = logEC50_value, slope = slope_value)
  myjacfun <- nlmrt::model2jacfun(drc_formula, pp)
  myjacobian <- myjacfun(pp, effect = curve.predict, concentration = x_values)
  error <- sqrt(diag(myjacobian %*% covmatrix %*% t(myjacobian)))
  dof <- nrow(data_for_fit) - 2
  if (dof < 1) fail("at least 3 concentrations are needed.")
  t_student <- qt(.975, df = dof)
  CI <- error*t_student
  upperCI <- curve.predict + CI
  lowerCI <- curve.predict - CI
  band <- data.frame(log_conc = x_values, fit = curve.predict, lower = lowerCI, upper = upperCI)

  range_txt <- sprintf("%s to %s", fmt_num(10^x_start), fmt_num(10^x_end))
  ## NtC upper
  up_idx <- which(upperCI >= 99.999999999)
  if (length(up_idx) == 0) {
    fail(sprintf("the upper confidence band never reaches 100%% in the fixed concentration range (%s). The original app fails here.", range_txt))
  }
  inds <- max(up_idx)
  value_atupper <- 100 - curve.predict[inds]
  NtC_upper <- if (value_atupper > 0) 10^(x_values[inds]) else 0

  ## NtC lower
  lo_idx <- which(lowerCI > 90)
  if (length(lo_idx) == 0) {
    fail(sprintf(paste("the lower confidence band never exceeds 90%% survival in the fixed concentration range (%s).",
                       "The original app fails here. If the first row is a control, tick 'First row is a control';",
                       "otherwise the Updated algorithm extends the range."), range_txt))
  }
  inds <- max(lo_idx)
  value_atlower <- 100 - curve.predict[inds]
  NtC_lower <- if (value_atlower > 0) 10^(x_values[inds]) else 0

  ## Choose NtC based on upper and lower CI criteria
  if (NtC_upper > NtC_lower) {
    NtC <- NtC_lower
    Effect_NtC <- value_atlower
    rule <- "lower"
  } else {
    ratio <- value_atlower/value_atupper
    if (is.nan(ratio)) fail("the predicted effect is 0 at both confidence criteria (0/0 in the narrow-CI rule).")
    if (narrow_ci_correction && ratio > 10) {
      Effect_NtC <- value_atlower/10
      NtC <- 10^(logEC50_value-log10(100/(100-Effect_NtC)-1)/slope_value)
      rule <- "narrow_ci_correction"
    } else {
      NtC <- NtC_upper
      Effect_NtC <- value_atupper
      rule <- "upper"
    }
  }
  NtC_model <- NtC

  ## Measured value criterion (script logic, including its index arithmetic)
  measured_applied <- FALSE
  if (measured_criterion && length(which(conc < NtC)) > 0) {
    if (!is.data.frame(effect_rep)) {
      fail("the measured value criterion needs at least two replicate columns (the original script fails with one).")
    }
    sub <- effect_rep[which(conc < NtC), ]
    inds <- which(sub <= 90)
    if (length(inds) > 0) {
      remained <- inds %% nrow(sub)
      remained <- replace(remained, remained == 0, nrow(sub))
      NtC <- conc[min(remained)]
      measured_applied <- TRUE
    }
  }
  toxic <- rowSums(as.matrix(data$effects) <= 90, na.rm = TRUE) > 0
  NtC_measured <- if (any(toxic)) min(conc[toxic]) else NA_real_

  list(NtC = unname(NtC), NtC_model = unname(NtC_model), Effect_NtC = unname(Effect_NtC), rule = rule,
       NtC_lower = unname(NtC_lower), NtC_upper = unname(NtC_upper),
       Effect_NtC_lower = unname(value_atlower), Effect_NtC_upper = unname(value_atupper),
       measured_applied = measured_applied, NtC_measured = NtC_measured,
       logEC50 = logEC50_value, slope = slope_value, EC50 = 10^logEC50_value,
       se = sqrt(diag(covmatrix)), df = dof, t = t_student, level = 0.95,
       fit_data = data.frame(log_conc = logconc, effect = effect_mean),
       fit_to = "concentration means", n_fit = nrow(data_for_fit),
       residuals = unname(residuals(curvefit)), converged = TRUE,
       band = band, range = c(x_start, x_end), n_extensions = 0, warnings = character())
}
