#' @importFrom stats coef vcov predict qt fitted residuals sd setNames
#' @importFrom utils packageVersion
#' @importFrom markdown mark_html
NULL

utils::globalVariables(c("conc", "log_conc", "fit", "lower", "upper", "effect", "replicate", "x", "label"))

# Signal an error that the app shows to the user as a message (not as a crash).
ntc_abort <- function(message, class = NULL) {
  stop(structure(class = c(class, "ntc_error", "error", "condition"),
                 list(message = message, call = NULL)))
}

# Model: survival (100% - effect) as a function of log10 concentration
drc_formula <- function() effect ~ 100 / (1 + 10^((logEC50 - logconc) * slope))

hill <- function(logconc, logEC50, slope) 100 / (1 + 10^((logEC50 - logconc) * slope))

`%||%` <- function(a, b) if (is.null(a)) b else a

fmt_num <- function(x, digits = 4) {
  ifelse(is.na(x), "NA", formatC(signif(x, digits), digits = digits, format = "fg", flag = "#"))
}

#' Versions of the NtC package, bendr and R
#'
#' @return named character vector
#' @export
ntc_versions <- function() {
  desc <- utils::packageDescription("NtC")
  c(NtC = as.character(packageVersion("NtC")),
    NtC_date = desc$Date %||% NA_character_,
    bendr = as.character(packageVersion("bendr")),
    R = paste(R.version$major, R.version$minor, sep = "."))
}
