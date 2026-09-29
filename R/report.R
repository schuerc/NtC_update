#' Write a self-contained HTML reproducibility report
#'
#' Contains the results, all settings, algorithm and package versions, fit diagnostics,
#' the plot and the input data as used.
#'
#' @param result an `ntc_result`
#' @param file output path (.html)
#' @param units concentration units
#' @return the path, invisibly
#' @export
ntc_report <- function(result, file, units = "units") {
  png <- tempfile(fileext = ".png")
  on.exit(unlink(png))
  ggplot2::ggsave(png, plot_ntc(result, units), width = 8, height = 5, dpi = 110)
  img <- base64enc::dataURI(file = png, mime = "image/png")

  tbl <- function(df) {
    htmltools::tags$table(class = "t",
      htmltools::tags$tr(lapply(names(df), htmltools::tags$th)),
      lapply(seq_len(nrow(df)), function(i) htmltools::tags$tr(lapply(df[i, ], function(v) htmltools::tags$td(as.character(v))))))
  }
  s <- result$settings
  settings <- data.frame(
    setting = c("Algorithm", "First row is a control", "Measured value criterion", "Narrow-CI correction",
                "Replicates", "Confidence level", "EC50/EC10 CIs", "Concentration units"),
    value = c(s$algorithm, s$control, s$measured_criterion, s$narrow_ci_correction,
              s$replicates, s$level, s$ec_ci, units))
  fmt <- result$data$format
  input <- data.frame(item = c("File", "Separator", "Decimal mark", "Header row"),
                      value = c(fmt$file %||% "", fmt$separator %||% "", fmt$decimal %||% "", as.character(fmt$header %||% "")))
  pkgs <- c("NtC", "bendr", "minpack.lm", "nlmrt", "nls2", "nlstools", "ggplot2", "shiny")
  versions <- data.frame(component = c("R", pkgs),
                         version = c(result$versions[["R"]], vapply(pkgs, function(p) tryCatch(as.character(packageVersion(p)), error = function(e) "not installed"), "")))

  doc <- htmltools::tags$html(
    htmltools::tags$head(htmltools::tags$meta(charset = "utf-8"), htmltools::tags$title("NtC report"),
      htmltools::tags$style(paste(
        "body{font-family:Segoe UI,Arial,sans-serif;max-width:900px;margin:2em auto;padding:0 1em;color:#222}",
        "table.t{border-collapse:collapse;margin:.5em 0 1.5em}",
        ".t th,.t td{border:1px solid #ccc;padding:3px 8px;text-align:left;font-size:14px}",
        ".t th{background:#f2f2f2}.fail{color:#b30000;font-weight:bold}img{max-width:100%}"))),
    htmltools::tags$body(
      htmltools::tags$h1("Non-toxic concentration (NtC) report"),
      htmltools::tags$p(sprintf("Created %s with NtC %s (%s), bendr %s, R %s.", result$time, result$versions[["NtC"]],
                                result$versions[["NtC_date"]], result$versions[["bendr"]], result$versions[["R"]])),
      htmltools::tags$p("Method: Stadnicka-Michalak et al. (2018), ALTEX, doi:10.14573/altex.1701231."),
      if (result$status != "ok") htmltools::tags$p(class = "fail", result$message),
      htmltools::tags$h2("Results"), tbl(ntc_table(result)),
      htmltools::tags$img(src = img),
      htmltools::tags$h2("Settings"), tbl(settings),
      htmltools::tags$h2("Fit diagnostics"), tbl(result$diagnostics),
      htmltools::tags$h2("Input data as used"), tbl(input), tbl(as.data.frame(result$data)),
      htmltools::tags$h2("Software versions"), tbl(versions)))
  writeLines(c("<!DOCTYPE html>", as.character(doc)), file, useBytes = TRUE)
  invisible(file)
}
