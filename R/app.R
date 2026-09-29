#' Run the NtC Shiny app
#'
#' @param ... passed to [shiny::shinyApp()] options, e.g. `launch.browser`
#' @return a shiny app object
#' @export
ntc_app <- function(...) {
  shiny::shinyApp(ntc_ui(), ntc_server, options = list(...))
}

text_file <- function(name) system.file("text", name, package = "NtC")

hint <- function(...) shiny::helpText(..., style = "margin-top:-6px;margin-bottom:12px;font-size:12.5px")

ntc_ui <- function() {
  v <- ntc_versions()
  shiny::navbarPage(
    title = "Determination of chemical Non-toxic Concentrations (NtC)",
    id = "tabs", collapsible = TRUE,
    header = shiny::tags$head(shiny::tags$style(shiny::HTML(paste(
      ".msg-error{background:#fdecea;border-left:4px solid #c62828;padding:8px 12px;margin-bottom:12px}",
      ".msg-note{background:#eef5fb;border-left:4px solid #2c7fb8;padding:8px 12px;margin-bottom:12px}",
      ".diag-warning td{color:#8a4b00}",
      ".footer{color:#777;font-size:12px;margin-top:24px;border-top:1px solid #eee;padding-top:8px}",
      "table{font-size:13.5px}")))),
    shiny::tabPanel("NtC",
      shiny::sidebarLayout(
        shiny::sidebarPanel(width = 4,
          shiny::fileInput("file", "Choose CSV file", accept = c(".csv", ".txt", ".tsv", "text/csv", "text/plain")),
          hint("First column: concentrations; other columns: replicates (survival in %). Separator and decimal mark are detected. (Instructions \u00a71)"),
          shiny::actionLink("demo", "Or load the synthetic example data"),
          shiny::tags$br(), shiny::tags$br(),
          shiny::uiOutput("mapping"),
          shiny::checkboxInput("control", "First row is a control (not used for the fit)", FALSE),
          hint("Ticked automatically if the first concentration is 0. (\u00a71)"),
          shiny::radioButtons("algorithm", "Algorithm",
                              c("Updated" = "updated", "Original (Stadnicka-Michalak et al.)" = "original")),
          hint("Original = published app incl. its quirks; Updated = corrected version. (\u00a77)"),
          shiny::conditionalPanel("input.algorithm == 'updated'",
            shiny::checkboxInput("independent", "Replicates are independent experiments", FALSE),
            hint("Leave unticked for wells on one plate. Changes the confidence band and the NtC. (\u00a76)"),
            shiny::tags$details(shiny::tags$summary("Advanced settings"),
              shiny::numericInput("level", "Confidence level", 0.95, min = 0.8, max = 0.99, step = 0.01),
              shiny::conditionalPanel("input.independent",
                shiny::checkboxInput("bootstrap", "Bootstrap CIs for EC50/EC10 (1000 resamples, slower)", FALSE)),
              hint("The NtC always uses the delta-method band. (\u00a73)"))),
          shiny::conditionalPanel("input.algorithm == 'original'",
            hint("Fixed settings of the original method: 95% confidence band, fit to concentration means, df = number of concentrations - 2, fixed concentration range.")),
          shiny::checkboxInput("measured", "Apply measured values criterion", TRUE),
          hint("NtC \u2264 lowest tested concentration with \u2265 10% effect in any replicate. (\u00a75)"),
          shiny::checkboxInput("narrow", "Apply correction for narrow confidence intervals", TRUE),
          hint("Uses 1/10 of the lower-CI effect when the CI is very narrow. (\u00a74)"),
          shiny::textInput("units", "Concentration units", "units"),
          shiny::tags$hr(),
          shiny::helpText("Method and settings are explained in the Instructions tab. See Stadnicka-Michalak et al.,",
                          shiny::tags$a(href = "https://www.ncbi.nlm.nih.gov/labs/articles/28653737/", target = "_blank",
                                        "https://www.ncbi.nlm.nih.gov/labs/articles/28653737/"))),
        shiny::mainPanel(width = 8,
          shiny::uiOutput("messages"),
          shiny::uiOutput("results_ui"),
          shiny::uiOutput("preview_ui"),
          shiny::div(class = "footer", shiny::textOutput("footer", inline = TRUE))))),
    shiny::tabPanel("Instructions", shiny::div(style = "max-width:900px", shiny::includeMarkdown(text_file("instructions.md")))),
    shiny::tabPanel("Changes vs. original app", shiny::div(style = "max-width:900px", shiny::includeMarkdown(system.file("NEWS.md", package = "NtC")))),
    footer = shiny::div(class = "container-fluid footer",
                        sprintf("NtC %s (%s) \u00b7 bendr %s \u00b7 R %s", v[["NtC"]], v[["NtC_date"]], v[["bendr"]], v[["R"]])))
}

ntc_server <- function(input, output, session) {
  as_message <- function(e) if (inherits(e, "ntc_error")) conditionMessage(e) else paste("Unexpected error:", conditionMessage(e))

  source_file <- shiny::reactiveVal(NULL)
  shiny::observeEvent(input$file, source_file(list(path = input$file$datapath, name = input$file$name)))
  shiny::observeEvent(input$demo, source_file(list(path = system.file("extdata", "demo_semicolon.csv", package = "NtC"),
                                                  name = "demo_semicolon.csv (synthetic)")))

  raw <- shiny::reactive({
    src <- source_file()
    shiny::req(src)
    tryCatch(read_ntc_file(src$path, name = src$name), error = function(e) structure(list(error = as_message(e)), class = "ntc_failed"))
  })

  # column mapping and control-row default when a new file is loaded
  output$mapping <- shiny::renderUI({
    r <- raw()
    if (inherits(r, "ntc_failed")) return(NULL)
    cols <- names(r$data)
    shiny::tagList(
      shiny::selectInput("conc_col", "Concentration column", cols, cols[1]),
      shiny::checkboxGroupInput("rep_cols", "Replicate columns", cols[-1], cols[-1], inline = TRUE))
  })
  shiny::observeEvent(list(raw(), input$conc_col), {
    r <- raw()
    if (inherits(r, "ntc_failed")) return()
    cc <- input$conc_col %||% names(r$data)[1]
    if (!cc %in% names(r$data)) return()
    conc <- r$data[[cc]]
    shiny::updateCheckboxInput(session, "control", value = isTRUE(min(conc, na.rm = TRUE) == 0))
  })

  prepared <- shiny::reactive({
    r <- raw()
    if (inherits(r, "ntc_failed")) return(r)
    shiny::req(input$conc_col)
    tryCatch(ntc_prepare(r, conc_col = input$conc_col, rep_cols = input$rep_cols, control = input$control),
             error = function(e) structure(list(error = as_message(e)), class = "ntc_failed"))
  })

  result <- shiny::reactive({
    d <- prepared()
    if (inherits(d, "ntc_failed")) return(d)
    boot <- isTRUE(input$independent) && isTRUE(input$bootstrap)
    level <- input$level
    if (!is.numeric(level) || is.na(level) || level < 0.8 || level > 0.99) level <- 0.95
    shiny::withProgress(message = if (boot) "Fitting (bootstrap, may take a minute)" else "Fitting", value = 0.5, {
      tryCatch(ntc(d, algorithm = input$algorithm, measured_criterion = input$measured,
                   narrow_ci_correction = input$narrow,
                   replicates = if (isTRUE(input$independent)) "independent" else "same_plate",
                   level = level, ec_ci = if (boot) "bootstrap" else "delta"),
               error = function(e) structure(list(error = as_message(e)), class = "ntc_failed"))
    })
  })

  output$messages <- shiny::renderUI({
    if (is.null(source_file())) {
      return(shiny::div(class = "msg-note", "Upload a CSV file or load the synthetic example data to start."))
    }
    r <- result()
    if (inherits(r, "ntc_failed")) return(shiny::div(class = "msg-error", shiny::strong("The file cannot be used: "), r$error))
    if (r$status != "ok") return(shiny::div(class = "msg-error", shiny::strong("No NtC: "), r$message))
    NULL
  })

  ok_result <- shiny::reactive({
    r <- result()
    shiny::req(!inherits(r, "ntc_failed"))
    r
  })

  output$results_ui <- shiny::renderUI({
    r <- ok_result()
    shiny::tagList(
      if (r$status == "ok") shiny::h3(sprintf("NtC = %s %s", fmt_num(r$fit$NtC), input$units), style = "margin-top:0"),
      shiny::plotOutput("plot", height = "430px"),
      shiny::fluidRow(
        shiny::column(6, shiny::h4("Results"), shiny::tableOutput("results")),
        shiny::column(6, shiny::h4("Fit diagnostics"), shiny::tableOutput("diagnostics"))),
      shiny::div(style = "margin:6px 0 18px",
        shiny::downloadButton("dl_png", "Plot (PNG)"), shiny::downloadButton("dl_pdf", "Plot (PDF)"),
        shiny::downloadButton("dl_csv", "Results (CSV)"), shiny::downloadButton("dl_report", "Report (HTML)")))
  })
  output$plot <- shiny::renderPlot(plot_ntc(ok_result(), input$units), res = 96)
  output$results <- shiny::renderTable(ntc_table(ok_result())[-1, ], striped = TRUE, na = "")
  output$diagnostics <- shiny::renderTable({
    d <- ok_result()$diagnostics
    data.frame(` ` = ifelse(d$level == "warning", "\u26a0", "\u2139"), message = d$message, check.names = FALSE)
  }, striped = TRUE)

  output$preview_ui <- shiny::renderUI({
    d <- prepared()
    if (inherits(d, "ntc_failed")) return(NULL)
    f <- d$format
    shiny::tagList(
      shiny::h4("Data preview"),
      shiny::helpText(sprintf("%s: separator %s, decimal mark '%s', %s. %s", d$file, f$separator %||% "-", f$decimal %||% "-",
                              if (isTRUE(f$header)) "header row" else "no header row", paste(d$notes, collapse = " "))),
      shiny::tableOutput("preview"))
  })
  output$preview <- shiny::renderTable(as.data.frame(prepared()), digits = 3)

  output$footer <- shiny::renderText({
    r <- result()
    if (inherits(r, "ntc_failed") || is.null(r)) return("")
    sprintf("Calculated with the %s algorithm, NtC %s, bendr %s, %s.", r$settings$algorithm, r$versions[["NtC"]], r$versions[["bendr"]], r$time)
  })

  stem <- function() sprintf("NtC_%s_%s", tools::file_path_sans_ext(basename(ok_result()$data$file)), ok_result()$settings$algorithm)
  output$dl_png <- shiny::downloadHandler(function() paste0(stem(), ".png"),
    function(file) ggplot2::ggsave(file, plot_ntc(ok_result(), input$units), width = 8, height = 5, dpi = 150))
  output$dl_pdf <- shiny::downloadHandler(function() paste0(stem(), ".pdf"),
    function(file) ggplot2::ggsave(file, plot_ntc(ok_result(), input$units), width = 8, height = 5))
  output$dl_csv <- shiny::downloadHandler(function() paste0(stem(), ".csv"),
    function(file) utils::write.csv(ntc_table(ok_result()), file, row.names = FALSE))
  output$dl_report <- shiny::downloadHandler(function() paste0(stem(), "_report.html"),
    function(file) ntc_report(ok_result(), file, input$units))
}
