# Reading and validating NtC input files ------------------------------------

# number patterns (decimal dot / decimal comma), optional exponent
re_int <- "^[-+]?[0-9]+$"
re_dot <- "^[-+]?([0-9]+\\.?[0-9]*|\\.[0-9]+)([eE][-+]?[0-9]+)?$"
re_comma <- "^[-+]?([0-9]+,?[0-9]*|,[0-9]+)([eE][-+]?[0-9]+)?$"
na_strings <- c("", "NA", "na", "N/A", "n/a", "-", "NaN", "nan")

separators <- c(tab = "\t", semicolon = ";", comma = ",", whitespace = " ")

split_line <- function(line, sep) {
  if (sep == " ") {
    f <- strsplit(trimws(line), "[ \t]+")[[1]]
  } else {
    f <- strsplit(line, sep, fixed = TRUE)[[1]]
    # strsplit drops one trailing empty field; keep the count consistent
    if (endsWith(line, sep)) f <- c(f, "")
  }
  f <- trimws(f)
  gsub('^"(.*)"$', "\\1", f)
}

is_number_like <- function(f, allow_comma = TRUE) {
  f %in% na_strings | grepl(re_dot, f) | (allow_comma & grepl(re_comma, f))
}

# Pick the separator that splits every line into the same number (>= 2) of fields
# and makes all data fields number-like. Order: tab, semicolon, comma, whitespace.
detect_separator <- function(lines) {
  for (nm in names(separators)) {
    sep <- separators[[nm]]
    fields <- lapply(lines, split_line, sep = sep)
    n <- lengths(fields)
    if (any(n != n[1]) || n[1] < 2) next
    # the first line may be a header
    data_fields <- unlist(fields[-1])
    if (!all(is_number_like(data_fields, allow_comma = sep != ","))) next
    return(list(name = nm, sep = sep, fields = fields))
  }
  ntc_abort(paste(
    "The column separator could not be detected. Every line must have the same number of",
    "columns, separated by tab, semicolon, comma or spaces, and all values below the header",
    "must be numbers."), "ntc_parse_error")
}

#' Read an NtC input file
#'
#' Detects the column separator (tab, semicolon, comma, whitespace), the decimal mark
#' (dot or comma) and whether the first line is a header. Values that cannot be read
#' unambiguously are rejected with an explanatory error.
#'
#' @param path path to the file
#' @param name file name used in messages (default: basename of path)
#'
#' @return a list with `data` (numeric data frame, columns as in the file) and `format`
#' @export
read_ntc_file <- function(path, name = basename(path)) {
  if (!file.exists(path)) ntc_abort(sprintf("File not found: %s", path), "ntc_parse_error")
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  parse_ntc_lines(lines, name = name)
}

#' @rdname read_ntc_file
#' @param lines character vector, one element per line
#' @export
parse_ntc_lines <- function(lines, name = "input") {
  lines <- sub("^\ufeff", "", lines)          # byte order mark
  lines <- sub("[ \t\r]+$", "", lines)
  lines <- lines[nzchar(trimws(lines))]
  if (length(lines) < 2) {
    ntc_abort(sprintf("%s: the file has fewer than two non-empty lines.", name), "ntc_parse_error")
  }
  det <- detect_separator(lines)
  fields <- det$fields

  first_numeric <- all(is_number_like(fields[[1]], allow_comma = det$sep != ","))
  has_header <- !first_numeric
  header <- if (has_header) fields[[1]] else c("conc", paste0("rep", seq_len(length(fields[[1]]) - 1)))
  header[header == ""] <- paste0("V", which(header == ""))
  header <- make.unique(header, sep = "_")
  body <- if (has_header) fields[-1] else fields
  line_no <- seq_along(body) + has_header
  m <- do.call(rbind, body)

  # decimal mark
  vals <- m[!(m %in% na_strings)]
  dec_comma <- det$sep != "," && any(grepl(",", vals, fixed = TRUE))
  dec_dot <- any(grepl(".", vals, fixed = TRUE))
  if (dec_comma && dec_dot) {
    ex_c <- vals[grepl(",", vals, fixed = TRUE)][1]
    ex_d <- vals[grepl(".", vals, fixed = TRUE)][1]
    ntc_abort(sprintf(paste(
      "%s: the file mixes decimal marks (e.g. '%s' and '%s'). Use either a dot or a comma",
      "as decimal mark throughout the file."), name, ex_c, ex_d), "ntc_parse_error")
  }
  dec <- if (dec_comma) "," else "."
  if (dec_comma) m <- gsub(",", ".", m, fixed = TRUE)

  num <- suppressWarnings(matrix(as.numeric(m), nrow = nrow(m)))
  num[m %in% na_strings] <- NA
  bad <- which(is.na(num) & !(m %in% na_strings), arr.ind = TRUE)
  if (nrow(bad) > 0) {
    ntc_abort(sprintf("%s, line %d, column '%s': '%s' is not a number.", name,
                      line_no[bad[1, 1]], header[bad[1, 2]], body[[bad[1, 1]]][bad[1, 2]]),
              "ntc_parse_error")
  }
  data <- as.data.frame(num)
  names(data) <- header
  list(data = data,
       format = list(file = name, separator = det$name, decimal = dec, header = has_header,
                     n_lines = length(line_no), line_numbers = line_no))
}

#' Prepare data for the NtC calculation
#'
#' Selects the concentration and replicate columns, sorts by concentration, handles the
#' control row and validates the values.
#'
#' @param x result of [read_ntc_file()], or a data frame
#' @param conc_col name or index of the concentration column (default: first column)
#' @param rep_cols names or indices of the replicate columns (default: all other columns)
#' @param control is the first row (lowest concentration after sorting) a control that is not
#'   used for the fit? `"auto"` (default): yes if its concentration is 0.
#'
#' @return an object of class `ntc_data`
#' @export
ntc_prepare <- function(x, conc_col = 1, rep_cols = NULL, control = "auto") {
  if (is.data.frame(x)) x <- list(data = x, format = list(file = "data frame"))
  raw <- x$data
  name <- x$format$file %||% "input"
  cols <- names(raw)
  conc_col <- if (is.numeric(conc_col)) cols[conc_col] else conc_col
  if (is.null(rep_cols)) rep_cols <- setdiff(cols, conc_col)
  if (is.numeric(rep_cols)) rep_cols <- cols[rep_cols]
  if (!conc_col %in% cols) ntc_abort(sprintf("Concentration column '%s' not found.", conc_col), "ntc_data_error")
  rep_cols <- setdiff(rep_cols, conc_col)
  if (length(rep_cols) == 0) ntc_abort("Select at least one replicate column.", "ntc_data_error")
  if (!all(rep_cols %in% cols)) {
    ntc_abort(sprintf("Replicate column(s) not found: %s", paste(setdiff(rep_cols, cols), collapse = ", ")), "ntc_data_error")
  }

  conc <- raw[[conc_col]]
  eff <- as.matrix(raw[, rep_cols, drop = FALSE])
  colnames(eff) <- rep_cols
  notes <- character()

  line_of <- function(i) if (!is.null(x$format$line_numbers)) sprintf("line %d", x$format$line_numbers[i]) else sprintf("row %d", i)
  if (anyNA(conc)) ntc_abort(sprintf("%s, %s: the concentration is missing.", name, line_of(which(is.na(conc))[1])), "ntc_data_error")
  if (any(conc < 0)) ntc_abort(sprintf("%s, %s: negative concentration (%s).", name, line_of(which(conc < 0)[1]), conc[conc < 0][1]), "ntc_data_error")
  if (anyDuplicated(conc)) {
    ntc_abort(sprintf("%s: concentration %s occurs more than once. Use one row per concentration and one column per replicate.",
                      name, conc[duplicated(conc)][1]), "ntc_data_error")
  }
  out <- which(!is.na(eff) & (eff > 150 | eff < -50), arr.ind = TRUE)
  if (nrow(out) > 0) {
    ntc_abort(sprintf(paste(
      "%s, %s, column '%s': %s is outside the plausible range for survival in %% (-50 to 150).",
      "Values must be survival in %% of the control (100 = no effect). A common cause is a decimal",
      "mark lost when exporting from a spreadsheet (e.g. 81196 instead of 81,196)."),
      name, line_of(out[1, 1]), rep_cols[out[1, 2]], eff[out[1, 1], out[1, 2]]), "ntc_data_error")
  }
  if (any(colSums(!is.na(eff)) == 0)) {
    ntc_abort(sprintf("Replicate column '%s' contains no values.", rep_cols[colSums(!is.na(eff)) == 0][1]), "ntc_data_error")
  }

  ord <- order(conc)
  sorted <- is.unsorted(conc)
  if (sorted) notes <- c(notes, "Concentrations were not in increasing order and have been sorted.")
  conc <- conc[ord]
  eff <- eff[ord, , drop = FALSE]

  if (identical(control, "auto")) control <- conc[1] == 0
  control <- isTRUE(control)
  control_row <- NULL
  if (control) {
    control_row <- list(conc = conc[1], effects = eff[1, ])
    notes <- c(notes, sprintf("First row (concentration %s) treated as control and not used for the fit.", conc[1]))
    conc <- conc[-1]
    eff <- eff[-1, , drop = FALSE]
  }
  if (any(conc == 0)) {
    ntc_abort(paste(
      "A concentration of 0 cannot be used on the log scale. If the row with concentration 0 is the",
      "control, tick 'First row is a control'."), "ntc_data_error")
  }
  empty <- rowSums(!is.na(eff)) == 0
  if (any(empty)) {
    notes <- c(notes, sprintf("Concentration(s) without any value ignored: %s.", paste(conc[empty], collapse = ", ")))
    conc <- conc[!empty]
    eff <- eff[!empty, , drop = FALSE]
  }
  if (length(conc) < 3) {
    ntc_abort(sprintf("At least 3 concentrations with values are needed for the fit (found %d).", length(conc)), "ntc_data_error")
  }

  structure(list(conc = conc, effects = eff, control = control, control_row = control_row,
                 conc_col = conc_col, rep_cols = rep_cols, format = x$format, notes = notes,
                 file = name),
            class = "ntc_data")
}

#' @export
as.data.frame.ntc_data <- function(x, ...) {
  d <- data.frame(concentration = x$conc, x$effects, check.names = FALSE)
  if (!is.null(x$control_row)) {
    ctrl <- data.frame(concentration = x$control_row$conc, t(x$control_row$effects), check.names = FALSE)
    names(ctrl) <- names(d)
    d <- rbind(ctrl, d)
    d <- cbind(role = c("control", rep("data", length(x$conc))), d)
  }
  d
}

#' @export
print.ntc_data <- function(x, ...) {
  cat(sprintf("NtC input: %s, %d concentrations x %d replicate column(s)%s\n", x$file, length(x$conc),
              ncol(x$effects), if (x$control) ", control row excluded" else ""))
  print(as.data.frame(x), row.names = FALSE)
  if (length(x$notes)) cat(paste0("Note: ", x$notes, collapse = "\n"), "\n")
  invisible(x)
}
