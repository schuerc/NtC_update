ext <- function(f) system.file("extdata", f, package = "NtC")

test_that("all demo file variants give the same data", {
  ref <- ntc_prepare(read_ntc_file(ext("demo_semicolon.csv")))
  for (f in c("demo_comma.csv", "demo_tab.txt", "demo_semicolon_decimal_comma.csv", "demo_unsorted.csv")) {
    d <- ntc_prepare(read_ntc_file(ext(f)))
    expect_equal(d$conc, ref$conc, info = f)
    expect_equal(unname(d$effects), unname(ref$effects), info = f)
  }
})

test_that("separator, decimal mark and header are detected", {
  fmt <- function(lines) parse_ntc_lines(lines)$format
  expect_equal(fmt(c("c;a;b", "1;90;91", "2;50;52", "4;10;12"))[c("separator", "decimal")], list(separator = "semicolon", decimal = "."))
  expect_equal(fmt(c("c;a;b", "1,5;90,1;91", "2;50;52", "4;10;12"))$decimal, ",")
  expect_equal(fmt(c("c,a,b", "1.5,90.1,91", "2,50,52", "4,10,12"))$separator, "comma")
  expect_equal(fmt(c("c\ta\tb", "1\t90\t91", "2\t50\t52"))$separator, "tab")
  expect_equal(fmt(c("1 90 91", "2 50 52", "4 10 12"))$separator, "whitespace")
  expect_false(fmt(c("1 90 91", "2 50 52", "4 10 12"))$header)
  expect_true(fmt(c("﻿conc;a", "1;90", "2;50"))$header)       # byte order mark
})

test_that("control row: automatic for concentration 0, optional otherwise", {
  lines <- c("c;a;b", "0;100;99", "1;95;96", "3;50;55", "10;5;4", "30;1;2")
  d <- ntc_prepare(parse_ntc_lines(lines))
  expect_true(d$control)
  expect_equal(d$conc, c(1, 3, 10, 30))
  expect_error(ntc_prepare(parse_ntc_lines(lines), control = FALSE), "concentration of 0", class = "ntc_data_error")
  d2 <- ntc_prepare(parse_ntc_lines(lines[-2]), control = TRUE)
  expect_equal(d2$conc, c(3, 10, 30))
})

test_that("unsorted concentrations are sorted with a note", {
  d <- ntc_prepare(parse_ntc_lines(c("c;a", "10;5", "1;95", "3;50")))
  expect_equal(d$conc, c(1, 3, 10))
  expect_equal(unname(d$effects[, 1]), c(95, 50, 5))
  expect_match(d$notes, "sorted")
})

test_that("invalid files are rejected with a clear message", {
  expect_error(ntc_prepare(read_ntc_file(ext("invalid_lost_decimal.csv"))), "decimal", class = "ntc_data_error")
  expect_error(parse_ntc_lines(c("c;a", "1,5;90.2", "2;50")), "mixes decimal marks", class = "ntc_parse_error")
  expect_error(parse_ntc_lines(c("c;a", "1;90", "2;fifty")), "separator|not a number", class = "ntc_parse_error")
  expect_error(parse_ntc_lines(c("c;a;b", "1;90", "2;50;40")), "separator", class = "ntc_parse_error")
  expect_error(ntc_prepare(parse_ntc_lines(c("c;a", "1;90", "1;80", "2;50"))), "more than once", class = "ntc_data_error")
  expect_error(ntc_prepare(parse_ntc_lines(c("c;a", "1;90", "2;50"))), "At least 3", class = "ntc_data_error")
  expect_error(parse_ntc_lines("c;a"), "fewer than two", class = "ntc_parse_error")
})

test_that("missing values and column mapping", {
  lines <- c("conc;x;a;b", "1;7;95;NA", "3;7;50;;", "10;7;5;4")
  lines[3] <- "3;7;50;"
  d <- ntc_prepare(parse_ntc_lines(lines), conc_col = "conc", rep_cols = c("a", "b"))
  expect_equal(colnames(d$effects), c("a", "b"))
  expect_true(is.na(d$effects[1, "b"]))
  expect_error(ntc_prepare(parse_ntc_lines(lines), rep_cols = character()), "at least one replicate")
})
