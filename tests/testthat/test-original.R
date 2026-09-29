# Regression test: the "original" algorithm must reproduce the reference script
# inst/reference/StadnickaNtC.R. The script always drops the first row, so the comparison
# uses control = TRUE.

run_script <- function(csv) {
  skip_if_not_installed("cowplot")
  script <- readLines(system.file("reference", "StadnickaNtC.R", package = "NtC"))
  script <- sub("PATH_TO_FILE", normalizePath(csv, winslash = "/"), script, fixed = TRUE)
  tmp <- tempfile(fileext = ".R")
  writeLines(script, tmp)
  env <- new.env()
  suppressWarnings(suppressMessages(sys.source(tmp, envir = env)))
  env
}

write_semicolon <- function(df) {
  f <- tempfile(fileext = ".csv")
  utils::write.table(df, f, sep = ";", row.names = FALSE, quote = FALSE)
  f
}

compare_with_script <- function(df, measured = TRUE) {
  csv <- write_semicolon(df)
  env <- run_script(csv)
  d <- ntc_prepare(read_ntc_file(csv), control = TRUE)
  r <- ntc_original(d, measured_criterion = measured)
  expect_equal(r$logEC50, unname(env$logEC50_value), tolerance = 1e-8)
  expect_equal(r$slope, unname(env$slope_value), tolerance = 1e-8)
  expect_equal(r$NtC_upper, unname(env$NtC_upper), tolerance = 1e-10)
  expect_equal(r$NtC_lower, unname(env$NtC_lower), tolerance = 1e-10)
  expect_equal(r$NtC, unname(env$NtC), tolerance = 1e-10)
  expect_equal(r$Effect_NtC, unname(env$Effect_NtC), tolerance = 1e-8)
  r
}

test_that("original matches StadnickaNtC.R on the demo data (control row)", {
  df <- read.csv(system.file("extdata", "demo_semicolon.csv", package = "NtC"), sep = ";", check.names = FALSE)
  compare_with_script(df)
})

test_that("where StadnickaNtC.R fails, the original algorithm fails with a message", {
  # upper band never reaches 100% in the fixed range: the script stops with an error
  df <- read.csv(system.file("extdata", "demo_no_control_below90.csv", package = "NtC"), sep = ";")
  csv <- write_semicolon(df)
  expect_error(run_script(csv))
  d <- ntc_prepare(read_ntc_file(csv), control = TRUE)
  expect_error(ntc_original(d), "upper confidence band never reaches 100%", class = "ntc_algorithm_error")
})

test_that("original matches StadnickaNtC.R when the measured value criterion applies", {
  df <- data.frame(conc = c(0, 0.3, 1, 3, 10, 30), a = c(100, 93, 100, 65, 9, 3),
                   b = c(100, 85, 100, 67, 17, 3), c = c(100, 92, 86, 61, 21, 8))
  r <- compare_with_script(df)
  expect_true(r$measured_applied)
  expect_equal(r$NtC, 0.3)
  expect_gt(r$NtC_model, 0.3)
})

test_that("original: toggles and failures", {
  d <- ntc_prepare(read_ntc_file(system.file("extdata", "demo_semicolon.csv", package = "NtC")))
  r0 <- ntc_original(d, measured_criterion = FALSE, narrow_ci_correction = FALSE)
  expect_true(r0$NtC %in% c(r0$NtC_lower, r0$NtC_upper))
  expect_false(r0$measured_applied)

  one_rep <- ntc_prepare(data.frame(conc = c(0.1, 0.3, 1, 3, 10, 30), a = c(99, 97, 88, 50, 10, 2)))
  expect_error(ntc_original(one_rep), "two replicate columns", class = "ntc_algorithm_error")
  r <- ntc(one_rep, "original")
  expect_equal(r$status, "failed")
  expect_match(r$message, "Original algorithm")
})
