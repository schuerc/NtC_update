demo <- function(f = "demo_semicolon.csv") ntc_prepare(read_ntc_file(system.file("extdata", f, package = "NtC")))

test_that("updated algorithm on the demo data", {
  r <- ntc(demo())
  expect_equal(r$status, "ok")
  expect_true(r$fit$NtC > 0 && r$fit$NtC < 3)
  expect_equal(r$fit$df, 5)                 # 7 concentration means - 2
  expect_equal(r$fit$fit_to, "concentration means")
  ri <- ntc(demo(), replicates = "independent")
  expect_equal(ri$fit$df, 19)               # 21 values - 2
  expect_false(isTRUE(all.equal(r$fit$NtC, ri$fit$NtC)))
})

test_that("survival never above 90%: the range is extended and an NtC is found", {
  r <- ntc(demo("demo_no_control_below90.csv"))
  expect_equal(r$status, "ok")
  expect_true(is.finite(r$fit$NtC) && r$fit$NtC < 0.5)
  expect_true(any(grepl("extrapolation", r$diagnostics$message)))
})

test_that("measured value criterion with one concentration below the NtC", {
  # only the lowest concentration lies below the model NtC; one replicate there shows > 10% effect
  df <- data.frame(conc = c(1, 3, 10, 30, 100), a = c(99, 97, 60, 10, 1), b = c(100, 85, 55, 12, 2))
  r_off <- ntc(df, measured_criterion = FALSE)
  r_on <- ntc(df)
  expect_equal(r_on$status, "ok")
  if (sum(df$conc < r_off$fit$NtC) >= 1 && any(df[df$conc < r_off$fit$NtC, -1] <= 90)) {
    expect_true(r_on$fit$measured_applied)
    expect_equal(r_on$fit$NtC, min(df$conc[df$conc < r_off$fit$NtC & apply(df[, -1] <= 90, 1, any)]))
  } else {
    expect_equal(r_on$fit$NtC, r_off$fit$NtC)
  }
})

test_that("flat curve and non-decreasing data fail with a message", {
  flat <- data.frame(conc = c(1, 3, 10, 30, 100), a = c(51, 49, 50, 52, 48), b = c(50, 52, 49, 50, 51))
  r <- ntc(flat)
  expect_equal(r$status, "failed")
  expect_match(r$message, "Updated algorithm")
  incr <- data.frame(conc = c(1, 3, 10, 30, 100), a = c(5, 20, 50, 80, 95), b = c(8, 25, 45, 85, 97))
  expect_match(ntc(incr)$message, "does not decrease")
})

test_that("non-convergence and noise do not crash", {
  noise <- data.frame(conc = c(1, 3, 10, 30, 100), a = c(20, 95, 5, 80, 40), b = c(60, 10, 90, 30, 70))
  for (alg in c("updated", "original")) {
    r <- ntc(noise, alg)
    expect_equal(r$status, "failed")
    expect_true(nzchar(r$message))
  }
})

test_that("settings: toggles, level, bootstrap", {
  d <- demo()
  r0 <- ntc(d, measured_criterion = FALSE, narrow_ci_correction = FALSE)
  expect_true(r0$fit$NtC %in% c(r0$fit$NtC_lower, r0$fit$NtC_upper))
  r90 <- ntc(d, level = 0.90)
  expect_equal(r90$fit$level, 0.90)
  expect_error(ntc_updated(d, ec_ci = "bootstrap"), "independent replicates")
  set.seed(1)
  rb <- ntc(d, replicates = "independent", ec_ci = "bootstrap")
  ri <- ntc(d, replicates = "independent")
  expect_equal(rb$fit$NtC, ri$fit$NtC)      # NtC always from the delta-method band
  expect_false(isTRUE(all.equal(rb$fit$EC50_ci, ri$fit$EC50_ci)))
})

test_that("outputs: table, plot, report", {
  r <- ntc(demo())
  expect_s3_class(plot_ntc(r, "mg/L"), "ggplot")
  expect_true("NtC" %in% ntc_table(r)$quantity)
  f <- tempfile(fileext = ".html")
  ntc_report(r, f, "mg/L")
  html <- paste(readLines(f, encoding = "UTF-8"), collapse = "\n")
  expect_match(html, "data:image/png;base64")
  expect_match(html, "Software versions")
  failed <- ntc(data.frame(conc = c(1, 3, 10, 30, 100), a = c(51, 49, 50, 52, 48), b = c(50, 52, 49, 50, 51)))
  expect_s3_class(plot_ntc(failed), "ggplot")
  expect_no_error(ntc_report(failed, tempfile(fileext = ".html")))
})
