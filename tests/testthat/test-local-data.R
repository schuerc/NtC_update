# Checks on confidential data that is NOT part of the repository. Put the files into
# <repo>/local-data/ (git-ignored) or point NTC_LOCAL_DATA to a folder. Skipped otherwise.
# No values or results of these files may be written into this file.

local_files <- function() {
  dir <- Sys.getenv("NTC_LOCAL_DATA", testthat::test_path("..", "..", "local-data"))
  list.files(dir, pattern = "\\.(csv|txt)$", full.names = TRUE)
}

test_that("original algorithm matches StadnickaNtC.R on local data (first row dropped)", {
  files <- local_files()
  skip_if(length(files) == 0, "no local data")
  skip_if_not_installed("cowplot")
  for (csv in files) {
    d <- tryCatch(ntc_prepare(read_ntc_file(csv), control = TRUE), ntc_error = function(e) NULL)
    if (is.null(d)) next
    script <- readLines(system.file("reference", "StadnickaNtC.R", package = "NtC"))
    script <- sub("PATH_TO_FILE", normalizePath(csv, winslash = "/"), script, fixed = TRUE)
    tmp <- tempfile(fileext = ".R")
    writeLines(script, tmp)
    env <- new.env()
    ok <- tryCatch({ suppressWarnings(suppressMessages(sys.source(tmp, envir = env))); TRUE }, error = function(e) FALSE)
    if (!ok || !exists("NtC", envir = env)) next    # e.g. comma-separated files the script cannot read
    r <- ntc_original(d)
    expect_equal(r$NtC, unname(env$NtC), tolerance = 1e-10, info = basename(csv))
  }
})
