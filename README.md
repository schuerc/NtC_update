# NtC_update

An updated version of the algorithm and Shiny app to calculate the non-toxic concentration (NtC)
according to Stadnicka-Michalak et al. (2018), ALTEX, doi:10.14573/altex.1701231.

The R package is called **NtC**. It offers

- the **original** algorithm, reproduced from the reference script of the method
  (`inst/reference/StadnickaNtC.R`) including its quirks, with clear messages where the original fails;
- an **updated** algorithm based on the corrected [bendr](https://github.com/schuerc/bendr) package:
  concentration range extended until both confidence criteria are found, choice between replicates
  on one plate and independent experiments, EC50/EC10 with confidence intervals, fit diagnostics;
- robust file import (separator, decimal mark and header detected; control row; unsorted
  concentrations; implausible values rejected with a message);
- a Shiny app with data preview, plot, results, diagnostics, downloads (PNG, PDF, CSV) and an HTML
  reproducibility report.

All differences to the original app are listed in [NEWS.md](NEWS.md) (shown in the app's tab
"Changes vs. original app"); the method and all settings are explained in
[inst/text/instructions.md](inst/text/instructions.md) (the app's Instructions tab).

## Installation

```r
# install.packages("remotes")
remotes::install_github("schuerc/NtC_update")
```

This also installs `bendr` from https://github.com/schuerc/bendr.

## Usage

```r
library(NtC)
ntc_app()                                   # start the app

f <- system.file("extdata", "demo_semicolon.csv", package = "NtC")
r <- ntc(f)                                 # updated algorithm, replicates on one plate
r
ntc(f, algorithm = "original")              # original algorithm
ntc(f, replicates = "independent")          # replicates from independent experiments
plot_ntc(r, units = "mg/L")
ntc_report(r, "report.html", units = "mg/L")
```

The files in `inst/extdata` are synthetic example data.

## Input format

First column: concentrations. Further columns: replicates, as survival in % of the control
(100 = no effect). Tab, semicolon, comma or whitespace separated; dot or comma as decimal mark.
A row with concentration 0 is treated as the control. Details: Instructions tab / `inst/text/instructions.md`.

## Background analysis

`analysis/` compares the CI and NtC calculation of the original script and BendR
([ci_comparison_original_vs_bendr.md](analysis/ci_comparison_original_vs_bendr.md)) and decomposes why
they give different NtCs ([ntc_decomposition.md](analysis/ntc_decomposition.md)), on synthetic data.

## Development

`renv.lock` records the package versions the app was tested with (including bendr from GitHub);
`renv::restore()` recreates that library.

```r
devtools::test()    # tests; test-app.R needs shinytest2 and Chrome
```

`tests/testthat/test-local-data.R` compares the original algorithm with the reference script on data
that is not part of the repository: put CSV files into `local-data/` (git-ignored) or set
`NTC_LOCAL_DATA`; the test is skipped otherwise.

## License

GPL-3. The reference script `inst/reference/StadnickaNtC.R` is included for the regression tests.
