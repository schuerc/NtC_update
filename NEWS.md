# Changes vs. the original app

This page lists every difference between this app and the original NtC app of
Stadnicka-Michalak et al. (2018, ALTEX, doi:10.14573/altex.1701231). The source code of the original
app is not available; the reference is the R script of the method (`inst/reference/StadnickaNtC.R`)
and the published app's help text. Version history: see the end of this page.

## Two algorithms

The app offers both algorithms; the one used is shown under the results, in the downloads and in the
report.

- **Original (Stadnicka-Michalak et al.)**: reproduces the published calculation, including its
  quirks. Where the original would fail or produce a meaningless value, the app shows a message that
  explains why instead of crashing.
- **Updated** (default): based on the corrected implementation in the `bendr` package
  (fork: https://github.com/schuerc/bendr).

## Differences in the calculation

| Aspect | Original | Updated |
|---|---|---|
| Curve fit | `nlxb`, then `nls2` evaluated at that point | Levenberg-Marquardt (`nlsLM`), grid search as fallback |
| Data fitted | Concentration means | Concentration means (replicates on the same plate, default) or all replicate values (independent experiments) |
| Degrees of freedom | number of concentrations - 2 (see below) | same plate: number of concentrations - 2; independent: number of values - 2 |
| Confidence level | 95% | 95% by default, adjustable |
| Concentration range for the CI criteria | Fixed: from concentration 1 (or 3 x log of the lowest concentration if it is <= 1) to 2 x log of the highest | Extended 1 log unit at a time below the tested range until both CI criteria are found |
| Upper CI criterion | upper band >= 99.999999999% | upper band > 100% |
| Criterion not found in the range | Script fails; published app fails | Range extended; if still not found after 10 log units: message, no NtC |
| Upper crossing outside the range | Treated as effect 0, which triggers the narrow-CI correction and gives a too low NtC | Cannot happen (range extended) |
| 0/0 in the narrow-CI rule | Script fails | Message, no NtC |
| Measured value criterion | Index arithmetic of the script; fails with a single replicate column | Same rule, robust with any number of replicates |
| EC50 / EC10 | EC50 only | EC50 and EC10 with confidence intervals (delta method, or bootstrap for independent replicates) |
| Bootstrap | - | Optional, for EC50/EC10/slope only; the NtC always uses the delta-method band, because bootstrap bands never exceed 100% |
| Fit diagnostics | - | Number of points, df, R², residual SD, partial curve, extrapolated NtC, unequal scatter, convergence |

**Why the Updated algorithm can give a different NtC.** On test data, the largest effects come from
(1) the concentration range: the original can miss where the upper band crosses 100%, which then
triggers the narrow-CI correction; (2) the degrees of freedom: treating replicates as independent
narrows the band. A wider band raises NtC_upperCI and lowers NtC_lowerCI, so the NtC can move in
either direction.

## Control row

- The original R script always treats the **first data row** as the control and drops it, even if
  its concentration is not 0. The published app apparently did not: it failed on files without a
  control row for which the script gives a result.
- Here, the checkbox **"First row is a control"** decides this for both algorithms. It is ticked
  automatically when the first concentration is 0. A concentration of 0 that is not excluded is
  rejected (log10(0) cannot be fitted).
- Degrees of freedom of the original algorithm: the script uses `nrow(file) - 3`, which equals
  "number of fitted concentrations - 2" only when one row is dropped as control. The app always uses
  the number of fitted concentrations - 2; with the control box ticked the result equals the script.

## Input files

| Aspect | Original | This app |
|---|---|---|
| Separator | Script: semicolon only; app help text: comma, semicolon or whitespace | Tab, semicolon, comma or whitespace, detected automatically |
| Decimal mark | Dot | Dot or comma, detected automatically |
| Header row | Required | Optional (detected) |
| Column order | First column concentrations, in increasing order | Concentration and replicate columns can be chosen; rows are sorted |
| Implausible values | Not checked | Rejected with a message, e.g. survival far above 100% because the decimal mark was lost on export |
| Mixed decimal marks, text in number cells, duplicate concentrations | Not checked | Rejected with a message naming line, column and value |
| Data preview | - | Shown before the results |

## Output

- Plot with the fitted curve, confidence band, data, NtC, NtC_lowerCI and NtC_upperCI.
- Results table with all intermediate values and the rule that determined the NtC.
- Downloads: plot (PNG, PDF), results (CSV), reproducibility report (HTML: results, settings,
  diagnostics, plot, input data as used, software versions).
- Version of the app, of bendr and of R shown in the footer and in all downloads.

## User interface

- Short hints next to each control; the full explanation is in the **Instructions** tab.
- New controls: algorithm, first row is a control, replicates independent, column mapping, confidence
  level, bootstrap for EC50/EC10.
- Unchanged: measured values criterion and narrow-CI correction (both on by default), concentration units.

# Version history

## NtC 0.1.0 (2026-09-29)

- First version of the rebuilt app: original and updated algorithm, robust file import, data preview,
  fit diagnostics, downloads and report. Uses bendr 0.4.1 (https://github.com/schuerc/bendr).
