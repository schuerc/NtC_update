# Instructions

This tab holds the full description of the method and of all settings. The main page only shows
short hints. Method reference: Stadnicka-Michalak et al., https://www.ncbi.nlm.nih.gov/labs/articles/28653737/

## 1. File format

The uploaded file should be a .csv file.

- **First column:** concentrations.
- **Other columns:** replicate measurements, as survival in % of the control (100% = no effect,
  i.e. the values are *100% − effect*).
- A header row with column names is expected.

Handled automatically (Updated algorithm; the data preview shows how the file was read):

- **Separator:** comma, semicolon, tab or whitespace.
- **Decimal mark:** dot or comma (e.g. `1.25` or `1,25`).
- **Order:** concentrations do not need to be sorted; they are sorted increasing.
- **Control row:** see "First row is a control" below.
- **Missing values:** empty cells are ignored.

Files are rejected with a message when values cannot be read unambiguously, for example:

- survival values far above 100% (e.g. `81196` instead of `81,196`: the decimal mark was lost on
  export from a spreadsheet);
- different decimal conventions within one file;
- non-numeric entries, fewer than 3 concentrations, or no replicate column.

**First row is a control (checkbox, both algorithms):** when ticked, the first row (lowest
concentration) is treated as the control and is used neither for the fit nor for the measured value
criterion. It is ticked automatically when the first concentration is 0, and unticked otherwise; you
can change it, e.g. when your control is listed with a small non-zero concentration. A concentration of
0 cannot be fitted on the log scale, so a file with a 0 concentration and the box unticked is rejected.

Note on the original method: the original R script of the method (`StadnickaNtC.R`) always drops
the first row, whether or not its concentration is 0. The published app apparently did not. For a
file without a control row this can decide whether the Original algorithm finds an NtC (box ticked,
as the script) or fails (box unticked: the lower band may never exceed 90% within the fixed
concentration range of the original method, as happened in the published app).

## 2. Algorithm description

The fitted sigmoidal curve is plotted together with 95% confidence intervals and the measured data.
To determine the NtC, the following parameters are calculated:

- **NtC_lowerCI:** the highest fitted chemical concentration that is lower than the concentration at
  the intersection of the lower confidence interval and the 10% effect (90% survival).
- **NtC_upperCI:** the highest fitted chemical concentration that is lower than the concentration at
  the intersection of the upper confidence interval and the 0% effect (100% survival).
- **NtC_measured:** the lowest tested chemical concentration that caused at least 10% effect in any
  of the biological replicates.

The NtC is chosen as follows:

1. If NtC_upperCI > NtC_lowerCI, the NtC is NtC_lowerCI.
2. Otherwise the NtC is NtC_upperCI, unless the narrow confidence interval correction applies (section 4).
3. If the measured value criterion is switched on and NtC_measured is below this NtC, the NtC is
   NtC_measured (section 5).

The confidence band is evaluated on a fine concentration grid (0.001 log10 units). In the Updated
algorithm the grid is extended below the tested range until both intersections are found; the
Original algorithm uses a fixed range and may miss an intersection (see "Changes vs. original app").

## 3. Sigmoidal dose-response fit

The following equation is fitted to the measured data:

    100% − effect = 100 / (1 + 10^((logEC50 − log10(concentration)) · slope))

Top and bottom are fixed at 100% and 0%. The two fitted parameters are logEC50 and the slope
(negative, because survival decreases with concentration).

The 95% confidence band is a pointwise band on the fitted curve, calculated with the delta method:
the parameter uncertainty (covariance matrix of the fit) is propagated to the curve and multiplied
by the t-quantile for the degrees of freedom (section 6).

Updated algorithm, optional: **bootstrap** confidence intervals (1000 resamples of the replicates
within each concentration). Only available for independent replicates.

## 4. Narrow confidence interval correction

In case of very narrow confidence intervals, caused e.g. by a non-optimal concentration range:

If Effect_NtC_lowerCI > 10 · Effect_NtC_upperCI, then Effect_NtC = Effect_NtC_lowerCI / 10,
where Effect_NtC_i (%) is the effect caused by NtC_i.

The NtC is then calculated from Effect_NtC with the fitted sigmoidal dose-response curve.

## 5. Measured value criterion

NtC_measured is the lowest tested chemical concentration that caused at least 10% effect
(survival ≤ 90%) in any of the biological replicates. When the criterion is applied and NtC_measured
is lower than the NtC from the fitted curve, the NtC is set to NtC_measured.

## 6. Replicates: same plate or independent experiments? (Updated algorithm only)

The width of the confidence band, and therefore the NtC, depends on how much independent information
the replicates carry.

- **Same plate (default):** the replicate columns are wells on one plate, measured in the same run.
  They share cells, dilution series and incubation, so their scatter does not show the variation
  between experiments (pseudo-replicates). The curve is fitted to the concentration means; degrees
  of freedom = *number of concentrations − 2*.
- **Independent experiments** (tick "Replicates are independent experiments"): each replicate column
  comes from a separate experiment (different day, cell passage, stock solution). The curve is
  fitted to all values; degrees of freedom = *number of values − 2*. Bootstrap confidence intervals
  are only available in this case.

If unsure, keep "same plate". Treating wells on one plate as independent overstates the precision of
the curve.

The choice can move the NtC in either direction. A wider band lowers NtC_lowerCI but raises
NtC_upperCI (the upper band stays above 100% up to higher concentrations), so "same plate" can give a
higher NtC than "independent".

The Original algorithm always fits the concentration means. Its degrees of freedom were number of
rows in the file − 3, which equals number of fitted concentrations − 2 only when one row is dropped as
control; the app uses number of fitted concentrations − 2 in both cases.

## 7. Original vs. Updated algorithm

- **Original (Stadnicka-Michalak et al.):** reproduces the published app, including its quirks.
  Where it would fail (e.g. the lower band never exceeds 90% in its fixed range), a message explains why.
- **Updated:** based on the corrected BendR implementation. Every difference, with its effect on the
  example data, is listed in the tab "Changes vs. original app".
