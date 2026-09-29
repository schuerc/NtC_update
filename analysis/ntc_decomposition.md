# Why the original algorithm and BendR give different NtCs

Script: `analysis/ntc_decomposition.R` (run from the repository root). Data:
`analysis/data/synthetic_no_control.csv`, a **synthetic** dataset (6 concentrations x 3
replicates, no control row, about 80% survival at the lowest concentration) chosen because it shows
the same effects that were first seen on real data. Results of all 32 combinations:
`analysis/synthetic_no_control_factorial.csv`.

One parametrised implementation reproduces both the original algorithm (`StadnickaNtC.R`) and BendR
0.3.0 exactly. Five differences are switched on and off:

| Factor | Original | BendR 0.3.0 |
|---|---|---|
| `drop_first` | first data row always dropped as "control" | all rows used |
| `fit_to` | concentration means | all replicate values |
| `df_rule` | `nrow(file) - 3` | number of values - 2 |
| `upper_thr` | upper band >= 99.999999999 | upper band > 100 |
| `grid` | fixed range from concentration 1 (or 3 x log of the lowest) to 2 x log of the highest | +/-20% of the tested range, extended only if the lower criterion fails |

## Path from the original to BendR 0.3.0 (one change at a time)

| Step | NtC | Chosen by | df |
|---|---|---|---|
| Original | 1.005 | lower CI | 3 |
| + BendR range | 0.999 | lower CI | 3 |
| + upper threshold `> 100` | 0.999 | lower CI | 3 |
| + fit all replicates instead of means | 1.127 | lower CI | 3 |
| + df = number of values - 2 | 1.255 | upper CI | 13 |
| + keep the first row (= BendR 0.3.0) | **0.251** | lower CI / 10 | 16 |
| BendR logic with a wide enough range (= bendr >= 0.4.0) | **0.706** | upper CI | 16 |

## Findings

1. **Range artifact in BendR 0.3.0 (bug, fixed in bendr 0.4.0).** The upper band exceeds 100% only
   below the +/-20% range, and the range was only extended when the *lower* criterion failed. The
   upper criterion was then "not found", treated as effect 0; the ratio lower/upper became infinite
   and the narrow-CI correction (/10) fired. Here this lowers the NtC by a factor of 2.8
   (0.706 -> 0.251). The original has the same latent problem when the crossing lies below its
   fixed range.
2. **The original only works because it drops the first row.** With all rows kept, the lower band
   never exceeds 90% within the original's fixed range (it starts at concentration 1) and the script
   fails, in every combination that uses the original range. This matches the observation that the
   published app failed on a file for which the script returned a result.
3. **Degrees of freedom and fitted data move the NtC in either direction.** A narrower band (more df,
   replicates instead of means) lowers NtC_upperCI and raises NtC_lowerCI; which one decides depends
   on the data. Here the NtC rises from 1.005 to 1.255; on other data it fell by a factor of about 2.
   The replication unit (same plate vs independent experiments) therefore matters and is a setting in
   the app.
4. **The upper threshold never matters** (identical in all 32 combinations); **the range as such
   matters little** once both crossings lie inside it.
5. **Fitting means vs replicates does not give equal standard errors** (here SE(logEC50) 0.047 vs
   0.037 on the same rows): the two residual variances estimate different things (lack of fit of the
   means vs pooled replicate scatter).

## Consequences for the updated algorithm

- The concentration range is extended until *both* criteria are found; "not found in the range" is
  never treated as effect 0 (bendr >= 0.4.0).
- The replication unit is an explicit choice: same plate (concentration means, df = number of
  concentrations - 2, default) or independent experiments (all values, df = number of values - 2).
- Whether the first row is a control is an explicit choice (checkbox), ticked automatically when its
  concentration is 0.
