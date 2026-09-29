# Decompose the NtC difference between the original algorithm (StadnickaNtC.R) and BendR
# (bendr >= 0.4.0 logic, fixed range) on one dataset. One parametrised implementation;
# each factor switches one difference between the two.
#
# Run from the repository root (default: the synthetic dataset in analysis/data/):
#   Rscript analysis/ntc_decomposition.R [file.csv]
# The file must be semicolon-separated with a header; first column concentrations.

suppressMessages({
  library(minpack.lm)
})

args = commandArgs(trailingOnly = TRUE)
input = if (length(args)) args[1] else "analysis/data/synthetic_no_control.csv"
if (!exists("raw", inherits = FALSE) || !is.data.frame(raw)) raw = read.csv(input, sep = ";", check.names = FALSE)

hill = function(x, L, s) 100 / (1 + 10^((L - x) * s))
# gradient of hill() w.r.t. (logEC50, slope)
hill_grad = function(x, L, s) {
  u = (L - x) * s
  dfdu = -100 * log(10) * 10^u / (1 + 10^u)^2
  cbind(logEC50 = dfdu * s, slope = dfdu * (L - x))
}

ntc = function(drop_first = TRUE,     # original: data[2:nrow(data),] (row 1 = control)
               fit_to = "means",      # "means" (original) or "replicates" (BendR)
               df_rule = "file",      # "file": nrow(file)-3 (original); "fit": n_fit - 2 (BendR)
               upper_thr = "orig",    # "orig": >= 99.999999999; "bendr": > 100
               grid = "orig",         # "orig": 0 or 3*logconc1 .. 2*logconc_max, step 0.001; "bendr": +-20% + extend_ci;
                              # "wide": 3 log units below the lowest conc (reference, not in either implementation)
               details = FALSE) {
  d = if (drop_first) raw[-1, ] else raw
  conc = d[[1]]
  logconc = log10(conc)
  reps = as.matrix(d[, -1])

  fdat = if (fit_to == "means") {
    data.frame(x = logconc, y = rowMeans(reps))
  } else {
    data.frame(x = rep(logconc, ncol(reps)), y = c(reps))
  }
  fit = nlsLM(y ~ 100 / (1 + 10^((logEC50 - x) * slope)), data = fdat,
              start = c(logEC50 = log10(conc[ceiling(length(conc) / 2)]), slope = -2))
  L = unname(coef(fit)["logEC50"]); s = unname(coef(fit)["slope"])
  V = vcov(fit)
  dof = if (df_rule == "file") nrow(raw) - 3 else nrow(fdat) - 2
  tq = qt(0.975, dof)

  band = function(x) {
    J = hill_grad(x, L, s)
    se = sqrt(rowSums((J %*% V) * J))
    p = hill(x, L, s)
    data.frame(x = x, p = p, lo = p - tq * se, up = p + tq * se)
  }

  if (grid == "orig") {
    x_start = if (logconc[1] > 0) 0 else 3 * logconc[1]
    b = band(seq(x_start, logconc[length(logconc)] * 2, 0.001))
    extended = FALSE
  } else if (grid == "wide") {
    # 3 log units below the lowest tested concentration: both criteria can be found
    b = band(seq(logconc[1] - 3, logconc[length(logconc)] * 2, 0.001))
    extended = FALSE
  } else {
    r = range(logconc); ext = diff(r) * 0.2
    b = band(seq(r[1] - ext, r[2] + ext, length.out = 1000))
    extended = !any(b$lo > 90)
    if (extended) {
      lox = min(L - log10(100 / 99.9 - 1) / s, min(b$x))
      b = band(seq(lox, max(b$x), length.out = 1000))
    }
  }

  up_idx = if (upper_thr == "orig") which(b$up >= 99.999999999) else which(b$up > 100)
  lo_idx = which(b$lo > 90)
  if (length(lo_idx) == 0) {
    return(list(NtC = NA, NtC_model = NA, branch = "lower CI never > 90", NtC_lower = NA, NtC_upper = NA,
                eff_lower = NA, eff_upper = NA, EC50 = 10^L, slope = s, se_logEC50 = sqrt(V[1, 1]),
                se_slope = sqrt(V[2, 2]), dof = dof, t = tq, grid_from = 10^min(b$x), extended = extended))
  }

  if (length(up_idx) > 0) {
    i = max(up_idx); eff_up = 100 - b$p[i]
    NtC_up = if (eff_up > 0) 10^b$x[i] else 0
  } else { eff_up = 0; NtC_up = 0 }
  i = max(lo_idx); eff_lo = 100 - b$p[i]
  NtC_lo = if (eff_lo > 0) 10^b$x[i] else 0

  if (NtC_up > NtC_lo) {
    NtC = NtC_lo; branch = "lower"
  } else if (eff_lo / eff_up > 10) {
    eff = eff_lo / 10
    NtC = 10^(L - log10(100 / (100 - eff) - 1) / s); branch = "lower/10"
  } else {
    NtC = NtC_up; branch = "upper"
  }
  NtC_model = NtC
  below = conc < NtC
  if (any(below) && any(reps[below, , drop = FALSE] <= 90)) {
    NtC = min(conc[below][apply(reps[below, , drop = FALSE] <= 90, 1, any)]); branch = paste(branch, "+ measured")
  }

  out = list(NtC = NtC, NtC_model = NtC_model, branch = branch, NtC_lower = NtC_lo, NtC_upper = NtC_up,
             eff_lower = eff_lo, eff_upper = eff_up, EC50 = 10^L, slope = s,
             se_logEC50 = sqrt(V[1, 1]), se_slope = sqrt(V[2, 2]), dof = dof, t = tq,
             grid_from = 10^min(b$x), extended = extended)
  if (details) out$band = b
  out
}

fmt = function(r) {
  with(r, sprintf("NtC %6.3f  [%-18s]  NtC_lo %6.3f  NtC_up %6.3f  eff_lo %5.2f  eff_up %5.2f  EC50 %5.3f  slope %6.3f  se(logEC50) %.4f  se(slope) %.3f  df %2d  t %.2f  grid from %.3f%s",
                   NtC, branch, NtC_lower, NtC_upper, eff_lower, eff_upper, EC50, slope, se_logEC50, se_slope, dof, t, grid_from,
                   if (extended) " (extended)" else ""))
}

if (isTRUE(getOption("decomp.functions_only"))) stop("functions loaded", call. = FALSE)

orig  = list(drop_first = TRUE,  fit_to = "means",      df_rule = "file", upper_thr = "orig",  grid = "orig")
bendr = list(drop_first = FALSE, fit_to = "replicates", df_rule = "fit",  upper_thr = "bendr", grid = "bendr")

cat("== Reproduction ==\n")
cat("original :", fmt(do.call(ntc, orig)), "\n")
cat("BendR    :", fmt(do.call(ntc, bendr)), "\n\n")

cat("== One factor at a time, starting from ORIGINAL ==\n")
for (f in names(orig)) {
  a = orig; a[[f]] = bendr[[f]]
  cat(sprintf("%-10s -> %-10s: %s\n", f, bendr[[f]], fmt(do.call(ntc, a))))
}
cat("\n== One factor at a time, starting from BENDR ==\n")
for (f in names(bendr)) {
  a = bendr; a[[f]] = orig[[f]]
  cat(sprintf("%-10s -> %-10s: %s\n", f, orig[[f]], fmt(do.call(ntc, a))))
}

cat("\n== Full factorial (32 combinations) ==\n")
grid_all = expand.grid(lapply(names(orig), function(f) c(orig[[f]], bendr[[f]])), stringsAsFactors = FALSE)
names(grid_all) = names(orig)
grid_all$NtC = NA; grid_all$branch = NA
for (k in seq_len(nrow(grid_all))) {
  r = do.call(ntc, as.list(grid_all[k, names(orig)]))
  grid_all$NtC[k] = r$NtC; grid_all$branch[k] = r$branch
}
print(grid_all[order(grid_all$NtC), ], row.names = FALSE, digits = 3)
write.csv(grid_all, file.path(dirname(input), "..", paste0(tools::file_path_sans_ext(basename(input)), "_factorial.csv")), row.names = FALSE)

# main effects (mean NtC with factor at BendR level minus at original level)
cat("\n== Main effects on NtC (mean over the other factors) ==\n")
for (f in names(orig)) {
  m1 = mean(grid_all$NtC[grid_all[[f]] == bendr[[f]]], na.rm = TRUE); m0 = mean(grid_all$NtC[grid_all[[f]] == orig[[f]]], na.rm = TRUE)
  cat(sprintf("%-10s  original %.3f  BendR %.3f  diff %+.3f\n", f, m0, m1, m1 - m0))
}

cat("
== Path from original to BendR (one change at a time, cumulative) ==
")
path = list(c("grid", "bendr"), c("upper_thr", "bendr"), c("fit_to", "replicates"), c("df_rule", "fit"), c("drop_first", FALSE))
a = orig
cat(sprintf("%-28s %s
", "original", fmt(do.call(ntc, a))))
for (st in path) {
  a[[st[1]]] = if (st[1] == "drop_first") as.logical(st[2]) else st[2]
  cat(sprintf("%-28s %s
", paste("+", st[1], "=", st[2]), fmt(do.call(ntc, a))))
}

cat("
== Same configurations with a wide grid (upper crossing inside the grid) ==
")
for (cfg in list(original = orig, bendr = bendr)) {
  a = cfg; a$grid = "wide"
  cat(fmt(do.call(ntc, a)), "
")
}
