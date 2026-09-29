#' Plot the fitted curve, confidence band, data and NtC
#'
#' @param result an `ntc_result`
#' @param units concentration units for the axis label
#' @return a ggplot object
#' @export
plot_ntc <- function(result, units = "units") {
  dat <- result$data
  f <- result$fit
  pts <- data.frame(conc = rep(dat$conc, ncol(dat$effects)), effect = c(dat$effects))
  pts <- pts[!is.na(pts$effect), ]
  xlab <- sprintf("Concentration (%s), log scale", units)

  lo <- min(dat$conc) / 3
  if (!is.null(f) && is.finite(f$NtC) && f$NtC > 0) lo <- min(lo, f$NtC / 3)
  hi <- max(dat$conc) * 2

  p <- ggplot2::ggplot() +
    ggplot2::geom_hline(yintercept = c(0, 90, 100), colour = "grey75", linewidth = 0.3, linetype = c(1, 3, 1)) +
    ggplot2::scale_x_log10(xlab, limits = c(lo, hi)) +
    ggplot2::scale_y_continuous("100% - effect (%)", breaks = seq(0, 100, 20)) +
    ggplot2::coord_cartesian(ylim = c(-15, 115)) +
    ggplot2::theme_bw(base_size = 13)

  if (!is.null(f)) {
    band <- f$band[10^f$band$log_conc >= lo & 10^f$band$log_conc <= hi, ]
    band$conc <- 10^band$log_conc
    p <- p +
      ggplot2::geom_ribbon(data = band, ggplot2::aes(x = conc, ymin = lower, ymax = upper),
                           fill = "#2c7fb8", alpha = 0.18) +
      ggplot2::geom_line(data = band, ggplot2::aes(x = conc, y = fit), colour = "#2c7fb8", linewidth = 0.9)
    crit <- data.frame(x = c(f$NtC_lower, f$NtC_upper), label = c("NtC_lowerCI", "NtC_upperCI"))
    crit <- crit[is.finite(crit$x) & crit$x >= lo & crit$x <= hi, ]
    if (nrow(crit)) {
      p <- p + ggplot2::geom_vline(data = crit, ggplot2::aes(xintercept = x), colour = "grey50", linetype = 3) +
        ggplot2::geom_text(data = crit, ggplot2::aes(x = x, y = 2, label = label), colour = "grey40",
                           size = 3.2, angle = 90, hjust = 0, vjust = -0.4)
    }
    if (is.finite(f$NtC) && f$NtC >= lo) {
      y_ntc <- hill(log10(f$NtC), f$logEC50, f$slope)
      p <- p + ggplot2::geom_vline(xintercept = f$NtC, colour = "#d95f02", linetype = 2) +
        ggplot2::annotate("point", x = f$NtC, y = y_ntc, shape = 4, size = 5, stroke = 1.3, colour = "#d95f02") +
        ggplot2::annotate("text", x = f$NtC, y = 110, label = sprintf("NtC = %s %s", fmt_num(f$NtC), units),
                          colour = "#d95f02", hjust = -0.05, size = 4)
    }
  }
  p <- p + ggplot2::geom_point(data = pts, ggplot2::aes(x = conc, y = effect), size = 2)
  title <- if (result$settings$algorithm == "original") "Original algorithm (Stadnicka-Michalak et al.)" else "Updated algorithm"
  p + ggplot2::ggtitle(title, if (result$status != "ok") "No NtC: see message" else NULL)
}
