# =============================================================================
# Figure 1 - radial acceleration relation of galaxy clusters under HMG.
#
# Left panel : modelled vs. observed a_tot/a_N with the two-parameter fit.
# Right panel: acceleration excess vs. Newton-Hubble speed ratio, with the
#              MOND-like, general and cluster HMG model bands.
# Output: outputs/Fig1.pdf
# =============================================================================

g_galaxy     = 0.5 * pi
gblack_hole  = 0.5 * pi
gempty_space = pi / 3

# Global one-parameter model curve at the median epsilon_H.
eps_H0 = median(eps_H, na.rm = T) / 2
vE2vH2 = (rar.escape_newton / rar.escape_hubble)^2
dens_div_dens = (rar.escape_newton / rar.escape_hubble)^2 / (eps_H0^2 + (rar.escape_newton / rar.escape_hubble)^2)
gamma_m1 = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(gblack_hole)^2 - sin(gempty_space)^2)))
gamma_01 = gamma_m1 / cos(gamma_m1)

# Sampling grids for the model curves.
vE2vH2_sort  = seq(min(vE2vH2), max(vE2vH2), length.out = 300)
vEvH_sort    = sqrt(vE2vH2_sort)
vE2vH2_sort2 = 10^seq(0, 8, length.out = 500)
vEvH_sort2   = sqrt(vE2vH2_sort2)

# Cluster-model curve with median / upper / lower epsilon_H (from stage Fig4/Fig5).
dens_div_dens = vE2vH2_sort / (mean(eps_clusC[-3, 1])^2 + vE2vH2_sort)
gamma_Mmed = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(gblack_hole)^2 - sin(gempty_space)^2)))
gamma_0med = gamma_Mmed / cos(gamma_Mmed)
dens_div_dens = vE2vH2_sort / (mean(eps_clusC[-3, 3])^2 + vE2vH2_sort)
gamma_Mup = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(gblack_hole)^2 - sin(gempty_space)^2)))
gamma_0up = gamma_Mup / cos(gamma_Mup)
dens_div_dens = vE2vH2_sort / (mean(eps_clusC[-3, 2])^2 + vE2vH2_sort)
gamma_Mdow = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(gblack_hole)^2 - sin(gempty_space)^2)))
gamma_0dow = gamma_Mdow / cos(gamma_Mdow)

# General-model band: envelope over an epsilon_H range.
eps_min = mean(eps_clusB[-3, 2])
eps_max = mean(eps_clusB[-3, 3])
eps_seq = seq(eps_min, eps_max, length.out = 10)
gamma_0mdd_seq  = array(NA, dim = c(length(vEvH_sort), length(eps_seq)))
gamma_0mdd2_seq = array(NA, dim = c(length(vEvH_sort2), length(eps_seq)))
for (iseq in 1:length(eps_seq))
{
  dens_div_dens = abs(vEvH_sort^2 - eps_seq[iseq]^2) / (eps_seq[iseq]^2 + vEvH_sort^2)
  gamma_M = asin(sqrt(sin(gempty_space)^2 - dens_div_dens * (sin(gempty_space)^2 - sin(g_galaxy)^2)))
  gamma_0mdd_seq[, iseq] = gamma_M / cos(gamma_M)

  dens_div_dens2 = abs(vEvH_sort2^2 - eps_seq[iseq]^2) / (eps_seq[iseq]^2 + vEvH_sort2^2)
  gamma_M2 = asin(sqrt(sin(gempty_space)^2 - dens_div_dens2 * (sin(gempty_space)^2 - sin(g_galaxy)^2)))
  gamma_0mdd2_seq[, iseq] = gamma_M2 / cos(gamma_M2)
}
gamma_0mdd  = apply(gamma_0mdd_seq, 1, quantile, 0.5)
gamma_0mup  = apply(gamma_0mdd_seq, 1, quantile, 0.9)
gamma_0mdo  = apply(gamma_0mdd_seq, 1, quantile, 0.1)
gamma_0mdd2 = apply(gamma_0mdd2_seq, 1, quantile, 0.5)
gamma_0mup2 = apply(gamma_0mdd2_seq, 1, quantile, 0.9)
gamma_0mdo2 = apply(gamma_0mdd2_seq, 1, quantile, 0.1)


{
  eps0 = mean(eps_clusC[, 1])
  eps0 = quantile(c(eps_clusB[, 1]), 0.5, na.rm = TRUE)

  pdf(paste0("outputs/Fig1.pdf"), width = 10, height = 5.05)
  {
    cex0 = 1.5
    gamma_00 = gamma_01
    gempty_space = pi / 3
    xmin = 2

    # ---- Left panel: modelled vs. observed a_tot/a_N ----
    par(fig = c(0, 0.5, 0.0941, 1), mar = c(3, 3, 0.5, 0.5))
    plot(((1 + rar$z) * c0) / (T0 * gamma_01 * rar.acc_newton), quot_accel, col = "white",
         xlab = "", ylab = "", axes = FALSE, xlim = c(xmin, 30), ylim = c(xmin, 30), log = "xy")
    clusters = as.character(unique(rar$Name))
    clus_col = rainbow(length(clusters))
    mtext(side = 1, expression(paste("", italic(a)[Tot] / italic(a)[N], "  modelled")), line = 1.8)
    mtext(side = 2, expression(paste("", italic(a)[Tot] / italic(a)[N], "  observed")), line = 1.6)
    axis(side = 1, c(2, 3, 5, 7, 10, 15, 20, 30), padj = -0.8)
    axis(side = 2, c(2, 3, 5, 7, 10, 15, 20, 30), padj = 0.8)
    abline(v = c(2, 3, 5, 7, 10, 15, 20, 30), h = c(2, 3, 5, 7, 10, 15, 20, 30), lty = 2, col = "gray95")
    abline(0, 1, lwd = 2, col = rgb(0.8, 0.4, 0.4, 1), lty = 2)
    box()

    for (iclu in 1:length(clusters))
    {
      eps_H0    = eps_clusC[iclu, 1]
      eps_H0_do = eps_clusC[iclu, 2]
      eps_H0_up = eps_clusC[iclu, 3]
      g_galaxy  = gal_clusC[iclu, 1] * pi

      dens_div_dens = (rar.escape_newton / rar.escape_hubble)^2 / (eps_H0^2 + (rar.escape_newton / rar.escape_hubble)^2)
      gamma_m0 = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(g_galaxy)^2 - sin(gempty_space)^2)))
      gamma_00 = gamma_m0 / cos(gamma_m0)

      eps_H0 = eps_H0_do
      dens_div_dens = (rar.escape_newton / rar.escape_hubble)^2 / (eps_H0^2 + (rar.escape_newton / rar.escape_hubble)^2)
      gamma_m2 = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(g_galaxy)^2 - sin(gempty_space)^2)))
      gamma_do = gamma_m2 / cos(gamma_m2)

      eps_H0 = eps_H0_up
      dens_div_dens = (rar.escape_newton / rar.escape_hubble)^2 / (eps_H0^2 + (rar.escape_newton / rar.escape_hubble)^2)
      gamma_m2 = asin(sqrt(sin(g_galaxy)^2 - dens_div_dens * (sin(g_galaxy)^2 - sin(gempty_space)^2)))
      gamma_up = gamma_m2 / cos(gamma_m2)

      lgclus = rar$Name == clusters[iclu]
      points((((1 + rar$z) * c0) / (T0 * gamma_00 * rar.acc_newton))[lgclus], quot_accel[lgclus], pch = 20, cex = cex0, col = clus_col[iclu])
      segments((((1 + rar$z) * c0) / (T0 * gamma_00 * rar.acc_newton))[lgclus], quot_accel_do[lgclus], (((1 + rar$z) * c0) / (T0 * gamma_00 * rar.acc_newton))[lgclus], quot_accel_up[lgclus], col = clus_col[iclu])
      segments((((1 + rar$z) * c0) / (T0 * gamma_do * rar.acc_newton))[lgclus], quot_accel[lgclus], (((1 + rar$z) * c0) / (T0 * gamma_up * rar.acc_newton))[lgclus], quot_accel[lgclus], col = clus_col[iclu])
    }

    eps00 = eps0 / sqrt(2)

    # ---- Right panel: acceleration excess vs. Newton-Hubble speed ratio ----
    par(fig = c(0.5, 1, 0.0941, 1), new = T, mar = c(3, 3, 0.5, 0.5))
    plot(vEvH / eps0 * eps00, diff_accel, col = "white", xlab = "", ylab = "", axes = FALSE, log = "xy", ylim = c(0.004, 1.5))
    box()
    clusters = as.character(unique(rar$Name))
    clus_col = rainbow(length(clusters))

    ones = rep(1, length(vE2vH2_sort))
    polygon_(sqrt(vE2vH2_sort) / eps0 * eps00, ones * cos(0.466 * pi) / (0.466 * pi), ones * cos(0.477 * pi) / (0.477 * pi), ones * cos(0.456 * pi) / (0.456 * pi), add = TRUE, col = rgb(0.12, 0.65, 0.99, 0.08))
    polygon_(sqrt(vE2vH2_sort) / eps0 * eps00, 1 / gamma_0mdd, 1 / gamma_0mdo, 1 / gamma_0mup, add = TRUE, col = rgb(0.75, 0.75, 0.75, 0.3))
    polygon_(sqrt(vE2vH2_sort) / eps0 * eps00, 1 / gamma_0med, 1 / gamma_0dow, 1 / gamma_0up, add = TRUE, col = rgb(0.8, 0.4, 0.4, 0.3))
    mtext(side = 1, expression(v[N] / v[H]), line = 1.8)
    mtext(side = 2, expression(paste("[ ", italic(a)[Tot] - italic(a)[N], " ]/[ c/t ]")), line = 1.6)
    axis(side = 1, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200),
         paste(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)), padj = -0.8)
    abline(v = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100), h = c(1, 2, 5, 10, 20, 50, 100) / 100, lty = 2, col = "gray95")
    axis(side = 2, c(1, 5) / 100, padj = 0.8)
    axis(side = 2, 2 / 100, padj = 0.8)
    axis(side = 2, c(0.1, 0.2, 0.5), padj = 0.8)
    axis(side = 2, 1, padj = 0.8)

    for (iclu in 1:length(clusters))
    {
      lgclus = rar$Name == clusters[iclu]
      eps_H0    = eps_clusB[iclu, 1]
      eps_H0_do = eps_clusB[iclu, 2]
      eps_H0_up = eps_clusB[iclu, 3]
      eps00 = eps_H0 / sqrt(2)
      points(vEvH[lgclus] / eps_H0 * eps00, diff_accel[lgclus], pch = 20, cex = cex0, col = clus_col[iclu])
      segments(vEvH[lgclus] / eps_H0 * eps00, diff_accel_do[lgclus], vEvH[lgclus] / eps_H0 * eps00, diff_accel_up[lgclus], col = clus_col[iclu])
    }

    legend("bottomright", ncol = 2, legend = c("MOND-like approach", "General model", "Cluster approach"), cex = 0.90,
           col = c(rgb(0.12, 0.65, 0.99, 0.08), rgb(0.75, 0.75, 0.75, 0.3), rgb(0.8, 0.4, 0.4, 0.3)), lty = c(1, 1, 1), lwd = c(6, 6.5, 6), bty = "n")
    legend("bottomright", ncol = 2, legend = c("MOND-like approach", "General model", "Cluster approach"), cex = 0.90,
           col = c(rgb(0.12, 0.65, 0.99, 0.2), "gray70", rgb(0.8, 0.4, 0.4, 1)), lty = c(1, 1, 2), lwd = c(1.8, 1.8, 1.8), bty = "n")

    abline(h = c(0.5, 1), col = "gray65", lwd = 1, lty = 2)
    text(2.15, 0.55, "Empty-space limit", cex = 0.6, col = "gray75")
    text(2.05, 1.10, "Causality limit", cex = 0.6, col = "gray75")

    eps00 = eps0 / sqrt(2)
    lines(sqrt(vE2vH2_sort) / eps0 * eps00, ones * cos(0.466 * pi) / (0.466 * pi), lwd = 2, col = rgb(0.12, 0.65, 0.99, 0.2), lty = 1)
    lines(sqrt(vE2vH2_sort) / eps0 * eps00, 1 / gamma_0mdd, lwd = 2, col = "gray70", lty = 1)
    lines(sqrt(vE2vH2_sort) / eps0 * eps00, 1 / gamma_0med, lwd = 2, col = rgb(0.8, 0.4, 0.4, 1), lty = 2)

    # ---- Cluster colour legend across the bottom ----
    par(fig = c(0, 1, 0, 0.0991), mar = c(0, 0, 0, 0), new = T, oma = c(0, 1, 0, 1))
    plot.new()
    legend("right", legend = clusters, pch = 20, col = clus_col, ncol = 10, pt.cex = cex0, bty = "n")
    legend("left", legend = expression(bold("Clusters:")), pch = 20, col = "white", bty = "n")
  }
  dev.off()
}
