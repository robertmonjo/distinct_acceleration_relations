# =============================================================================
# Figure 5 + Table 1 - eps_H vs. neighbourhood density and (eps_H, gamma_cen)
# relation for galaxies and clusters.
#
# Also writes the cluster fit tables consumed by the manuscript (Table 1) and
# the supplementary variant.
# Outputs: outputs/Fig5.pdf, outputs/Table1_general_model.txt,
#          outputs/Table1_cluster_model.txt,
#          outputs/Suppl_table_general_model_variant.txt
# =============================================================================

# -----------------------------------------------------------------------------
# Confidence interval of a fitted parameter from its chi-square profile.
#
# Well-constrained clusters use the classical criterion: the grid points whose
# chi-square falls below the absolute threshold `thr` (the qchisq level for the
# number of data points). A cluster whose best fit already exceeds that
# threshold fails the chi-square test (p-value > 0.95); no grid point qualifies
# and the classical range() would return (Inf, -Inf). For those cases a mixed
# criterion is applied instead: the chi-square is rescaled so that the best fit
# has a reduced chi-square of one (standard error inflation for a poor fit) and
# the same absolute threshold is applied,  chi2 < thr * min(chi2) / dof.
conf_interval <- function(chi2, grid, thr, dof)
{
  below <- chi2 < thr
  if (any(below, na.rm = TRUE)) return(range(grid[below]))
  chi2_min <- min(chi2, na.rm = TRUE)
  range(grid[chi2 < thr * chi2_min / dof])
}

# -----------------------------------------------------------------------------
# Best-fit epsilon_H (and gamma_cen) per cluster from the chi-square grids, plus
# per-cluster redshift, mass, radius and escape-speed summaries.
# -----------------------------------------------------------------------------
for (iclu in 1:10)
{
  dof = rar.num[iclu] - 1
  eps_clusA[iclu, 1]   = seq_eps[order(clus_chi2A[iclu, ])[1]]
  eps_clusA[iclu, 2:3] = conf_interval(clus_chi2A[iclu, ], seq_eps, p95[iclu], dof)
  eps_clusB[iclu, 1]   = seq_eps[order(clus_chi2B[iclu, ])[1]]
  eps_clusB[iclu, 2:3] = conf_interval(clus_chi2B[iclu, ], seq_eps, p95[iclu], dof)
  eps_clusb[iclu, 1]   = seq_eps[order(clus_chi2b[iclu, ])[1]]
  eps_clusb[iclu, 2:3] = conf_interval(clus_chi2b[iclu, ], seq_eps, p95[iclu], dof)

  best_gal_C = order(apply(clus_chi2C[iclu, , ], 2, min))[1]
  best_eps_C = order(apply(clus_chi2C[iclu, , ], 1, min))[1]
  eps_clusC[iclu, 1]   = seq_eps[order(clus_chi2C[iclu, , best_gal_C])[1]]
  eps_clusC[iclu, 2:3] = conf_interval(clus_chi2C[iclu, , best_gal_C], seq_eps, p67[iclu], dof)
  gal_clusC[iclu, 1]   = seq_gal[order(clus_chi2C[iclu, best_eps_C, ])[1]]
  gal_clusC[iclu, 2:3] = conf_interval(clus_chi2C[iclu, best_eps_C, ], seq_gal, p67[iclu], dof)

  clus_z[iclu]       = max(rar$z[rar$Name == clusters[iclu]])
  clus_mass[iclu]    = max(rar.mass[rar$Name == clusters[iclu]])
  clus_massmin[iclu] = min(rar.mass[rar$Name == clusters[iclu]])
  clus_rad[iclu]     = max(rar.radius[rar$Name == clusters[iclu]])
  clus_radmin[iclu]  = min(rar.radius[rar$Name == clusters[iclu]])
  clus_ve[iclu, 1]   = mean(rar.escape_newton[rar$Name == clusters[iclu]])
  clus_ve[iclu, 2:3] = range(rar.escape_newton[rar$Name == clusters[iclu]])
  clus_vH[iclu, 1]   = mean(rar.escape_hubble[rar$Name == clusters[iclu]])
  clus_vH[iclu, 2:3] = range(rar.escape_hubble[rar$Name == clusters[iclu]])
  clus_vevH2[iclu, 1]   = mean((rar.escape_newton[rar$Name == clusters[iclu]] / rar.escape_hubble[rar$Name == clusters[iclu]])^2)
  clus_vevH2[iclu, 2:3] = range((rar.escape_newton[rar$Name == clusters[iclu]] / rar.escape_hubble[rar$Name == clusters[iclu]])^2)
}

# -----------------------------------------------------------------------------
# Table 1 - one-parameter (general) and two-parameter (cluster) fits per cluster.
# -----------------------------------------------------------------------------
table_clusB = cbind(eps_clusB, round(pchisq(apply(clus_chi2B[, ], 1, min), rar.num - 1), 2))
colnames(table_clusB) = c("eps_med", "eps_low", "eps_upp", "p_value")
write.table(table_clusB, "outputs/Table1_general_model.txt", quote = FALSE)

table_clusb = cbind(eps_clusb, round(pchisq(apply(clus_chi2b[, ], 1, min), rar.num - 1), 2))
colnames(table_clusb) = c("eps_med", "eps_low", "eps_upp", "p_value")
write.table(table_clusb, "outputs/Suppl_table_general_model_variant.txt", quote = FALSE)

table_clusC = cbind(eps_clusC, round(gal_clusC, 3), round(pchisq(apply(clus_chi2C[, , ], 1, min), rar.num - 1), 2))
colnames(table_clusC) = c("eps_med", "eps_low", "eps_upp", "ggal_med", "ggal_low", "ggal_upp", "p_value")
write.table(table_clusC, "outputs/Table1_cluster_model.txt", quote = FALSE)


# -----------------------------------------------------------------------------
# Figure 5.
#   Left : eps_H vs. square root of the relative neighbourhood density.
#   Right: cos(gamma_cen) vs. log(1/eps_H).
# -----------------------------------------------------------------------------
{
  pdf(paste0("outputs/Fig5.pdf"), width = 7, height = 4.5)
  {
    par(fig = c(0, 0.48, 0.1, 1), oma = c(5.9, 3, 0.5, 3), mar = c(0, 0, 0, 0))
    plot.new()
    rect(-10, -10, 10, 10, col = rgb(1, 0.9, 0.9, 0.1))
    par(fig = c(0, 0.48, 0.1, 1), oma = c(5.9, 3, 0.5, 3), mar = c(0, 0, 0, 0), new = TRUE)

    # ---- One-parameter model: galaxies and clusters ----
    mcgaugh.VO = aggregate(mcgaugh.V, by = list(gal_name_row), mean, na.rm = T)
    galaxy_eps_k_err  = apply(cbind(galaxy_eps_0.47, galaxy_eps_0.48, galaxy_eps_k), 1, sd) / sqrt(2)
    gamma_cen_k_err = apply(cbind(0.47 * pi, 0.48 * pi, gamma_cen_k), 1, sd) / 2

    xlab = expression(paste(sqrt(rho[typ] / rho[vac])))
    ylab = expression(paste(epsilon[H]))

    x0 = sqrt(mcgaugh.mass / (4 / 3 * pi * (mcgaugh.R$x * kpc * 4)^3) / (3 / (8 * pi * GN * T0^2)))
    plot(x0, galaxy_eps_0.48 - 5 / 6, lwd = 1.5, xlab = xlab, ylab = ylab, xlim = c(3, 180), ylim = c(3, 180), log = "xy", pch = 0, axes = FALSE,
         col = rainbow(length(gal_counts))[order(mcgaugh.VO$x, decreasing = TRUE)[-55]])
    segments(x0, galaxy_eps_0.48 - galaxy_eps_k_err, x0, galaxy_eps_0.48 + galaxy_eps_k_err,
             col = rainbow(length(gal_counts))[order(mcgaugh.VO$x, decreasing = TRUE)[-55]])

    x0 = sqrt(clus_massmin / (4 / 3 * pi * ((clus_radmin)^3)) / (3 / (8 * pi * GN * (1 + clus_z)^(-1) * T0^2)))
    points(x0, eps_clusB[, 1], pch = 20, cex = 1.5, col = rainbow(length(clusters)))
    segments(x0, eps_clusB[, 2], x0, eps_clusB[, 3], col = rainbow(length(clusters)))

    box()
    abline(0, 1)
    axis(side = 1, c(5, 10, 20, 50, 100), padj = -0.8)
    axis(side = 2, c(5, 10, 20, 50, 100), padj = 0.8)
    abline(v = c(5, 10, 20, 50, 100), h = c(5, 10, 20, 50, 100), lty = 2, col = "gray95")
    mtext(side = 1, xlab, line = 1.99)
    mtext(side = 2, ylab, line = 1.6)

    # ---- Two-parameter model: galaxies and clusters ----
    par(fig = c(0.52, 1, 0.1, 1), oma = c(5.9, 3, 0.5, 3), mar = c(0, 0, 0, 0), new = TRUE)
    plot.new()
    rect(-10, -10, 10, 10, col = rgb(1, 0.9, 0.9, 0.1))
    par(fig = c(0.52, 1, 0.1, 1), oma = c(5.9, 3, 0.5, 3), mar = c(0, 0, 0, 0), new = TRUE)

    xlab = expression(paste(log(1 / epsilon[H])))
    ylab = expression(paste(cos(gamma[cen])))

    x0 = -log(galaxy_eps_k)[gal_counts > 1] - abs(seq(0, 0.05, length.out = sum(gal_counts > 1)))
    xer1 = -log(galaxy_eps_k * 0.3)[gal_counts > 1]
    xer1[is.na(xer1)] = x0
    plot(x0, cos(gamma_cen_k)[gal_counts > 1], pch = 0, xlab = xlab, ylab = ylab, axes = FALSE, ylim = c(0, 0.16),
         col = rainbow(length(gal_counts))[order(mcgaugh.VO$x, decreasing = TRUE)[-55]], lwd = 1.5)
    segments(xer1, cos(gamma_cen_k)[gal_counts > 1], -log(galaxy_eps_k + galaxy_eps_k_err / 2)[gal_counts > 1], cos(gamma_cen_k)[gal_counts > 1],
             col = rainbow(length(gal_counts))[order(mcgaugh.VO$x, decreasing = TRUE)[-55]])
    segments(x0, cos(gamma_cen_k - gamma_cen_k_err)[gal_counts > 1], x0, cos(gamma_cen_k + gamma_cen_k_err)[gal_counts > 1],
             col = rainbow(length(gal_counts))[order(mcgaugh.VO$x, decreasing = TRUE)[-55]])
    box()
    abline(a = cos(0.46 * pi), b = 0.02)

    points(-log(eps_clusC[-3, 1]), cos(gal_clusC[-3, 1] * pi), pch = 20, cex = 1.5, col = rainbow(length(clusters)))
    segments(-log(eps_clusC[-3, 1]), cos(gal_clusC[-3, 2] * pi), -log(eps_clusC[-3, 1]), cos(gal_clusC[-3, 3] * pi), col = rainbow(length(clusters)))
    segments(-log(eps_clusC[-3, 2]), cos(gal_clusC[-3, 1] * pi), -log(eps_clusC[-3, 3]), cos(gal_clusC[-3, 1] * pi), col = rainbow(length(clusters)))

    axis(side = 1, seq(-5, 0, 1), padj = -0.8)
    axis(side = 4, seq(0, 0.14, 0.02), padj = -0.8)
    mtext(side = 1, xlab, line = 1.8)
    mtext(side = 4, ylab, line = 1.6)
    abline(v = seq(-5, 0, 1), h = seq(0, 0.14, 0.02), lty = 2, col = "gray95")

    # ---- Legends ----
    par(fig = c(0, 1, 0, 0.181), mar = c(0, 0, 0, 0), oma = c(0.01, 0.1, 1, 0.1), new = TRUE)
    plot.new()
    legend("bottomleft", legend = clusters, pch = 20, col = clus_col, ncol = 2, pt.cex = cex0 * 0.61, cex = cex0 * 0.45, bty = "n")
    legend("bottomright", legend = names(gal_names_up)[order(mcgaugh.VO$x, decreasing = TRUE)[-55]], pch = 0, col = rainbow(length(gal_counts))[-55], ncol = 10, pt.cex = cex0 * 0.48, cex = cex0 * 0.365, bty = "n")

    par(fig = c(0, 1, 0, 0.225), mar = c(0, 0, 0, 0), oma = c(0, 0.1, 0, 0.45), new = T)
    plot.new()
    legend("topleft", cex = 0.801, legend = expression(bold("Clusters:")), pch = 20, col = "white", bty = "n")
    legend("topright", cex = 0.801, legend = expression(bold("Galaxies:")), pch = 20, col = "white", bty = "n", pt.cex = 0.001)
  }
  dev.off()
}
