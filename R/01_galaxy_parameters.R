# =============================================================================
# Per-galaxy best-fit parameters and supplementary diagnostics
#
# Selects, for each galaxy, the best-fit eps_H (one- and two-parameter models)
# and central angle from the chi-square grids built during setup, writes the
# parameter table, and draws the supplementary 60-panel diagnostic figures.
# =============================================================================

# -----------------------------------------------------------------------------
# Supplementary figure: escape-Hubble ratio vs. acceleration excess per galaxy,
# with the fitted one- and two-parameter HMG curves.
# -----------------------------------------------------------------------------
{
  pdf(paste0("outputs/Suppl_fig_galaxies.pdf"), width = 7, height = 7)
  {
    par(mfrow = c(10, 6), mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))

    mcgaugh.Rg = rep(NA, length(mcgaugh$R))
    for (i in 1:length(gal_counts))
    {
      lg_gal = gal_name_row == names(gal_counts)[i]
      mcgaugh.Rg[lg_gal] = max(mcgaugh$R[lg_gal])

      mcgaugh.Vk = sqrt(mcgaugh$Vst^2 + mcgaugh$Vgas^2)
      mcgaugh.V  = (mcgaugh.Vk^2 * kms^2 + 2 * mcgaugh.Vk^2 * kms^2 * mcgaugh$R * kpc * c0 / (9 * T0))^0.25 / kms

      # Keplerian (baryonic) acceleration and derived quantities.
      a_newton = mcgaugh.Vk[lg_gal]^2 / (mcgaugh$R[lg_gal] * kpc) * kms^2
      aobs0 = mcgaugh$Vobs[lg_gal]^2 / (mcgaugh$R[lg_gal] * kpc) * kms^2
      vel_ratio_sq = mcgaugh$Vobs[lg_gal]^2 / mcgaugh.Vk[lg_gal]^2
      mass_bar = mcgaugh.Vk[lg_gal]^2 * (mcgaugh$R[lg_gal] * kpc) * kms^2
      gal_vevH = (sqrt(2) * mcgaugh.Vk[lg_gal] * kms * T0) / (mcgaugh$R[lg_gal] * kpc)
      gammaU = pi / 3

      g_ratio_btfr = mean(1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2), na.rm = T)
      g_ratio = 1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2 - 1)
      g_ratio = round(g_ratio, 2)
      g_sys1 = acos.((pi / 2) / g_ratio)
      g_sys1[g_sys1 < 1] = NA

      # Best-fit parameters for this galaxy from the stage-0 chi-square grids.
      gal_eps0.47 = seq_galeps[order(gal_chi2_0.47[i, ])[1]]
      gal_eps0.48 = seq_galeps[order(gal_chi2_0.48[i, ])[1]]
      gal_eps0.5  = seq_galeps[order(gal_chi2_0.5[i, ])[1]]
      gal_eps  = seq_galeps[order(apply(gal_chi22[i, , ], 1, min))[1]]
      gamma_cen  = seq_gammagc[order(apply(gal_chi22[i, , ], 2, min))[1]]
      gamempty = pi / 3
      gamma_cen_k[i]   = gamma_cen
      galaxy_eps_k[i]    = gal_eps
      galaxy_eps_0.47[i] = gal_eps0.47
      galaxy_eps_0.48[i] = gal_eps0.48
      galaxy_eps_0.5[i]  = gal_eps0.5

      # Two-parameter model prediction.
      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)

      # One-parameter model prediction (gamma_cen fixed at 0.48 pi).
      gamma_cen = 0.48 * pi
      gal_eps0 = gal_eps0.48
      quotient_vevh = abs(gal_vevH^2 - gal_eps0^2) / (gal_eps0^2 + gal_vevH^2)
      g_sys0_pred = asin((sin(gammaU)^2 + (sin(gamma_cen)^2 - sin(gammaU)^2) * quotient_vevh)^0.5)

      if (i != 58)
      {
        g1_pred = g_sys1_pred / cos(g_sys1_pred)
        g0_pred = g_sys0_pred / cos(g_sys0_pred)

        plot(xlab = "", ylab = "", gal_vevH / galaxy_eps_0.48[i], 1 / g_ratio, log = "x", ylim = c(0.01, 0.3), xlim = c(0.5, 1000), pch = 20, axes = FALSE)
        lines(gal_vevH / galaxy_eps_0.48[i], 1 / g0_pred, col = "indianred")
        lines(gal_vevH / galaxy_eps_0.48[i], 1 / g1_pred, col = "blue")
        text(100, 0.27, names(gal_counts)[i], cex = 0.75)
        box()

        if ((i - 1) %% 6 == 0 & i != 61)
          axis(side = 2, seq(0.0, 0.2, 0.1), hadj = 0.9, las = 2)
        if ((i %% 6 == 0 & i != 60) | i == 61)
          axis(side = 4, seq(0.0, 0.2, 0.1), hadj = 0.1, las = 2)
        if (i <= 6)
          axis(side = 3, c(1, 2, 10, 20, 50, 100, 200, 500), padj = 0.8)
        if (i >= 55)
          axis(side = 1, c(1, 2, 10, 20, 50, 100, 200, 500), padj = -0.8)
      }
    }

    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))
    plot.new()
    mtext(side = 1, expression(paste("Escape-Hubble ratio [", v[E] / (epsilon[H] * v[H]), "]")), line = 1.9)
    mtext(side = 3, expression(paste("Escape-Hubble ratio [", v[E] / (epsilon[H] * v[H]), "]")), line = 1.2)
    mtext(side = 4, expression(paste("[ ", italic(a)[Tot] - italic(a)[N], " ]/[ c/t ]")), line = 2.5)
    mtext(side = 2, expression(paste("[ ", italic(a)[Tot] - italic(a)[N], " ]/[ c/t ]")), line = 2.2)

    par(fig = c(0, 1, 0, 0.2), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 4, 3, 4))
    plot.new()
    legend("bottom", legend = c("Observations", "1-parameter model", "2-parameter model"), ncol = 3,
           col = c("black", "red", "blue"), lwd = c(0, 1, 1), lty = c(0, 1, 1), pch = 20, pt.cex = c(1, 0, 0), bty = "n")
  }
  dev.off()
}

# Global best-fit epsilon (one-parameter model) and parameter table.
galaxy_eps_0 = galaxy_eps_0.48
rapid = data.frame(Name = names(gal_counts)[-58], mod1_eps = galaxy_eps_0[-58],
                   mod2_eps = galaxy_eps_k[-58], mod2_gam = gamma_cen_k[-58] / pi)
write.table(rapid, row.names = FALSE, sep = "\t",
            file = paste0("outputs/Suppl_data_galaxy_parameters.txt"), quote = FALSE)


# -----------------------------------------------------------------------------
# Supplementary figure: rotation curves with the fitted HMG models.
# -----------------------------------------------------------------------------
{
  pdf(paste0("outputs/Suppl_fig_rotation_curves.pdf"), width = 7, height = 7)
  {
    par(mfrow = c(10, 6), mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))
    for (i in 1:length(gal_counts))
    {
      lg_gal = gal_name_row == names(gal_counts)[i]
      mcgaugh.Rg[lg_gal] = max(mcgaugh$R[lg_gal])

      mcgaugh.Vk = sqrt(mcgaugh$Vst^2 + mcgaugh$Vgas^2)
      mcgaugh.V  = (mcgaugh.Vk^2 * kms^2 + 2 * mcgaugh.Vk^2 * kms^2 * mcgaugh$R * kpc * c0 / (9 * T0))^0.25 / kms

      a_newton = mcgaugh.Vk[lg_gal]^2 / (mcgaugh$R[lg_gal] * kpc) * kms^2
      aobs0 = mcgaugh$Vobs[lg_gal]^2 / (mcgaugh$R[lg_gal] * kpc) * kms^2
      vel_ratio_sq = mcgaugh$Vobs[lg_gal]^2 / mcgaugh.Vk[lg_gal]^2
      mass_bar = mcgaugh.Vk[lg_gal]^2 * (mcgaugh$R[lg_gal] * kpc) * kms^2
      g_ratio = 1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2 - 1)
      gal_vevH = (sqrt(2) * mcgaugh.Vk[lg_gal] * kms * T0) / (mcgaugh$R[lg_gal] * kpc)
      gammaU = pi / 3

      gamempty = pi / 3
      gamma_cen = gamma_cen_k[i]
      gal_eps = galaxy_eps_k[i]

      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)

      gamma_cen = 0.48 * pi
      gal_eps0 = galaxy_eps_0.48[i]
      quotient_vevh = abs(gal_vevH^2 - gal_eps0^2) / (gal_eps0^2 + gal_vevH^2)
      g_sys0_pred = asin((sin(gammaU)^2 + (sin(gamma_cen)^2 - sin(gammaU)^2) * quotient_vevh)^0.5)

      if (i != 58)
      {
        g1_pred = g_sys1_pred / cos(g_sys1_pred)
        g0_pred = g_sys0_pred / cos(g_sys0_pred)
        v_obs   = mcgaugh.Vk[lg_gal] * (1 + 1 / a_newton * 2 * c0 / (g_ratio * T0))^(1 / 4)
        v0_pred = mcgaugh.Vk[lg_gal] * (1 + 1 / a_newton * 2 * c0 / (g0_pred * T0))^(1 / 4)
        v1_pred = mcgaugh.Vk[lg_gal] * (1 + 1 / a_newton * 2 * c0 / (g1_pred * T0))^(1 / 4)

        normalise = TRUE
        if (normalise)  denom = rep(1, length(lg_gal))
        if (!normalise) denom = mcgaugh$Vobs
        ylim = c(0.01, 1.1 * max(mcgaugh$Vobs) / max(denom))

        plot(xlab = "", ylab = "", mcgaugh$R[lg_gal], mcgaugh$Vobs[lg_gal] / denom[lg_gal], log = "x", ylim = ylim, xlim = c(0.5, 80), pch = 20, axes = FALSE)
        lines(mcgaugh$R[lg_gal], v0_pred / denom[lg_gal], col = "indianred")
        lines(mcgaugh$R[lg_gal], v1_pred / denom[lg_gal], col = "blue")
        text(40, 0.9 * ylim[2], names(gal_counts)[i], cex = 0.75)
        box()

        if (!normalise)
          if ((i - 1) %% 6 == 0 & i != 61)
            axis(side = 2, seq(0.0, 0.2, 0.1), hadj = 0.9, las = 2)
        if (!normalise)
          if ((i %% 6 == 0 & i != 60) | i == 61)
            axis(side = 4, seq(0.0, 0.2, 0.1), hadj = 0.12, las = 2)
        if ((i - 1) %% 6 == 0 & i != 61)
          axis(side = 2, seq(0, 250, 50), hadj = 0.9, las = 2)
        if ((i %% 6 == 0 & i != 60) | i == 61)
          axis(side = 4, seq(0, 250, 50), hadj = 0.12, las = 2)
        if (i <= 6)
          axis(side = 3, c(1, 2, 10, 20, 50, 100, 200, 500), padj = 0.8)
        if (i >= 55)
          axis(side = 1, c(1, 2, 10, 20, 50, 100, 200, 500), padj = -0.8)
      }
    }

    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))
    plot.new()
    mtext(side = 1, expression(paste(italic(R), " [kpc]")), line = 1.9)
    mtext(side = 3, expression(paste(italic(R), " [kpc]")), line = 1.2)
    mtext(side = 4, expression(paste(italic(v)[obs], " [km/s]")), line = 2.7)
    mtext(side = 2, expression(paste(italic(v)[obs], " [km/s]")), line = 2.4)

    par(fig = c(0, 1, 0, 0.2), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 4, 3, 4))
    plot.new()
    legend("bottom", legend = c("Observations", "1-parameter model", "2-parameter model"), ncol = 3,
           col = c("black", "red", "blue"), lwd = c(0, 1, 1), lty = c(0, 1, 1), pch = 20, pt.cex = c(1, 0, 0), bty = "n")
  }
  dev.off()
}


# -----------------------------------------------------------------------------
# Supplementary figure: predicted/observed acceleration ratio vs. g_N.
# -----------------------------------------------------------------------------
{
  pdf(paste0("outputs/Suppl_fig_interpolation.pdf"), width = 7, height = 7)
  {
    par(mfrow = c(10, 6), mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))
    for (i in 1:length(gal_counts))
    {
      lg_gal = gal_name_row == names(gal_counts)[i]
      mcgaugh.Rg[lg_gal] = max(mcgaugh$R[lg_gal])

      mcgaugh.Vk = sqrt(mcgaugh$Vst^2 + mcgaugh$Vgas^2)
      mcgaugh.V  = (mcgaugh.Vk^2 * kms^2 + 2 * mcgaugh.Vk^2 * kms^2 * mcgaugh$R * kpc * c0 / (9 * T0))^0.25 / kms

      a_newton = mcgaugh.Vk[lg_gal]^2 / (mcgaugh$R[lg_gal] * kpc) * kms^2
      aobs0 = mcgaugh$Vobs[lg_gal]^2 / (mcgaugh$R[lg_gal] * kpc) * kms^2
      vel_ratio_sq = mcgaugh$Vobs[lg_gal]^2 / mcgaugh.Vk[lg_gal]^2
      mass_bar = mcgaugh.Vk[lg_gal]^2 * (mcgaugh$R[lg_gal] * kpc) * kms^2
      g_ratio = 1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2 - 1)
      gal_vevH = (sqrt(2) * mcgaugh.Vk[lg_gal] * kms * T0) / (mcgaugh$R[lg_gal] * kpc)
      gammaU = pi / 3

      gamempty = pi / 3
      gamma_cen = gamma_cen_k[i]
      gal_eps = galaxy_eps_k[i]

      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)

      gamma_cen = 0.48 * pi
      gal_eps0 = galaxy_eps_0.48[i]
      quotient_vevh = abs(gal_vevH^2 - gal_eps0^2) / (gal_eps0^2 + gal_vevH^2)
      g_sys0_pred = asin((sin(gammaU)^2 + (sin(gamma_cen)^2 - sin(gammaU)^2) * quotient_vevh)^0.5)

      if (i != 58)
      {
        g1_pred = g_sys1_pred / cos(g_sys1_pred)
        g0_pred = g_sys0_pred / cos(g_sys0_pred)
        v_obs   = mcgaugh.Vk[lg_gal] * (1 + 1 / a_newton * 2 * c0 / (g_ratio * T0))^(1 / 4)
        v0_pred = mcgaugh.Vk[lg_gal] * (1 + 1 / a_newton * 2 * c0 / (g0_pred * T0))^(1 / 4)
        v1_pred = mcgaugh.Vk[lg_gal] * (1 + 1 / a_newton * 2 * c0 / (g1_pred * T0))^(1 / 4)

        ylim = c(-0.5, 0.5)
        plot(xlab = "", ylab = "", 0, 0, ylim = ylim, xlim = c(-12, -8), pch = 20, axes = FALSE)
        lines(log10(a_newton), log10(v0_pred^2 / mcgaugh$Vobs[lg_gal]^2), col = "indianred")
        lines(log10(a_newton), log10(v1_pred^2 / mcgaugh$Vobs[lg_gal]^2), col = "blue")
        text(-9, 0.9 * ylim[2], names(gal_counts)[i], cex = 0.75)
        box()

        if ((i - 1) %% 6 == 0 & i != 61)
          axis(side = 2, round(seq(-0.6, 0.1, 0.2), 1), hadj = 0.9, las = 2)
        if ((i %% 6 == 0 & i != 60) | i == 61)
          axis(side = 4, round(seq(-0.6, 0.1, 0.2), 1), hadj = 0.12, las = 2)
        if (i <= 6)
          axis(side = 3, seq(-10, -8, 1), padj = 0.8)
        if (i >= 55)
          axis(side = 1, seq(-10, -8, 1), padj = -0.8)
      }
    }

    par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))
    plot.new()
    mtext(side = 1, expression(paste(log[10](italic(g)[N]), " [", m / s^2, "]")), line = 2.2)
    mtext(side = 3, expression(paste(log[10](italic(g)[N]), " [", m / s^2, "]")), line = 1.12)
    mtext(side = 4, expression(paste(log[10](italic(g)[pred] / italic(g)[obs]))), line = 2.7)
    mtext(side = 2, expression(paste(log[10](italic(g)[pred] / italic(g)[obs]))), line = 2.5)

    par(fig = c(0, 1, 0, 0.2), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 4, 3, 4))
    plot.new()
    legend("bottom", legend = c("1-parameter model", "2-parameter model"), ncol = 3,
           col = c("red", "blue"), lwd = c(1, 1), lty = c(1, 1), pch = 20, pt.cex = c(1, 0, 0), bty = "n")
  }
  dev.off()
}
