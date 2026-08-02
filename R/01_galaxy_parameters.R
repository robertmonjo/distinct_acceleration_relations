# =============================================================================
# Stage 1 - per-galaxy best-fit parameters and supplementary diagnostics
#
# Selects, for each galaxy, the best-fit eps_H (one- and two-parameter models)
# and central angle from the chi-square grids built in stage 0, writes the
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

    smcGa07.Rg = rep(NA, length(smcGa07$R))
    for (i in 1:length(namesg))
    {
      lg_gal = namess == names(namesg)[i]
      smcGa07.Rg[lg_gal] = max(smcGa07$R[lg_gal])

      smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
      smcGa07.V  = (smcGa07.Vk^2 * kms^2 + 2 * smcGa07.Vk^2 * kms^2 * smcGa07$R * kpc * c0 / (9 * T0))^0.25 / kms

      # Keplerian (baryonic) acceleration and derived quantities.
      ak000 = smcGa07.Vk[lg_gal]^2 / (smcGa07$R[lg_gal] * kpc) * kms^2
      aobs0 = smcGa07$Vobs[lg_gal]^2 / (smcGa07$R[lg_gal] * kpc) * kms^2
      di000 = smcGa07$Vobs[lg_gal]^2 / smcGa07.Vk[lg_gal]^2
      mk000 = smcGa07.Vk[lg_gal]^2 * (smcGa07$R[lg_gal] * kpc) * kms^2
      gal_vevH = (sqrt(2) * smcGa07.Vk[lg_gal] * kms * T0) / (smcGa07$R[lg_gal] * kpc)
      gammaU = pi / 3

      gg000_btfr = mean(1 / ak000 * (2 * c0 / T0) / (di000^2), na.rm = T)
      gg000 = 1 / ak000 * (2 * c0 / T0) / (di000^2 - 1)
      gg000 = round(gg000, 2)
      g_sys1 = acos.((pi / 2) / gg000)
      g_sys1[g_sys1 < 1] = NA

      # Best-fit parameters for this galaxy from the stage-0 chi-square grids.
      gal_eps0.47 = seq_galeps[order(gal_xi2_0.47[i, ])[1]]
      gal_eps0.48 = seq_galeps[order(gal_xi2_0.48[i, ])[1]]
      gal_eps0.5  = seq_galeps[order(gal_xi2_0.5[i, ])[1]]
      gal_eps  = seq_galeps[order(apply(gal_xi22[i, , ], 1, min))[1]]
      gammagc  = seq_gammagc[order(apply(gal_xi22[i, , ], 2, min))[1]]
      gamempty = pi / 3
      gammagc_k[i]   = gammagc
      galeps_k[i]    = gal_eps
      galeps_0.47[i] = gal_eps0.47
      galeps_0.48[i] = gal_eps0.48
      galeps_0.5[i]  = gal_eps0.5

      # Two-parameter model prediction.
      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)

      # One-parameter model prediction (gamma_cen fixed at 0.48 pi).
      gammagc = 0.48 * pi
      gal_eps0 = gal_eps0.48
      quotient_vevh = abs(gal_vevH^2 - gal_eps0^2) / (gal_eps0^2 + gal_vevH^2)
      g_sys0_pred = asin((sin(gammaU)^2 + (sin(gammagc)^2 - sin(gammaU)^2) * quotient_vevh)^0.5)

      if (i != 58)
      {
        g1_pred = g_sys1_pred / cos(g_sys1_pred)
        g0_pred = g_sys0_pred / cos(g_sys0_pred)

        plot(gal_vevH / galeps_0.48[i], 1 / gg000, log = "x", ylim = c(0.01, 0.3), xlim = c(0.5, 1000), pch = 20, axes = FALSE)
        lines(gal_vevH / galeps_0.48[i], 1 / g0_pred, col = "indianred")
        lines(gal_vevH / galeps_0.48[i], 1 / g1_pred, col = "blue")
        text(100, 0.27, names(namesg)[i], cex = 0.75)
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
galeps_0 = galeps_0.48
rapid = data.frame(Name = names(namesg)[-58], mod1_eps = galeps_0[-58],
                   mod2_eps = galeps_k[-58], mod2_gam = gammagc_k[-58] / pi)
write.table(rapid, row.names = FALSE, sep = "\t",
            file = paste0("outputs/Suppl_data_galaxy_parameters.txt"), quote = FALSE)


# -----------------------------------------------------------------------------
# Supplementary figure: rotation curves with the fitted HMG models.
# -----------------------------------------------------------------------------
{
  pdf(paste0("outputs/Suppl_fig_rotation_curves.pdf"), width = 7, height = 7)
  {
    par(mfrow = c(10, 6), mar = c(0, 0, 0, 0), oma = c(5, 4, 3, 4))
    for (i in 1:length(namesg))
    {
      lg_gal = namess == names(namesg)[i]
      smcGa07.Rg[lg_gal] = max(smcGa07$R[lg_gal])

      smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
      smcGa07.V  = (smcGa07.Vk^2 * kms^2 + 2 * smcGa07.Vk^2 * kms^2 * smcGa07$R * kpc * c0 / (9 * T0))^0.25 / kms

      ak000 = smcGa07.Vk[lg_gal]^2 / (smcGa07$R[lg_gal] * kpc) * kms^2
      aobs0 = smcGa07$Vobs[lg_gal]^2 / (smcGa07$R[lg_gal] * kpc) * kms^2
      di000 = smcGa07$Vobs[lg_gal]^2 / smcGa07.Vk[lg_gal]^2
      mk000 = smcGa07.Vk[lg_gal]^2 * (smcGa07$R[lg_gal] * kpc) * kms^2
      gg000 = 1 / ak000 * (2 * c0 / T0) / (di000^2 - 1)
      gal_vevH = (sqrt(2) * smcGa07.Vk[lg_gal] * kms * T0) / (smcGa07$R[lg_gal] * kpc)
      gammaU = pi / 3

      gamempty = pi / 3
      gammagc = gammagc_k[i]
      gal_eps = galeps_k[i]

      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)

      gammagc = 0.48 * pi
      gal_eps0 = galeps_0.48[i]
      quotient_vevh = abs(gal_vevH^2 - gal_eps0^2) / (gal_eps0^2 + gal_vevH^2)
      g_sys0_pred = asin((sin(gammaU)^2 + (sin(gammagc)^2 - sin(gammaU)^2) * quotient_vevh)^0.5)

      if (i != 58)
      {
        g1_pred = g_sys1_pred / cos(g_sys1_pred)
        g0_pred = g_sys0_pred / cos(g_sys0_pred)
        v_obs   = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (gg000 * T0))^(1 / 4)
        v0_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g0_pred * T0))^(1 / 4)
        v1_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g1_pred * T0))^(1 / 4)

        entre_uno = TRUE
        if (entre_uno)  divido = rep(1, length(lg_gal))
        if (!entre_uno) divido = smcGa07$Vobs
        ylim = c(0.01, 1.1 * max(smcGa07$Vobs) / max(divido))

        plot(smcGa07$R[lg_gal], smcGa07$Vobs[lg_gal] / divido[lg_gal], log = "x", ylim = ylim, xlim = c(0.5, 80), pch = 20, axes = FALSE)
        lines(smcGa07$R[lg_gal], v0_pred / divido[lg_gal], col = "indianred")
        lines(smcGa07$R[lg_gal], v1_pred / divido[lg_gal], col = "blue")
        text(40, 0.9 * ylim[2], names(namesg)[i], cex = 0.75)
        box()

        if (!entre_uno)
          if ((i - 1) %% 6 == 0 & i != 61)
            axis(side = 2, seq(0.0, 0.2, 0.1), hadj = 0.9, las = 2)
        if (!entre_uno)
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
    for (i in 1:length(namesg))
    {
      lg_gal = namess == names(namesg)[i]
      smcGa07.Rg[lg_gal] = max(smcGa07$R[lg_gal])

      smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
      smcGa07.V  = (smcGa07.Vk^2 * kms^2 + 2 * smcGa07.Vk^2 * kms^2 * smcGa07$R * kpc * c0 / (9 * T0))^0.25 / kms

      ak000 = smcGa07.Vk[lg_gal]^2 / (smcGa07$R[lg_gal] * kpc) * kms^2
      aobs0 = smcGa07$Vobs[lg_gal]^2 / (smcGa07$R[lg_gal] * kpc) * kms^2
      di000 = smcGa07$Vobs[lg_gal]^2 / smcGa07.Vk[lg_gal]^2
      mk000 = smcGa07.Vk[lg_gal]^2 * (smcGa07$R[lg_gal] * kpc) * kms^2
      gg000 = 1 / ak000 * (2 * c0 / T0) / (di000^2 - 1)
      gal_vevH = (sqrt(2) * smcGa07.Vk[lg_gal] * kms * T0) / (smcGa07$R[lg_gal] * kpc)
      gammaU = pi / 3

      gamempty = pi / 3
      gammagc = gammagc_k[i]
      gal_eps = galeps_k[i]

      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)

      gammagc = 0.48 * pi
      gal_eps0 = galeps_0.48[i]
      quotient_vevh = abs(gal_vevH^2 - gal_eps0^2) / (gal_eps0^2 + gal_vevH^2)
      g_sys0_pred = asin((sin(gammaU)^2 + (sin(gammagc)^2 - sin(gammaU)^2) * quotient_vevh)^0.5)

      if (i != 58)
      {
        g1_pred = g_sys1_pred / cos(g_sys1_pred)
        g0_pred = g_sys0_pred / cos(g_sys0_pred)
        v_obs   = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (gg000 * T0))^(1 / 4)
        v0_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g0_pred * T0))^(1 / 4)
        v1_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g1_pred * T0))^(1 / 4)

        ylim = c(-0.5, 0.5)
        plot(0, 0, ylim = ylim, xlim = c(-12, -8), pch = 20, axes = FALSE)
        lines(log10(ak000), log10(v0_pred^2 / smcGa07$Vobs[lg_gal]^2), col = "indianred")
        lines(log10(ak000), log10(v1_pred^2 / smcGa07$Vobs[lg_gal]^2), col = "blue")
        text(-9, 0.9 * ylim[2], names(namesg)[i], cex = 0.75)
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
