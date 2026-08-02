# =============================================================================
# Figure 2 - distinct acceleration relations of galaxies and galaxy clusters.
#
# Part A (top)   : rotation curves and the predicted/observed acceleration ratio
#                  of galaxies, with the MOND interpolating functions overlaid.
# Part B (bottom): acceleration excess vs. Newton-Hubble ratio for clusters and
#                  galaxies, with the MOND-like and general HMG model bands.
# Output: outputs/Fig2.pdf
# =============================================================================

{
  pdf(paste0("outputs/Fig2.pdf"), width = 7, height = 8.1)
  {
    # ------------------------------- Part A -------------------------------
    {
      smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
      smcGa07.V  = (smcGa07.Vk^2 * kms^2 + 2 * smcGa07.Vk^2 * kms^2 * smcGa07$R * kpc * c0 / (9 * T0))^0.25 / kms
      smcGa07.VO = aggregate(smcGa07.V, by = list(namess), mean, na.rm = T)

      # --- Rotation curves with the fitted HMG models ---
      par(fig = c(0, 0.495, 0.540, 1), mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3))
      for (j in c(1:length(namesg)))
      {
        i = order(smcGa07.VO$x, decreasing = TRUE)[j]
        lg_gal = namess == names(namesg)[i]
        smcGa07.Rg[lg_gal] = max(smcGa07$R[lg_gal])

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

        if (j != 55)
        {
          g1_pred = g_sys1_pred / cos(g_sys1_pred)
          g0_pred = g_sys0_pred / cos(g_sys0_pred)
          v_obs   = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (gg000 * T0))^(1 / 4)
          v0_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g0_pred * T0))^(1 / 4)
          v1_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g1_pred * T0))^(1 / 4)

          entre_uno = TRUE
          if (entre_uno)  divido = rep(1, length(lg_gal))
          if (!entre_uno) divido = smcGa07$Vobs
          ylim = c(0.05, 325)

          if (j == 1)
          {
            plot(smcGa07$R[lg_gal], smcGa07$Vobs[lg_gal] / divido[lg_gal], log = "x", ylim = ylim, xlim = c(0.2, 70), pch = 0, lwd = 0.5, cex = 0.6381, axes = FALSE)
            box()
            abline(h = seq(0, 300, 50), lty = 2, col = "gray95")
            abline(v = c(0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500), lty = 2, col = "gray95")
            axis(side = 2, seq(0, 300, 50), padj = 0.8)
            axis(side = 3, c(1, 2, 5, 10, 20, 50, 100, 200, 500), padj = 0.8)
            axis(side = 3, c(0.1, 0.2, 0.5), padj = 0.8)
          }
          points(smcGa07$R[lg_gal], smcGa07$Vobs[lg_gal] / divido[lg_gal], pch = 0, cex = 0.6381, lwd = 0.5, col = rainbow(length(namesg))[j])
          lines(smcGa07$R[lg_gal], v1_pred / divido[lg_gal], col = rainbow(length(namesg))[j])
          lines(smcGa07$R[lg_gal], v0_pred / divido[lg_gal], col = rainbow(length(namesg))[j], lty = 2)
        }
      }

      legend("topleft", legend = c("Observed rotation curves", "1-parameter HMG   ", "2-parameter HMG"), col = c(rgb(0.6, 0.6, 0.6, 0.9), rgb(0.5, 0.5, 0.5, 0.9)),
             lwd = c(0, 1, 1), lty = c(0, 2, 1), pch = c(0, 0, 0), pt.cex = c(1, 0, 0), cex = 0.801, bty = "n")

      par(fig = c(0, 0.495, 0.540, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3))
      plot.new()
      mtext(side = 2, expression(paste(italic(v)[obs], " [km/s]")), line = 1.7)
      mtext(side = 3, expression(paste(italic(R), " [kpc]")), line = 1.6)

      # --- Predicted/observed acceleration ratio vs. g_N, with MOND curves ---
      interpol_ref    = seq(-12, -8, 0.5)
      interpol_funct1 = array(NA, dim = c(length(namesg), length(interpol_ref)))
      interpol_funct2 = array(NA, dim = c(length(namesg), length(interpol_ref)))

      par(fig = c(0.505, 1, 0.540, 1), mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3), new = TRUE)
      for (j in 1:length(namesg))
        if (j != 55)
        {
          i = order(smcGa07.VO$x, decreasing = TRUE)[j]
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

          g1_pred = g_sys1_pred / cos(g_sys1_pred)
          g0_pred = g_sys0_pred / cos(g_sys0_pred)
          v_obs   = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (gg000 * T0))^(1 / 4)
          v0_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g0_pred * T0))^(1 / 4)
          v1_pred = smcGa07.Vk[lg_gal] * (1 + 1 / ak000 * 2 * c0 / (g1_pred * T0))^(1 / 4)

          if (j == 1)
          {
            plot(0, 0, ylim = c(-0.29, 0.105), xlim = c(-12, -8.8), pch = 20, axes = FALSE)
            axis(side = 4, round(seq(-0.6, 0.4, 0.1), 1), padj = -0.8)
            axis(side = 3, seq(-12, -8, 1), padj = 0.8)
            box()
            abline(h = 0, lty = 3)
            abline(h = round(seq(-0.2, 0.4, 0.1), 1), lty = 2, col = "gray95")
            abline(v = seq(-12, -8, 1), lty = 2, col = "gray95")
          }
          lines(log10(ak000), log10(v1_pred^2 / smcGa07$Vobs[lg_gal]^2), col = rgb(0.5, 0.5, 0.5, 0.1))
          lines(log10(ak000), log10(v0_pred^2 / smcGa07$Vobs[lg_gal]^2), col = rgb(0.6, 0.6, 0.6, 0.1), lty = 2)

          for (k in 1:length(interpol_ref))
          {
            lgk = log10(ak000) >= interpol_ref[k] - 0.25 & log10(ak000) < interpol_ref[k] + 0.25
            interpol_funct1[i, k] = mean(log10(v0_pred^2 / smcGa07$Vobs[lg_gal]^2)[lgk])
            interpol_funct2[i, k] = mean(log10(v1_pred^2 / smcGa07$Vobs[lg_gal]^2)[lgk])
          }
        }

      # MOND interpolating functions with 1-sigma envelopes.
      lines(log10_gN_values, log10_nu_simple_values - log10_nu_obs_values, col = rgb(0.9, 0.9, 0.2, 0.8), lwd = 3.5)
      lines(log10_gN_values, log10_nu_simple_values - log10_nu_obs_values - sigma_dex_values, col = rgb(0.9, 0.9, 0.2, 0.8), lwd = 1)
      lines(log10_gN_values, log10_nu_simple_values - log10_nu_obs_values + sigma_dex_values, col = rgb(0.9, 0.9, 0.2, 0.8), lwd = 1)
      lines(log10_gN_values, log10_nu_MLS_values - log10_nu_obs_values, col = rgb(0.3, 0.3, 0.9, 0.75), lwd = 2.3)
      lines(log10_gN_values, log10_nu_MLS_values - log10_nu_obs_values - sigma_dex_values, col = rgb(0.3, 0.3, 0.9, 0.75), lwd = 1)
      lines(log10_gN_values, log10_nu_MLS_values - log10_nu_obs_values + sigma_dex_values, col = rgb(0.3, 0.3, 0.9, 0.75), lwd = 1)
      lines(log10_gN_values, log10_nu_standard_values - log10_nu_obs_values, col = rgb(0.8, 0.2, 0.2, 0.5), lwd = 2.3)
      lines(log10_gN_values, log10_nu_standard_values - log10_nu_obs_values - sigma_dex_values, col = rgb(0.8, 0.2, 0.2, 0.5), lwd = 1)
      lines(log10_gN_values, log10_nu_standard_values - log10_nu_obs_values + sigma_dex_values, col = rgb(0.8, 0.2, 0.2, 0.5), lwd = 1)
      lines(log10_gN_values, log10_nu_sharp_values - log10_nu_obs_values, col = rgb(0.2, 0.8, 0.2, 0.5), lwd = 2.3)
      lines(log10_gN_values, log10_nu_sharp_values - log10_nu_obs_values - sigma_dex_values, col = rgb(0.2, 0.8, 0.2, 0.5), lwd = 1)
      lines(log10_gN_values, log10_nu_sharp_values - log10_nu_obs_values + sigma_dex_values, col = rgb(0.2, 0.8, 0.2, 0.5), lwd = 1)

      lines(interpol_ref, apply(interpol_funct1, 2, quantile, 0.17, na.rm = T), col = "gray45", lwd = 1.3, lty = 2)
      lines(interpol_ref, apply(interpol_funct1, 2, quantile, 0.50, na.rm = T), col = "gray45", lwd = 2, lty = 2)
      lines(interpol_ref, apply(interpol_funct1, 2, quantile, 0.84, na.rm = T), col = "gray45", lwd = 1.3, lty = 2)
      lines(interpol_ref, apply(interpol_funct2, 2, quantile, 0.17, na.rm = T), col = "gray45", lwd = 1)
      lines(interpol_ref, apply(interpol_funct2, 2, quantile, 0.50, na.rm = T), col = "gray45", lwd = 3)
      lines(interpol_ref, apply(interpol_funct2, 2, quantile, 0.84, na.rm = T), col = "gray45", lwd = 1)

      legend("bottomleft", legend = c("1-parameter HMG", "2-parameter HMG", "MLS", "Simple", "Standard", "Sharp"), ncol = 2, lwd = 2, lty = c(2, 1, 1, 1, 1, 1),
             col = c("gray45", "gray45", rgb(0.3, 0.3, 0.9, 0.75), rgb(0.9, 0.9, 0.2, 0.8), rgb(0.8, 0.2, 0.2, 0.5), rgb(0.2, 0.8, 0.2, 0.5)), cex = 0.801, bty = "n")

      par(fig = c(0.505, 1, 0.555, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3))
      plot.new()
      mtext(side = 4, expression(paste(log[10](italic(a)[pred] / italic(a)[obs]))), line = 1.7)
      mtext(side = 3, expression(paste(log[10](italic(a)[N]), " [", m / s^2, "]")), line = 1.6)
    }

    # ------------------------------- Part B -------------------------------
    modo = "nomed"
    eps0    = quantile(c(eps_clusB[, 1]), 0.5, na.rm = TRUE)
    eps0_do = quantile(c(eps_clusB[, 1]), 0.1, na.rm = TRUE)
    eps0_up = quantile(c(eps_clusB[, 1]), 0.9, na.rm = TRUE)

    uno = rep(1, length(vE2vH2_sort2))
    par(fig = c(0.0, 0.495, 0.07, 0.530), mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3), new = TRUE)
    plot(rar.escape_newton / rar.escape_hubble / eps0, diff_accel, col = "white", xlab = "", ylab = "", axes = FALSE, log = "xy", ylim = c(0.005, 1.5), xlim = c(0.02, 100))
    polygon_(sqrt(vE2vH2_sort2) / eps0, uno * cos(0.466 * pi) / (0.466 * pi), uno * cos(0.477 * pi) / (0.477 * pi), uno * cos(0.456 * pi) / (0.456 * pi), add = TRUE, col = rgb(0.12, 0.65, 0.99, 0.08))
    polygon_(sqrt(vE2vH2_sort2) / eps0, 1 / gamma_0mdd2, 1 / gamma_0mdo2, 1 / gamma_0mup2, add = TRUE, col = rgb(0.75, 0.75, 0.75, 0.3))
    box()
    clusters = as.character(unique(rar$Name))
    clus_col = rainbow(length(clusters))
    mtext(side = 1, expression(paste("Newton-Hubble ratio [", sqrt(2) * v[N] / (epsilon[H] * v[H]), "]")), line = 1.8)
    mtext(side = 2, expression(paste("[ ", italic(a)[Tot] - italic(a)[N], " ]/[ c/t ]")), line = 1.6)
    axis(side = 1, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200), paste(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200)), padj = -0.8)
    axis(side = 1, c(1), paste(c(1)), padj = -0.8)
    axis(side = 2, c(0.5, 1), paste(c(0.5, 1)), padj = 0.8)
    axis(side = 2, c(1, 2, 5, 10, 20) / 100, padj = 0.8)
    abline(v = c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200), h = c(1, 2, 5, 10, 20, 50, 100) / 100, lty = 2, col = "gray95")
    abline(v = c(1), lty = 2, col = "gray78")
    abline(h = c(1), lty = 2, col = "gray78")
    abline(h = c(0.5), lty = 2, col = "gray78")

    for (iclu in 1:length(clusters))
    {
      eps_H0    = eps_clusB[iclu, 1]
      eps_H0_do = eps_clusB[iclu, 2]
      eps_H0_up = eps_clusB[iclu, 3]
      if (modo == "med")
      {
        eps_H0 = eps0; eps_H0_do = eps0_do; eps_H0_up = eps0_up
      }
      lgclus = rar$Name == clusters[iclu]
      points((rar.escape_newton / rar.escape_hubble)[lgclus] / eps_H0, diff_accel[lgclus], pch = 20, cex = cex0, col = clus_col[iclu])
      segments((rar.escape_newton / rar.escape_hubble)[lgclus] / eps_H0, diff_accel_do[lgclus], (rar.escape_newton / rar.escape_hubble)[lgclus] / eps_H0, diff_accel_up[lgclus], col = clus_col[iclu])
      segments(((rar.escape_newton - rar.acc_newton_err) / rar.escape_hubble)[lgclus] / eps_H0, diff_accel[lgclus], ((rar.escape_newton + rar.acc_newton_err) / rar.escape_hubble)[lgclus] / eps_H0, diff_accel[lgclus], col = clus_col[iclu])
    }

    nonan = !is.na(gal_xi2[, 1]) & galeps_ko < 150
    for (j in c(1:length(namesg))[nonan])
    {
      i = order(smcGa07.VO$x, decreasing = TRUE)[j]
      lg_gal = namess == names(namesg)[i]
      x = sqrt(2 * mk00 / (smcGa07$R^3 * kpc^3) * T0^2)
      x_med = quantile(x[lg_gal], 0.5, na.rm = T)
      x_up  = quantile(x[lg_gal], 0.67, na.rm = T)
      x_do  = quantile(x[lg_gal], 0.33, na.rm = T)
      g0_med = quantile(cos(gamma1[lg_gal]) / (gamma1[lg_gal]), 0.5, na.rm = T)
      g0_up  = quantile(cos(gamma1[lg_gal]) / (gamma1[lg_gal]), 0.67, na.rm = T)
      g0_do  = quantile(cos(gamma1[lg_gal]) / (gamma1[lg_gal]), 0.33, na.rm = T)
      eps0_gal = galeps_ko[i]
      if (modo == "med") eps0_gal = eps0
      points(x_med / eps0_gal, g0_med, col = rainbow(length(namesg))[j], pch = 0, cex = 0.6381)
      segments(x0 = x_med / eps0_gal, x1 = x_med / eps0_gal, y0 = g0_do, y1 = g0_up, col = rainbow(length(namesg))[j])
      segments(x0 = x_do / eps0_gal, x1 = x_up / eps0_gal, y0 = g0_med, y1 = g0_med, col = rainbow(length(namesg))[j])
    }

    lines(sqrt(vE2vH2_sort2) / eps0, uno * cos(0.466 * pi) / (0.466 * pi), lwd = 2, col = rgb(0.12, 0.65, 0.99, 0.2), lty = 1)
    lines(sqrt(vE2vH2_sort2) / eps0, 1 / gamma_0mdd2, lwd = 1.8, col = "gray70", lty = 1)
    rect(xleft = 3.5, ybottom = 0.037, xright = 69.9, ytop = 0.14, border = "lightblue3")

    # --- Zoomed inset of the transition regime ---
    par(fig = c(0.505, 1.0, 0.07, 0.530), mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3), new = TRUE)
    plot(rar.escape_newton / rar.escape_hubble / eps0, diff_accel, col = "white", xlab = "", ylab = "", axes = FALSE, log = "xy", ylim = c(0.04, 0.12), xlim = c(4, 69.9))
    polygon_(sqrt(vE2vH2_sort2) / eps0, uno * cos(0.466 * pi) / (0.466 * pi), uno * cos(0.477 * pi) / (0.477 * pi), uno * cos(0.456 * pi) / (0.456 * pi), add = TRUE, col = rgb(0.12, 0.65, 0.99, 0.08))
    polygon_(sqrt(vE2vH2_sort2) / eps0, 1 / gamma_0mdd2, 1 / gamma_0mdo2, 1 / gamma_0mup2, add = TRUE, col = rgb(0.75, 0.75, 0.75, 0.3))
    clusters = as.character(unique(rar$Name))
    clus_col = rainbow(length(clusters))
    mtext(side = 1, expression(paste("Newton-Hubble ratio [", sqrt(2) * v[N] / (epsilon[H] * v[H]), "]")), line = 1.8)
    mtext(side = 4, expression(paste("[ ", italic(a)[Tot] - italic(a)[N], " ]/[ c/t ]")), line = 1.7)
    axis(side = 1, c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 50, 100, 200), paste(c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 50, 100, 200)), padj = -0.8)
    axis(side = 4, c(1, 2, 4, 5, 7, 10, 20, 50, 100) / 100, padj = -0.8)
    abline(v = c(0.01, 0.02, 0.04, 0.05, 0.07, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100, 200), h = c(1, 2, 4, 5, 7, 10, 20, 50, 100) / 100, lty = 2, col = "gray95")
    box(col = "lightblue3")
    lines(sqrt(vE2vH2_sort2) / eps0, uno * cos(0.466 * pi) / (0.466 * pi), lwd = 2, col = rgb(0.12, 0.65, 0.99, 0.2), lty = 1)
    lines(sqrt(vE2vH2_sort2) / eps0, 1 / gamma_0mdd2, lwd = 2, col = "gray70", lty = 1)

    x_esc = rep(NA, length(namesg))
    y_acc = rep(NA, length(namesg))
    e_acc = rep(NA, length(namesg))
    for (j in c(1:length(namesg)))
    {
      i = order(smcGa07.VO$x, decreasing = TRUE)[j]
      lg_gal = namess == names(namesg)[i]
      x = sqrt(2 * mk00 / (smcGa07$R^3 * kpc^3) * T0^2)
      x_med = quantile(x[lg_gal], 0.5, na.rm = T)
      x_up  = quantile(x[lg_gal], 0.67, na.rm = T)
      x_do  = quantile(x[lg_gal], 0.33, na.rm = T)
      g0_med = quantile(cos(gamma1[lg_gal]) / (gamma1[lg_gal]), 0.5, na.rm = T)
      g0_up  = quantile(cos(gamma1[lg_gal]) / (gamma1[lg_gal]), 0.67, na.rm = T)
      g0_do  = quantile(cos(gamma1[lg_gal]) / (gamma1[lg_gal]), 0.33, na.rm = T)
      g0_err = (g0_up - g0_do) / 2
      x_err  = (x_up - x_do) / 2
      eps0_gal = galeps_ko[i]
      if (modo == "med") eps0_gal = eps0
      x_esc[i] = x_med / eps0_gal
      y_acc[i] = g0_med
      e_acc[i] = g0_err * x_err
      points(x_med / eps0_gal, g0_med, col = rainbow(length(namesg))[j], lwd = 1.7, pch = 0, cex = 0.6781)
      segments(x0 = x_med / eps0_gal, x1 = x_med / eps0_gal, y0 = g0_do, y1 = g0_up, col = rainbow(length(namesg))[j], lty = 1)
      segments(x0 = x_do / eps0_gal, x1 = x_up / eps0_gal, y0 = g0_med, y1 = g0_med, col = rainbow(length(namesg))[j], lty = 1)
    }

    legend("topright", legend = c("MOND-like HMG", "General HMG prediction"), lwd = c(6.8, 6.8, 0.0), lty = c(1, 1, 2), col = c(rgb(0.12, 0.65, 0.99, 0.08), rgb(0.75, 0.75, 0.75, 0.3), "gray65"), cex = 0.80, bty = "n")
    legend("topright", legend = c("MOND-like HMG", "General HMG prediction"), lwd = c(1.8, 1.8, 1.2), lty = c(1, 1, 2), col = c(rgb(0.12, 0.65, 0.99, 0.2), "gray75", "gray45"), cex = 0.80, bty = "n")

    # --- Legends ---
    par(fig = c(0, 1, 0, 0.101), mar = c(0, 0, 0, 0), oma = c(0.01, 0.1, 1, 0.1), new = TRUE)
    plot.new()
    legend("bottomleft", legend = clusters, pch = 20, col = clus_col, ncol = 2, pt.cex = cex0 * 0.61, cex = cex0 * 0.45, bty = "n")
    legend("bottomright", legend = names(namesG)[order(smcGa07.VO$x, decreasing = TRUE)[-55]], pch = 0, col = rainbow(length(namesg))[-55], ncol = 10, pt.cex = cex0 * 0.48, cex = cex0 * 0.365, bty = "n")

    par(fig = c(0, 1, 0, 0.125), mar = c(0, 0, 0, 0), oma = c(0, 0.1, 0, 0.45), new = T)
    plot.new()
    legend("topleft", cex = 0.801, legend = expression(bold("Clusters:")), pch = 20, col = "white", bty = "n")
    legend("topright", cex = 0.801, legend = expression(bold("Galaxies:")), pch = 20, col = "white", bty = "n", pt.cex = 0.001)

    par(fig = c(0, 1, 0.11, 0.530), mar = c(0, 0, 0, 0), oma = c(5, 3, 3, 3), new = T)
    Flecha(xmin = 0.465, xmax = 0.51, ymin = 0.29, ymax = 0.54, mode = "right", dx = 0.01, dy = 0.2, px = 0.43, py = 0.5, color = rgb(0.97, 0.97, 1), border = "lightblue3", lwd = 0.7, new = TRUE)
  }
  dev.off()
}
