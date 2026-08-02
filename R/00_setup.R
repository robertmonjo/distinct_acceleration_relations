# =============================================================================
# Stage 0 - packages, physical constants, model helpers and data preprocessing
#
# Monjo & Banik (2025), ApJ 992, 35. Sourced first by run_all.R; every later
# stage relies on the objects created here (shared session state).
# =============================================================================

library("jpeg")
library("png")
library("fields")
library("stringr")

cex0 <- 1.5   # global point/label scale used across the figures


# -----------------------------------------------------------------------------
# Helper: shaded arrow polygon (used to annotate the transition regime).
# -----------------------------------------------------------------------------
Flecha = function(xmin, xmax, ymin, ymax, mode, dx = 0.1, dy = 0.1,
                  px = 0.5, py = 0.5, color = "lightblue", border = "blue",
                  new = TRUE, ...)
{
  Dx = xmax - xmin
  Dy = ymax - ymin
  if (new) plot.new()

  if (mode == "right")
    polygon(c(xmin, xmin + (px - dx) * Dx, xmin + px * Dx, xmax, xmin + px * Dx, xmin + (px - dx) * Dx, xmin),
            c(ymin + dy * Dy, ymin + dy * Dy, ymin, (ymax + ymin) / 2, ymax, ymax - dy * Dy, ymax - dy * Dy), col = color, border = border, ...)
  if (mode == "left")
    polygon(c(xmax, xmax - (px - dx) * Dx, xmax - px * Dx, xmin, xmax - px * Dx, xmax - (px - dx) * Dx, xmax),
            c(ymax - dy * Dy, ymax - dy * Dy, ymax, (ymax + ymin) / 2, ymin, ymin + dy * Dy, ymin + dy * Dy), col = color, border = border, ...)
  if (mode == "up")
    polygon(c(xmin + dx * Dx, xmin + dx * Dx, xmin, (xmax + xmin) / 2, xmax, xmax - dx * Dx, xmax - dx * Dx),
            c(ymin, ymin + (py - dy) * Dy, ymin + py * Dy, ymax, ymin + py * Dy, ymin + (py - dy) * Dy, ymin), col = color, border = border, ...)
  if (mode == "down")
    polygon(c(xmin + dx * Dx, xmin + dx * Dx, xmin, (xmax + xmin) / 2, xmax, xmax - dx * Dx, xmax - dx * Dx),
            c(ymax, ymax - (py - dy) * Dy, ymax - py * Dy, ymin, ymax - py * Dy, ymax - (py - dy) * Dy, ymax), col = color, border = border, ...)
}


# -----------------------------------------------------------------------------
# Helper: projective angle gamma from the ratio gamma0 = gamma/cos(gamma).
# Numerical inversion of gamma0 = (2/pi) / pi_gg via fixed-point iteration.
# -----------------------------------------------------------------------------
acos. = function(pi_gg)
{
  ggg = (1 / pi_gg) * (pi / 2)
  gamma1 = atan(ggg)
  gamma1[ggg > 0 & !is.na(ggg)] = ggg[ggg > 0 & !is.na(ggg)]^0.68 / 1.34
  gamma2 = gamma1
  gamma3 = gamma1
  gamma4 = gamma1
  gamma2[ggg > 1 & !is.na(ggg)] = acos((gamma1 / ggg)[ggg > 1 & !is.na(ggg)])
  gamma3[ggg > 1 & !is.na(ggg)] = acos((gamma2 / ggg)[ggg > 1 & !is.na(ggg)])
  gamma4[ggg > 1 & !is.na(ggg)] = acos(cos(gamma3[ggg > 1 & !is.na(ggg)]) / 0.91)
  return(gamma4)
}


# -----------------------------------------------------------------------------
# Helper: filled band between two curves (confidence envelopes).
# -----------------------------------------------------------------------------
polygon_ = function(x, v, y1, y2, GRID = F, cols = "gray", dens = NULL, add = T)
{
  years = x; v1 = v
  xI = sort(x, decreasing = T)
  y2I = y2; for (i in 1:length(y2)) y2I[i] = y2[length(y2) - i + 1]
  lg = is.na(y1) == F; nye = length(x[lg])
  lgI = is.na(y2I) == F
  xx = c(x[lg], xI[lgI])
  yy = c(y1[lg], y2I[lgI])
  if (add == F)
    plot(xx, yy * NA, ylim = c(min(yy, na.rm = T), max(yy, na.rm = T)))
  polygon(xx, yy, col = cols, border = cols, anomaly = dens)
  if (GRID) grid(lwd = cex0 * c(1, 1, 2.7, 2.7)[f])
}


# -----------------------------------------------------------------------------
# Physical constants (SI units unless stated otherwise).
# -----------------------------------------------------------------------------
c0   = 3 * 10^8                          # speed of light [m/s]
Msol = 1.9891e30                         # solar mass [kg]
ua   = 1.496e+11                         # astronomical unit [m]
kpc  = 3261.8116478174 * 365 * 24 * 3600 * c0   # kiloparsec [m]
pc   = kpc / 1000                        # parsec [m]
kms  = 1000                              # km/s in [m/s]
T0   = 13.7 * 10^9 * 365 * 24 * 3600     # cosmic time t = 1/H0 [s]
GN   = 6.674e-11                         # gravitational constant [m^3 kg^-1 s^-2]


# -----------------------------------------------------------------------------
# Binned galaxy radial acceleration relation (used as a reference band later).
# -----------------------------------------------------------------------------
data = read.table("data/RAR_data_Indranil.txt")
log10_gN_values     = data[, 1]
log10_g_values      = data[, 2]
Dispersion          = data[, 3]
N_data              = data[, 4]
sigma_dex_values    = Dispersion / sqrt(N_data)
log10_nu_obs_values = log10_g_values - log10_gN_values

# MOND interpolating functions (drawn as reference curves in Figure 2).
a_0 = 1.2e-10
nu_simple   = 0.5 + sqrt(0.25 + a_0 / 10^log10_gN_values)
nu_standard = sqrt(0.5 + sqrt(0.25 + a_0 * a_0 / 100^log10_gN_values))
nu_sharp    = sqrt(a_0 / 10^log10_gN_values)
nu_sharp[nu_sharp < 1.0] = 1.0
nu_MLS_values = 1.0 / (1 - exp(-sqrt((10^log10_gN_values) / a_0)))
log10_nu_simple_values   = log10(nu_simple)
log10_nu_standard_values = log10(nu_standard)
log10_nu_sharp_values    = log10(nu_sharp)
log10_nu_MLS_values      = log10(nu_MLS_values)


# -----------------------------------------------------------------------------
# McGaugh et al. (2007, ApJ 659, 149) galaxy mass models.
#   mcgaugh  : per-radius rotation-curve table (R, Vobs, Vgas, Vst, ...)
#   mcgaugh_md   : companion mass-discrepancy table
# -----------------------------------------------------------------------------
mcgaugh_md  = read.table("data/S_McGaugh2007_mass_discrepancy.txt", sep = "\t", header = T)
mcgaugh = read.table("data/McGaugh2007.txt", sep = "\t", header = T)
mcgaugh_md  = mcgaugh_md[!is.na(mcgaugh$R), ]
gal_name_row  = mcgaugh$Name
md_names  = mcgaugh_md$Name
gal_counts  = table(gal_name_row)
gal_names_up  = gal_counts
names(gal_names_up) = str_to_upper(names(gal_counts))
names(gal_names_up)[c(1, 2, 3)] = c("F563-1", "F563-V2", "F568-V1")

# Baryonic (Vst + Vgas) circular speed and derived accelerations.
mcgaugh.Vk = sqrt(mcgaugh$Vst^2 + mcgaugh$Vgas^2)
mcgaugh.V  = (mcgaugh.Vk^2 * kms^2 + 2 * mcgaugh.Vk^2 * kms^2 * mcgaugh$R * kpc * c0 / (9 * T0))^0.25 / kms
a_newton = mcgaugh.Vk^2 / (mcgaugh$R * kpc) * kms^2       # Newtonian (baryonic) acceleration
aobs0 = mcgaugh$Vobs^2 / (mcgaugh$R * kpc) * kms^2     # observed acceleration
vel_ratio_sq = mcgaugh$Vobs^2 / mcgaugh.Vk^2                  # squared velocity ratio
mass_bar = mcgaugh.Vk^2 * (mcgaugh$R * kpc) * kms^2       # enclosed baryonic mass proxy
mass_bar_pt  = mass_bar

# Per-point projective-angle factor gamma1 (needed by the galaxy branch of Fig. 2).
a_newton_pt = mcgaugh.Vk^2 / (mcgaugh$R * kpc) * kms^2
vel_ratio_sq_pt = mcgaugh$Vobs^2 / mcgaugh.Vk^2
gg = 1 / a_newton_pt * (2 * c0 / T0) / (vel_ratio_sq_pt^2 - 1)
gg[gg > 30] = NA
gg[gg < 1]  = NA
a_newton_pt[(a_newton_pt / ((mcgaugh$R * kpc) / T0^2))^0.25 > 30] = NA
gg = round(gg, 2)
gamma1 = acos.((pi / 2) / gg)
gamma1[gamma1 < 1.4] = NA
g_ratio = 1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2 - 1)
gal_vevH = (sqrt(2) * mcgaugh.Vk * kms * T0) / (mcgaugh$R * kpc)

# Per-galaxy aggregates.
mcgaugh.R    = aggregate(mcgaugh$R, by = list(gal_name_row), max, na.rm = T)
mcgaugh.VN   = aggregate(mcgaugh.Vk, by = list(gal_name_row), max, na.rm = T)
mcgaugh.VO   = aggregate(mcgaugh.V, by = list(gal_name_row), mean, na.rm = T)
mcgaugh.mass = mcgaugh.VN$x^2 * mcgaugh.R$x * kpc * kms^2 / GN


# -----------------------------------------------------------------------------
# Per-galaxy fit of the neighbourhood parameter (chi-square over a grid of
# eps_H = seq_galeps and central angle gamma_cen = seq_gammagc). Results feed
# the galaxy-parameter stage and Figures 2 and 5.
# -----------------------------------------------------------------------------
seq_gammagc = seq(0.44, 0.5, 0.001) * pi
seq_galeps  = unique(c(1, seq(1, 10, 0.1), seq(10, 200, 1)))

gal_chi2_0.47 = array(NA, dim = c(length(names(gal_counts)), length(seq_galeps)),
                     dimnames = list(names(gal_counts), seq_galeps))
gal_chi2_0.48 = array(NA, dim = c(length(names(gal_counts)), length(seq_galeps)),
                     dimnames = list(names(gal_counts), seq_galeps))
gal_chi2_0.5  = array(NA, dim = c(length(names(gal_counts)), length(seq_galeps)),
                     dimnames = list(names(gal_counts), seq_galeps))
gal_chi22     = array(NA, dim = c(length(names(gal_counts)), length(seq_galeps), length(seq_gammagc)),
                     dimnames = list(names(gal_counts), seq_galeps, seq_gammagc))

mcgaugh.Rg = rep(NA, length(mcgaugh$R))
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
  gal_vevH = (sqrt(2) * mcgaugh.Vk[lg_gal] * kms * T0) / (mcgaugh$R[lg_gal] * kpc)
  gammaU = pi / 3

  g_ratio_btfr = mean(1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2), na.rm = T)
  g_ratio = 1 / a_newton * (2 * c0 / T0) / (vel_ratio_sq^2 - 1)
  g_ratio = round(g_ratio, 2)
  g_sys1 = acos.((pi / 2) / g_ratio)
  g_sys1[g_sys1 < 1] = NA

  # One-parameter (general) model at three fixed central angles.
  for (k in 1:length(seq_galeps))
  {
    gaempty = pi / 3
    gal_eps = seq_galeps[k]
    quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)

    gamma_cen = 0.47 * pi
    g_sys0_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
    g_ratio_pred = g_sys0_pred / cos(g_sys0_pred)
    gal_chi2_0.48[i, k] = mean(abs(log(g_ratio_pred) - log(g_ratio)), na.rm = T)

    gamma_cen = 0.48 * pi
    g_sys0_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
    g_ratio_pred = g_sys0_pred / cos(g_sys0_pred)
    gal_chi2_0.48[i, k] = mean(abs(log(g_ratio_pred) - log(g_ratio)), na.rm = T)

    gamma_cen = 0.50 * pi
    g_sys0_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
    g_ratio_pred = g_sys0_pred / cos(g_sys0_pred)
    gal_chi2_0.5[i, k] = mean(abs(log(g_ratio_pred) - log(g_ratio)), na.rm = T)
  }

  # Two-parameter (cluster-like) model over the (eps_H, gamma_cen) grid.
  for (k in 1:length(seq_galeps))
    for (j in 1:length(seq_gammagc))
    {
      gamma_cen = seq_gammagc[j]
      gaempty = pi / 3
      gal_eps = seq_galeps[k]
      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gamma_cen)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
      g_ratio_gal_pred = g_sys1_pred / cos(g_sys1_pred)
      gal_chi22[i, k, j] = mean(abs(log(g_ratio_gal_pred) - log(g_ratio)), na.rm = T)
    }
}

# Best-fit containers filled in the next stage.
galaxy_eps_0.47 = rep(NA, length(gal_counts))
galaxy_eps_0.48 = rep(NA, length(gal_counts))
galaxy_eps_0.5  = rep(NA, length(gal_counts))
galaxy_eps_k    = rep(NA, length(gal_counts))
gamma_cen_k   = rep(NA, length(gal_counts))
