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
#   smcGa07  : per-radius rotation-curve table (R, Vobs, Vgas, Vst, ...)
#   mcGa07   : companion mass-discrepancy table
# -----------------------------------------------------------------------------
mcGa07  = read.table("data/S_McGaugh2007_mass_discrepancy.txt", sep = "\t", header = T)
smcGa07 = read.table("data/McGaugh2007.txt", sep = "\t", header = T)
mcGa07  = mcGa07[!is.na(smcGa07$R), ]
namess  = smcGa07$Name
namesm  = mcGa07$Name
namesg  = table(namess)
namesG  = namesg
names(namesG) = str_to_upper(names(namesg))
names(namesG)[c(1, 2, 3)] = c("F563-1", "F563-V2", "F568-V1")

# Baryonic (Vst + Vgas) circular speed and derived accelerations.
smcGa07.Vk = sqrt(smcGa07$Vst^2 + smcGa07$Vgas^2)
smcGa07.V  = (smcGa07.Vk^2 * kms^2 + 2 * smcGa07.Vk^2 * kms^2 * smcGa07$R * kpc * c0 / (9 * T0))^0.25 / kms
ak000 = smcGa07.Vk^2 / (smcGa07$R * kpc) * kms^2       # Newtonian (baryonic) acceleration
aobs0 = smcGa07$Vobs^2 / (smcGa07$R * kpc) * kms^2     # observed acceleration
di000 = smcGa07$Vobs^2 / smcGa07.Vk^2                  # squared velocity ratio
mk000 = smcGa07.Vk^2 * (smcGa07$R * kpc) * kms^2       # enclosed baryonic mass proxy
mk00  = mk000

# Per-point projective-angle factor gamma1 (needed by the galaxy branch of Fig. 2).
ak00 = smcGa07.Vk^2 / (smcGa07$R * kpc) * kms^2
di00 = smcGa07$Vobs^2 / smcGa07.Vk^2
gg = 1 / ak00 * (2 * c0 / T0) / (di00^2 - 1)
gg[gg > 30] = NA
gg[gg < 1]  = NA
ak00[(ak00 / ((smcGa07$R * kpc) / T0^2))^0.25 > 30] = NA
gg = round(gg, 2)
gamma1 = acos.((pi / 2) / gg)
gamma1[gamma1 < 1.4] = NA
gg000 = 1 / ak000 * (2 * c0 / T0) / (di000^2 - 1)
gal_vevH = (sqrt(2) * smcGa07.Vk * kms * T0) / (smcGa07$R * kpc)

# Per-galaxy aggregates.
smcGa07.R    = aggregate(smcGa07$R, by = list(namess), max, na.rm = T)
smcGa07.VN   = aggregate(smcGa07.Vk, by = list(namess), max, na.rm = T)
smcGa07.VO   = aggregate(smcGa07.V, by = list(namess), mean, na.rm = T)
smcGa07.mass = smcGa07.VN$x^2 * smcGa07.R$x * kpc * kms^2 / GN


# -----------------------------------------------------------------------------
# Per-galaxy fit of the neighbourhood parameter (chi-square over a grid of
# eps_H = seq_galeps and central angle gamma_cen = seq_gammagc). Results feed
# the galaxy-parameter stage and Figures 2 and 5.
# -----------------------------------------------------------------------------
seq_gammagc = seq(0.44, 0.5, 0.001) * pi
seq_galeps  = unique(c(1, seq(1, 10, 0.1), seq(10, 200, 1)))

gal_xi2_0.47 = array(NA, dim = c(length(names(namesg)), length(seq_galeps)),
                     dimnames = list(names(namesg), seq_galeps))
gal_xi2_0.48 = array(NA, dim = c(length(names(namesg)), length(seq_galeps)),
                     dimnames = list(names(namesg), seq_galeps))
gal_xi2_0.5  = array(NA, dim = c(length(names(namesg)), length(seq_galeps)),
                     dimnames = list(names(namesg), seq_galeps))
gal_xi22     = array(NA, dim = c(length(names(namesg)), length(seq_galeps), length(seq_gammagc)),
                     dimnames = list(names(namesg), seq_galeps, seq_gammagc))

smcGa07.Rg = rep(NA, length(smcGa07$R))
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
  gal_vevH = (sqrt(2) * smcGa07.Vk[lg_gal] * kms * T0) / (smcGa07$R[lg_gal] * kpc)
  gammaU = pi / 3

  gg000_btfr = mean(1 / ak000 * (2 * c0 / T0) / (di000^2), na.rm = T)
  gg000 = 1 / ak000 * (2 * c0 / T0) / (di000^2 - 1)
  gg000 = round(gg000, 2)
  g_sys1 = acos.((pi / 2) / gg000)
  g_sys1[g_sys1 < 1] = NA

  # One-parameter (general) model at three fixed central angles.
  for (k in 1:length(seq_galeps))
  {
    gaempty = pi / 3
    gal_eps = seq_galeps[k]
    quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)

    gammagc = 0.47 * pi
    g_sys0_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
    gg000_pred = g_sys0_pred / cos(g_sys0_pred)
    gal_xi2_0.48[i, k] = mean(abs(log(gg000_pred) - log(gg000)), na.rm = T)

    gammagc = 0.48 * pi
    g_sys0_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
    gg000_pred = g_sys0_pred / cos(g_sys0_pred)
    gal_xi2_0.48[i, k] = mean(abs(log(gg000_pred) - log(gg000)), na.rm = T)

    gammagc = 0.50 * pi
    g_sys0_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
    gg000_pred = g_sys0_pred / cos(g_sys0_pred)
    gal_xi2_0.5[i, k] = mean(abs(log(gg000_pred) - log(gg000)), na.rm = T)
  }

  # Two-parameter (cluster-like) model over the (eps_H, gamma_cen) grid.
  for (k in 1:length(seq_galeps))
    for (j in 1:length(seq_gammagc))
    {
      gammagc = seq_gammagc[j]
      gaempty = pi / 3
      gal_eps = seq_galeps[k]
      quotient_vevh = abs(gal_vevH^2 - gal_eps^2) / (gal_eps^2 + gal_vevH^2)
      g_sys1_pred = asin((sin(gaempty)^2 + (sin(gammagc)^2 - sin(gaempty)^2) * quotient_vevh)^0.5)
      gg001_pred = g_sys1_pred / cos(g_sys1_pred)
      gal_xi22[i, k, j] = mean(abs(log(gg001_pred) - log(gg000)), na.rm = T)
    }
}

# Best-fit containers filled in the next stage.
galeps_0.47 = rep(NA, length(namesg))
galeps_0.48 = rep(NA, length(namesg))
galeps_0.5  = rep(NA, length(namesg))
galeps_k    = rep(NA, length(namesg))
gammagc_k   = rep(NA, length(namesg))
