# =============================================================================
# Parameter economy of the HMG general model (supplementary analysis).
#
# The neighbourhood parameter is closed by the local baryonic density,
#     eps_H^2 = rho_nei(s) / rho_vac + 1/6,    rho_nei(s) = M_bar / [(4/3) pi (s r)^3]
# (Monjo 2026, MNRAS 549, stag965, Eq. 4), so eps_H is fixed by the observed mass
# and size rather than fitted. This script quantifies, on the paper's own data,
# the fit and the Bayesian information criterion of three descriptions:
#   HMG-A : eps_H closed by the density (s = 4 for galaxies)   -> 0 fitted params/object
#   HMG-B : eps_H fitted per object                            -> 1 fitted param /object
#   MOND  : simple interpolating function, a0 fixed            -> 0 fitted params
# and reports, for the ten clusters, the per-cluster scale s that places each on
# the density relation. Baryons are identical for every model (no M/L fitting).
#
# Output: outputs/Suppl_parameter_economy.txt
# Usage:  Rscript supplement/parameter_economy.R
# =============================================================================
.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.file) == 1) setwd(dirname(dirname(normalizePath(.file))))  # repo root

pdf(NULL)
source("R/00_setup.R")

a0      <- 1.2e-10                 # MOND acceleration scale [m/s^2]
gamma_cen <- 0.48 * pi            # galaxy central angle used in the paper
gamma_U   <- pi / 3
rho_vac   <- 3 / (8 * pi * GN * T0^2)

# --- Galaxies -----------------------------------------------------------------
gnames <- names(gal_counts); Ngal <- length(gnames)
# density ratio rho_nei/rho_vac at s = 4, per galaxy (outer point aggregates)
rho_ratio_gal <- mcgaugh.mass / (4/3 * pi * (mcgaugh.R$x * kpc * 4)^3) / rho_vac
epsH_closed   <- sqrt(rho_ratio_gal + 1/6)          # Eq. 4 closure

# Chi-square per galaxy under HMG (given eps_H) or MOND, with an assumed velocity
# uncertainty sigma = sqrt((FRAC*Vobs)^2 + FLOOR^2) [km/s] (the paper's data carry
# no per-point errors; FRAC = 0.05, FLOOR = 5 km/s give chi2_nu ~ 1 for the closure).
FRAC <- 0.05; FLOOR <- 5
gal_chi2 <- function(i, model, epsH = NA) {
  lg <- gal_name_row == gnames[i]
  Vk <- sqrt(mcgaugh$Vst[lg]^2 + mcgaugh$Vgas[lg]^2)
  R  <- mcgaugh$R[lg]; Vobs <- mcgaugh$Vobs[lg]
  aN <- Vk^2 / (R * kpc) * kms^2
  if (model == "MOND") {
    D <- 0.5 + sqrt(0.25 + a0 / aN)
  } else {
    vevH <- (sqrt(2) * Vk * kms * T0) / (R * kpc)
    quot <- abs(vevH^2 - epsH^2) / (epsH^2 + vevH^2)
    gsys <- asin(sqrt(sin(gamma_U)^2 + (sin(gamma_cen)^2 - sin(gamma_U)^2) * quot))
    D    <- sqrt(1 + (2 * c0 / T0) / (aN * gsys / cos(gsys)))
  }
  Vmod <- Vk * sqrt(D); sig <- sqrt((FRAC * Vobs)^2 + FLOOR^2)
  ok <- is.finite(Vmod) & is.finite(Vobs) & Vobs > 0
  list(chi2 = sum(((Vmod[ok] - Vobs[ok]) / sig[ok])^2), n = sum(ok))
}

grid_eps <- unique(c(seq(0.5, 10, 0.1), seq(10, 200, 0.5)))
N <- 0; used <- 0; chi2A <- 0; chi2B <- 0; chi2M <- 0; epsB <- rep(NA, Ngal)
for (i in 1:Ngal) {
  if (!is.finite(epsH_closed[i])) next
  a <- gal_chi2(i, "HMG", epsH_closed[i]); if (a$n == 0) next
  used <- used + 1; N <- N + a$n
  chi2A <- chi2A + a$chi2
  chi2M <- chi2M + gal_chi2(i, "MOND")$chi2
  sB <- sapply(grid_eps, function(e) gal_chi2(i, "HMG", e)$chi2)
  chi2B <- chi2B + min(sB); epsB[i] <- grid_eps[which.min(sB)]
}
bic <- function(chi2, k) chi2 + k * log(N)                   # standard BIC, known errors

gal_tab <- data.frame(
  model = c("HMG eps_H closed by density (s=4)", "HMG eps_H fitted per galaxy",
            "MOND (a0 fixed = 1.2e-10)"),
  fitted_params = c(0, used, 0),
  k = c(1, used + 1, 0),
  chi2 = round(c(chi2A, chi2B, chi2M)),
  chi2_nu = round(c(chi2A/(N-1), chi2B/(N-used), chi2M/N), 2),
  BIC = round(c(bic(chi2A, 1), bic(chi2B, used + 1), bic(chi2M, 0))))
gal_tab$dBIC_vs_MOND <- round(gal_tab$BIC - gal_tab$BIC[3])
r_pred_fit <- cor(epsH_closed, epsB, use = "complete.obs")

# --- Clusters -----------------------------------------------------------------
rar <- read.table("data/clusterRAR.dat")
colnames(rar) <- c("Name","z","radius","log(gbar)","log(gtot)","err_low","err_up")
clusters <- as.character(unique(rar$Name))
acc_tot <- 10^rar$`log(gtot)`; acc_new <- 10^rar$`log(gbar)`
acc_up  <- 10^(rar$`log(gtot)`+rar$err_up); acc_do <- 10^(rar$`log(gtot)`-rar$err_low)
quot     <- acc_tot / acc_new
quot_err <- ((acc_up - acc_do) / 2) / acc_tot * quot
rar.radius <- rar$radius * kpc
rar.mass   <- acc_new * rar.radius^2 / GN
esc_N <- sqrt(2 * acc_new * rar.radius); esc_H <- rar.radius / T0 * (1 + rar$z)
pred_b <- function(epsH) {                 # cluster model b (gamma_cen = pi/2)
  vEvH <- esc_N / esc_H
  dd <- abs(vEvH^2 - epsH^2) / (epsH^2 + vEvH^2)
  gM <- asin(sqrt(sin(pi/3)^2 - dd * (sin(pi/3)^2 - sin(pi/2)^2)))
  ((1 + rar$z) * c0) / (T0 * gM / cos(gM) * acc_new)
}
grid_c <- unique(c(seq(0.5,10,0.1), seq(10,200,0.5), seq(200,3000,5)))
eps_fit <- sapply(clusters, function(cl) { lg <- rar$Name == cl
  grid_c[which.min(sapply(grid_c, function(e) sum((pred_b(e)[lg]-quot[lg])^2/quot_err[lg]^2)))] })
rho_ratio_clu_s1 <- sapply(clusters, function(cl) { lg <- rar$Name == cl
  min(rar.mass[lg]) / (4/3 * pi * min(rar.radius[lg])^3) / (3/(8*pi*GN*(1+max(rar$z[lg]))^(-1)*T0^2)) })
# scale centring each cluster on eps_H^2 = rho_nei(s)/rho_vac + 1/6
s_centred <- (rho_ratio_clu_s1 / (eps_fit^2 - 1/6))^(1/3)
clu_tab <- data.frame(cluster = clusters, eps_fit = round(eps_fit,1),
                      s_centred = round(s_centred,2))

# --- Report -------------------------------------------------------------------
sink("outputs/Suppl_parameter_economy.txt")
cat("Parameter economy of the HMG general model\n")
cat("eps_H^2 = rho_nei(s)/rho_vac + 1/6   (Monjo 2026, MNRAS 549, Eq. 4)\n\n")
cat(sprintf("GALAXIES  (McGaugh 2007: %d galaxies, %d rotation-curve points)\n", used, N))
cat(sprintf("sigma_V = sqrt((%.2f Vobs)^2 + %d^2) km/s ; BIC = chi2 + k ln N\n\n", FRAC, FLOOR))
print(gal_tab, row.names = FALSE)
cat(sprintf("\nr(eps_H closed, eps_H fitted) = %.3f\n", r_pred_fit))
cat(sprintf("\nCLUSTERS  (%d clusters, %d RAR points): per-cluster scale s for the density relation\n\n",
            length(clusters), nrow(rar)))
print(clu_tab, row.names = FALSE)
cat(sprintf("\ns: median = %.2f, range %.2f - %.2f  (galaxies use s = 4)\n",
            median(s_centred), min(s_centred), max(s_centred)))
sink()
cat(readLines("outputs/Suppl_parameter_economy.txt"), sep = "\n")
