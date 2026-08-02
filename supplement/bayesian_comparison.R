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
# Output: outputs/Suppl_bayesian_comparison.txt
# Usage:  Rscript supplement/bayesian_comparison.R
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
# rho_nei(s=1)/rho_vac per galaxy (outer-point aggregates); eps_H(s) follows Eq. 4
rho_ratio_s1 <- mcgaugh.mass / (4/3 * pi * (mcgaugh.R$x * kpc)^3) / rho_vac
epsH_of_s    <- function(s) sqrt(pmax(rho_ratio_s1 / s^3 + 1/6, 1e-6))

# Goodness per galaxy: (i) chi2 with the SPARC velocity uncertainty (Lelli,
# McGaugh & Schombert 2016: median 4.3% fractional error, 4.6 km/s floor);
# (ii) the error-model-free ln V residual for a profiled-variance BIC.
FRAC <- 0.043; FLOOR <- 4.6
gal_terms <- function(i, model, epsH = NA, gc = gamma_cen) {
  lg <- gal_name_row == gnames[i]
  Vk <- sqrt(mcgaugh$Vst[lg]^2 + mcgaugh$Vgas[lg]^2)
  R  <- mcgaugh$R[lg]; Vobs <- mcgaugh$Vobs[lg]
  aN <- Vk^2 / (R * kpc) * kms^2
  if (model == "MOND") {
    D <- 0.5 + sqrt(0.25 + a0 / aN)
  } else {
    vevH <- (sqrt(2) * Vk * kms * T0) / (R * kpc)
    quot <- abs(vevH^2 - epsH^2) / (epsH^2 + vevH^2)
    gsys <- asin(sqrt(pmin(pmax(sin(gamma_U)^2 + (sin(gc)^2 - sin(gamma_U)^2) * quot, 0), 1)))
    D    <- sqrt(1 + (2 * c0 / T0) / (aN * gsys / cos(gsys)))
  }
  Vmod <- Vk * sqrt(D); sig <- sqrt((FRAC * Vobs)^2 + FLOOR^2)
  ok <- is.finite(Vmod) & is.finite(Vobs) & Vobs > 0
  list(chi2 = sum(((Vmod[ok]-Vobs[ok])/sig[ok])^2),
       ssr  = sum((log(Vmod[ok])-log(Vobs[ok]))^2), n = sum(ok))
}
# total (chi2, ssr, N) over galaxies for a per-galaxy eps_H vector (or MOND)
gal_total <- function(epsH, model = "HMG") {
  C <- 0; S <- 0; n <- 0
  for (i in 1:Ngal) {
    if (model == "HMG" && !is.finite(epsH[i])) next
    a <- gal_terms(i, model, if (model == "HMG") epsH[i] else NA)
    if (a$n == 0) next
    C <- C + a$chi2; S <- S + a$ssr; n <- n + a$n
  }
  c(chi2 = C, ssr = S, n = n)
}

# HMG general model: a single GLOBAL density scale s (one fitted parameter),
# calibrated to all galaxies. The best-fit value confirms the s = 4 of the paper.
s_grid <- seq(1, 30, 0.02)
s_best <- s_grid[which.min(sapply(s_grid, function(s) gal_total(epsH_of_s(s))["chi2"]))]
gA <- gal_total(epsH_of_s(s_best))
N  <- as.numeric(gA["n"]); chi2A <- as.numeric(gA["chi2"]); ssrA <- as.numeric(gA["ssr"])
gM <- gal_total(NULL, "MOND"); chi2M <- as.numeric(gM["chi2"]); ssrM <- as.numeric(gM["ssr"])

# HMG with eps_H fitted freely per galaxy (most flexible variant)
grid_eps <- unique(c(seq(0.5, 10, 0.1), seq(10, 200, 0.5)))
chi2B <- 0; ssrB <- 0; epsB <- rep(NA, Ngal)
for (i in 1:Ngal) {
  if (!is.finite(rho_ratio_s1[i])) next
  if (gal_terms(i, "HMG", epsH_of_s(s_best)[i])$n == 0) next
  tb <- lapply(grid_eps, function(e) gal_terms(i, "HMG", e))
  chi2B <- chi2B + min(sapply(tb, `[[`, "chi2"))
  ssrB  <- ssrB  + min(sapply(tb, `[[`, "ssr"))
  epsB[i] <- grid_eps[which.min(sapply(tb, `[[`, "chi2"))]
}
used <- sum(!is.na(epsB))

# HMG with eps_H AND gamma_cen fitted freely per galaxy (2 params/galaxy: the
# analogue of the 248-parameter configuration penalised by Aydogdu 2026)
eg2 <- unique(c(seq(0.5, 10, 0.2), seq(10, 200, 1))); gg2 <- seq(0.44, 0.50, 0.005) * pi
chi2C <- 0; ssrC <- 0
for (i in 1:Ngal) {
  if (!is.finite(rho_ratio_s1[i]) || gal_terms(i, "HMG", eg2[1])$n == 0) next
  best <- Inf; bssr <- NA
  for (gc in gg2) for (e in eg2) { t <- gal_terms(i, "HMG", e, gc)
    if (t$chi2 < best) { best <- t$chi2; bssr <- t$ssr } }
  chi2C <- chi2C + best; ssrC <- ssrC + bssr
}

bicC <- function(chi2, k) chi2 + k * log(N)                 # BIC with SPARC errors
bicP <- function(ssr, k)  N * log(ssr / N) + k * log(N)     # BIC, profiled variance

gal_tab <- data.frame(
  model = c(sprintf("HMG global density scale s (best-fit %.2f)", s_best),
            "HMG eps_H fitted per galaxy",
            "HMG eps_H and gamma_cen per galaxy",
            "MOND (a0 fixed = 1.2e-10)"),
  fitted_params = c(1, used, 2*used, 0),
  k = c(1, used, 2*used, 0),
  chi2 = round(c(chi2A, chi2B, chi2C, chi2M)),
  chi2_nu = round(c(chi2A/(N-1), chi2B/(N-used), chi2C/(N-2*used), chi2M/N), 2),
  p_value = signif(c(pchisq(chi2A,N-1,lower.tail=FALSE), pchisq(chi2B,N-used,lower.tail=FALSE),
                     pchisq(chi2C,N-2*used,lower.tail=FALSE), pchisq(chi2M,N,lower.tail=FALSE)), 2),
  BIC = round(c(bicC(chi2A,1), bicC(chi2B,used), bicC(chi2C,2*used), bicC(chi2M,0))),
  BIC_prof = round(c(bicP(ssrA,1), bicP(ssrB,used), bicP(ssrC,2*used), bicP(ssrM,0))))
gal_tab$dBIC      <- round(gal_tab$BIC - gal_tab$BIC[4])
gal_tab$dBIC_prof <- round(gal_tab$BIC_prof - gal_tab$BIC_prof[4])
r_pred_fit <- cor(epsH_of_s(s_best), epsB, use = "complete.obs")

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
sink("outputs/Suppl_bayesian_comparison.txt")
cat("Bayesian comparison of HMG and MOND on galaxy rotation curves\n")
cat("eps_H^2 = rho_nei(s)/rho_vac + 1/6   (Monjo 2026, MNRAS 549, Eq. 4)\n")
cat(sprintf("HMG general model: one global density scale s (best-fit %.2f, k=1); MOND: a0 fixed (k=0).\n", s_best))
cat("61 galaxies with finite McGaugh parameters (3 single-point); main text quotes 60 (drops UGC 6923).\n\n")
cat(sprintf("GALAXIES  (McGaugh 2007: %d galaxies, %d rotation-curve points)\n", used, N))
cat(sprintf("chi2: sigma_V = sqrt((%.3f Vobs)^2 + %.1f^2) km/s [SPARC, Lelli+2016], BIC = chi2 + k lnN\n", FRAC, FLOOR))
cat("BIC_prof: profiled-variance BIC = N ln(SSR_lnV/N) + k lnN (no error model assumed)\n\n")
print(gal_tab, row.names = FALSE)
cat(sprintf("\nr(eps_H closed, eps_H fitted) = %.3f\n", r_pred_fit))
cat(sprintf("\nCLUSTERS  (%d clusters, %d RAR points): per-cluster scale s for the density relation\n\n",
            length(clusters), nrow(rar)))
print(clu_tab, row.names = FALSE)
cat(sprintf("\ns: median = %.2f, range %.2f - %.2f  (galaxies: global best-fit s = %.2f)\n",
            median(s_centred), min(s_centred), max(s_centred), s_best))
sink()
cat(readLines("outputs/Suppl_bayesian_comparison.txt"), sep = "\n")
