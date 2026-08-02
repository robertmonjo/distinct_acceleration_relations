# Regression test for the reproducible pipeline.
#
# Regenerates all outputs and checks the per-cluster best-fit eps_H against the
# published Table 1 values, within a tolerance set by the eps_H grid resolution.
# Exits with status 0 if every check passes, 1 otherwise.
#
# Usage (from anywhere):  Rscript tests/test_outputs.R

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.file) == 1) setwd(dirname(dirname(normalizePath(.file))))  # repo root

# Run the pipeline in the same order as run_all.R (sourced directly so the
# working directory stays at the repository root).
dir.create("outputs", showWarnings = FALSE)
pdf(NULL)
for (script in c("R/00_setup.R", "R/01_galaxy_parameters.R", "R/Fig4.R",
                 "R/Fig5.R", "R/Fig1.R", "R/Fig2.R"))
  source(script)

pass <- TRUE
check <- function(label, got, expected, tol) {
  ok <- is.finite(got) && abs(got - expected) <= tol
  cat(sprintf("  [%s] %-22s got = %-6s  expected = %-4s (tol %s)\n",
              if (ok) "PASS" else "FAIL", label, got, expected, tol))
  if (!ok) pass <<- FALSE
}

# --- Table 1, one-parameter general model (gamma_cen = pi/2) ---
general <- read.table("outputs/Table1_general_model.txt", header = TRUE)
published_general <- c(A0085 = 40, A1795 = 55, A2029 = 52, A2142 = 46, A3158 = 53,
                       A0262 = 96, A2589 = 60, A3571 = 57, A0576 = 44, A0496 = 50)
for (name in names(published_general))
  check(paste0("general ", name), general[name, "eps_med"], published_general[[name]], 1.5)

# --- Table 1, two-parameter cluster model ---
cluster <- read.table("outputs/Table1_cluster_model.txt", header = TRUE)
published_cluster <- c(A0085 = 31, A1795 = 39, A2029 = 37, A2142 = 41, A3158 = 37,
                       A0262 = 67, A2589 = 41, A3571 = 39, A0576 = 31, A0496 = 33)
for (name in names(published_cluster))
  check(paste0("cluster ", name), cluster[name, "eps_med"], published_cluster[[name]], 2)

# --- All expected artefacts were produced ---
artefacts <- c("Fig1.pdf", "Fig2.pdf", "Fig4.pdf", "Fig5.pdf",
               "Table1_general_model.txt", "Table1_cluster_model.txt")
for (f in artefacts)
  check(paste0("exists ", f), as.integer(file.exists(file.path("outputs", f))), 1L, 0)

cat(if (pass) "\nALL REGRESSION CHECKS PASSED\n" else "\nREGRESSION TEST FAILED\n")
quit(status = if (pass) 0L else 1L)
