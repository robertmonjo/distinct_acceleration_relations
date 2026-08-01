# Monjo & Banik (2025), ApJ 992, 35 — reproduce all paper figures from vendored data.
# Usage:  Rscript run_all.R      (run from anywhere; self-locates the repo root)
# Requires R packages: jpeg, png, fields, stringr

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.file) == 1) setwd(dirname(normalizePath(.file)))
dir.create("outputs", showWarnings = FALSE)

# Send stray diagnostic screen plots to a null device; paper figures use explicit pdf() devices.
pdf(NULL)

# Sourced in computation-dependency order (not figure order): Fig4/Fig5 (RAR fits) precede Fig1/Fig2 by data dependency.
source("R/00_setup.R")
source("R/01_galaxy_parameters.R")
source("R/Fig4.R")
source("R/Fig5.R")
source("R/Fig1.R")
source("R/Fig2.R")

cat("\nDone. Figures and tables written to ./outputs/\n")
