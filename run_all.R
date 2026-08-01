# Monjo & Banik (2025), ApJ 992, 35 — reproduce all paper figures from vendored data.
# Usage:  Rscript run_all.R      (run from anywhere; self-locates the repo root)
# Requires R packages: jpeg, png, fields, stringr

.args <- commandArgs(trailingOnly = FALSE)
.file <- sub("^--file=", "", .args[grep("^--file=", .args)])
if (length(.file) == 1) setwd(dirname(normalizePath(.file)))
dir.create("outputs", showWarnings = FALSE)

# Send stray diagnostic screen plots to a null device; paper figures use explicit pdf() devices.
pdf(NULL)

source("R/00_setup.R")
source("R/01_galaxy_parameters.R")
source("R/FigC2.R")
source("R/FigC3.R")
source("R/Fig1.R")
source("R/Fig2.R")

cat("\nDone. Figures and tables written to ./outputs/\n")
