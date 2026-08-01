# Distinct acceleration relations of galaxies and galaxy clusters

Reproduction code for the figures and tables of:

> **Monjo, R. & Banik, I. (2025).** *Distinct acceleration relations of galaxies and
> galaxy clusters from hyperconical modified gravity.* The Astrophysical Journal **992**, 35.
> DOI: [10.3847/1538-4357/adfcc0](https://doi.org/10.3847/1538-4357/adfcc0)

The pipeline is **self-contained**: all input data are vendored under [`data/`](data/) and every
path is relative. No network access or external files are required.

## Requirements

- R ≥ 4.0
- R packages: `jpeg`, `png`, `fields`, `stringr`

```r
install.packages(c("jpeg", "png", "fields", "stringr"))
```

## Quick start

From the repository root:

```sh
Rscript run_all.R
```

Figures and tables are written to [`outputs/`](outputs/). The script self-locates its own
directory, so it can be launched from anywhere.

## Repository layout

```
run_all.R                 # entry point: sources the ordered scripts below
R/
  00_setup.R              # packages, constants, models, data loading, preprocessing
  01_galaxy_parameters.R  # galaxy epsilon fits -> Fig00_galaxies_parameters.txt
  FigC2.R                 # Figure C2
  FigC3.R                 # Figure C3 + cluster tables
  Fig1.R                  # Figure 1
  Fig2.R                  # Figure 2
data/                     # vendored input datasets (see provenance below)
outputs/                  # generated figures and tables (git-ignored)
```

The scripts share a single R session (sourced in order by `run_all.R`); each builds on the
state of the previous one.

## Scripts, figures and tables

| Script | Paper output | File(s) in `outputs/` |
|--------|--------------|-----------------------|
| `Fig1.R`  | Figure 1  — cluster radial acceleration relation | `Fig01_cluster_RAR.pdf` |
| `Fig2.R`  | Figure 2  — galaxies vs. clusters | `Fig02_cluster_vs_galaxies.pdf` |
| `FigC2.R` | Figure C2 — fitting of the RAR | `FigC2_fitting_RAR.pdf` |
| `FigC3.R` | Figure C3 + tables | `FigC3_fitting_and_relations.pdf`, `table_cluster_RAR_B.txt`, `table_cluster_RAR_B2.txt`, `table_cluster_RAR_C.txt` |
| `01_galaxy_parameters.R` | intermediate galaxy fits | `Fig00_galaxies_parameters.txt` and `Fig00*` diagnostics |

## Data provenance

All files in `data/` are third-party datasets redistributed here for reproducibility. Please
cite the original sources when using them.

| File | Source |
|------|--------|
| `MassModels_Lelli2016c.mrt.txt` | SPARC mass models — Lelli, McGaugh & Schombert (2016), AJ 152, 157 ([SPARC database](http://astroweb.cwru.edu/SPARC/)) |
| `SPARC_BTFR_Lelli2019.txt` | SPARC baryonic Tully–Fisher relation — Lelli et al. (2019) |
| `McGaugh2007.txt`, `S_McGaugh2007_mass_discrepancy.txt` | Galaxy mass models — McGaugh (2007), ApJ 659, 149 |
| `RAR_data_Indranil.txt` | Binned galaxy radial acceleration relation (I. Banik) |
| `clusterRAR.dat`, `table_cluster_RAR_B.txt` | Galaxy-cluster radial acceleration relation compiled for this work |

## Citation

If you use this code, please cite the paper above.

## License

Code is released for reproducibility of the published results. Bundled datasets remain under
the terms of their original providers.
