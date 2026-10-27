# Supplementary Material

Code and data to reproduce results in "Near-boundary asymptotics for profile score statistics" by Matias Shedden and Karl Oskar Ekvall.

## Code

The root directory contains the R code. Run each script from the root directory.

Tables:

- `tab_cover.R`: Table 1.
- `zip_coverage_sims.R`: simulations for Tables 2 and 3 (uses `zip_coverage_engine.R` and `zip_skewness.R`).
- `zip_coverage_table.R`: Tables 2 and 3, and the numbers quoted in Section 4.1, from the simulation output. Run after `zip_coverage_sims.R`.

Figures:

- `fig_power.R`: Figure 1.
- `fig_zip_kde.R`: Figure 2.

## Results

The results are in `results/`, and the R code writes there.

## Requirements

Base `R` only.
