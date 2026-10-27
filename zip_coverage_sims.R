## zip_coverage_sims.R ----------------------------------------------------------
##
## Coverage of nominal 95% regions for the Poisson mean in the zero-inflated
## Poisson example, with the zero-inflation probability a nuisance parameter
## at or near its boundary. One cell per (mean, zero-inflation, sample size).
##
## Usage, from zip_code/:
##     Rscript zip_coverage_sims.R
##
## Writes ZIP_RESDIR/zip_coverage_cells.Rds and
## ZIP_RESDIR/zip_coverage_dat.Rds, which zip_coverage_table.R reads to make
## Tables 2 and 3.

source("zip_coverage_engine.R")

N_SIM   <- as.integer(Sys.getenv("ZIP_NSIM", "100000"))
res_dir <- Sys.getenv("ZIP_RESDIR", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

## --- The cells ---------------------------------------------------------------
## Mean 0.5 gives tilde-beta_n near -1.8 (|kappa| near 0.88);
## and mean 1.25 gives tilde-beta_n near -1.0 (|kappa| near 0.71, the kappa = -0.70
## column of Table 1). The zero-inflation values range from the boundary to a
## standardized distance h_n above two at the largest sample size.
cells <- expand.grid(n = c(25L, 100L, 400L), lam = c(0, 0.02, 0.05, 0.2),
                     psi = c(0.5, 1.25), KEEP.OUT.ATTRS = FALSE)
cells <- cells[c("psi", "lam", "n")]
cells$seed <- 201000L + seq_len(nrow(cells))
cells$key  <- sprintf("psi%g_lam%g_n%d_B%d", cells$psi, cells$lam, cells$n, N_SIM)

## --- Run ---------------------------------------------------------------------
cache_file <- file.path(res_dir, "zip_coverage_cells.Rds")
cache <- if (file.exists(cache_file)) readRDS(cache_file) else list()

cat(sprintf("%d cells, %d replications each, results in %s/\n",
            nrow(cells), N_SIM, res_dir))
for (i in seq_len(nrow(cells))) {
  key <- cells$key[i]
  if (!is.null(cache[[key]])) { cat("  on disk:", key, "\n"); next }
  set.seed(cells$seed[i])
  cache[[key]] <- run_cell(cells$psi[i], cells$lam[i], cells$n[i], N_SIM)
  saveRDS(cache, cache_file)
}

## --- Pool ----------------------------------------------------------------------
dat <- do.call(rbind, cache[cells$key])
rownames(dat) <- NULL
saveRDS(dat, file.path(res_dir, "zip_coverage_dat.Rds"))
cat(sprintf("\nwrote %s: %d rows\n",
            file.path(res_dir, "zip_coverage_dat.Rds"), nrow(dat)))
