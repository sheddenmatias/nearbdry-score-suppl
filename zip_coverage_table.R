## zip_coverage_table.R ---------------------------------------------------------
##
## Tables 2 and 3 for the zero-inflated Poisson example, from the pooled
## results of zip_coverage_sims.R. Writes the table to
## ZIP_RESDIR/tab_zip_coverage.tex, and the numbers quoted in the paper
## to ZIP_RESDIR/zip_coverage_macros.tex.
##
## Column headers:
##   Extended   C_n, extended-set maximizer (13), Theorem 8       (SCR_U)
##   Efficient  C_n^e, constrained maximizer (17), Theorem 2      (EFF_C)
##   Naive      C_n, constrained maximizer (17)                   (SCR_C, ONE_C)
##   Approx     coverage under R_n^*, Corollary 1 (two-sided, pred2;
##              lower bound, pred1)
##   LF         C_n^lf, least favorable critical value, Proposition 7 (LF)
##   Joint      box C^s_n with critical value m_alpha(hat-Omega_n),
##              Proposition 4 (BOX_M)
## a = tilde-beta_n(theta).
##
## Usage, from zip_code/:   Rscript zip_coverage_table.R

res_dir <- Sys.getenv("ZIP_RESDIR", "results")
dat <- readRDS(file.path(res_dir, "zip_coverage_dat.Rds"))

## One row per cell, statistics as columns.
ids <- c("psi", "lam", "n", "a", "h", "pred2", "pred1", "mbar")
wide <- reshape(dat[c(ids, "stat", "prop")], idvar = ids,
                timevar = "stat", direction = "wide")
names(wide) <- sub("^prop\\.", "", names(wide))
wide <- wide[order(wide$psi, wide$lam, wide$n), ]

f3 <- function(x) sprintf("%.3f", x)
f2 <- function(x) sprintf("%.2f", x)
fg <- function(x) format(x, trim = TRUE, scientific = FALSE)

rows <- character(0)
for (l in unique(wide$psi)) {
  w <- wide[wide$psi == l, ]
  rows <- c(rows, sprintf("\\multicolumn{11}{l}{$\\psi = %s$} \\\\", fg(l)))
  rows <- c(rows, sprintf(
    "%d & %s & $%s$ & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
    w$n, fg(w$lam), f2(w$a), f2(w$h),
    f3(w$SCR_U), f3(w$EFF_C), f3(w$SCR_C), f3(w$pred2),
    f3(w$ONE_C), f3(w$pred1), f3(w$LF)))
  if (l != tail(unique(wide$psi), 1)) rows <- c(rows, "\\addlinespace")
}

tab <- c(
  "\\begin{tabular}{rrrrrrrrrrr}",
  "\\toprule",
  " & & & & \\multicolumn{4}{c}{Two-sided} & \\multicolumn{3}{c}{Lower bound} \\\\",
  "\\cmidrule(lr){5-8} \\cmidrule(lr){9-11}",
  paste("$n$ & $\\lambda$ & $\\tilde\\beta_n(\\theta)$ & $h_n(\\theta)$ &",
        "Extended & Efficient & Naive & Approx & Naive & Approx & LF \\\\"),
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}")
writeLines(tab, file.path(res_dir, "tab_zip_coverage.tex"))

## Table 3: simultaneous coverage of boxes for both coordinates, from the
## same samples. The |kappa| column is |kappa_n(theta)| = |a| / (1 + a^2)^{1/2},
## also the absolute correlation of the inverse information.
rows <- character(0)
for (l in unique(wide$psi)) {
  w <- wide[wide$psi == l, ]
  rows <- c(rows, sprintf("\\multicolumn{8}{l}{$\\psi = %s$} \\\\", fg(l)))
  rows <- c(rows, sprintf(
    "%d & %s & %s & %s & %s & %s & %s & %s \\\\",
    w$n, fg(w$lam), f2(abs(w$a) / sqrt(1 + w$a^2)), f2(w$mbar),
    f3(w$SCR_LAM), f3(w$BOX_P), f3(w$BOX_B), f3(w$BOX_M)))
  if (l != tail(unique(wide$psi), 1)) rows <- c(rows, "\\addlinespace")
}
tab2 <- c(
  "\\begin{tabular}{rrrrrrrr}",
  "\\toprule",
  " & & & & \\multicolumn{1}{c}{Marginal} & \\multicolumn{3}{c}{Box} \\\\",
  "\\cmidrule(lr){5-5} \\cmidrule(lr){6-8}",
  paste("$n$ & $\\lambda$ & $|\\kappa|$ & $\\bar m$ & $\\lambda$ &",
        "Product & Bonferroni & Joint \\\\"),
  "\\midrule", rows, "\\bottomrule", "\\end{tabular}")
writeLines(tab2, file.path(res_dir, "tab_zip_box.tex"))

## Numbers quoted in Section 4.1 of the paper. Ranges are over all cells.
rng <- function(s) range(dat$prop[dat$stat == s])
nsim <- unique(dat$reps[dat$stat == "SCR_U"])
mac <- c(
  sprintf("\\newcommand{\\zipcovnsim}{%s}",
          format(max(nsim), big.mark = "{,}", scientific = FALSE)),
  sprintf("\\newcommand{\\zipvalidmaxdev}{%.3f}",
          max(abs(dat$prop[dat$stat %in% c("SCR_U", "EFF_C")] - 0.95))),
  sprintf("\\newcommand{\\zipnaivemin}{%.3f}", rng("SCR_C")[1]),
  sprintf("\\newcommand{\\zipnaiveonemin}{%.3f}", rng("ONE_C")[1]),
  sprintf("\\newcommand{\\zipwaldmin}{%.3f}", rng("WLD")[1]),
  sprintf("\\newcommand{\\zipwaldmax}{%.3f}", rng("WLD")[2]),
  sprintf("\\newcommand{\\ziplrtmin}{%.3f}", rng("LRT")[1]),
  sprintf("\\newcommand{\\ziplrtmax}{%.3f}", rng("LRT")[2]),
  sprintf("\\newcommand{\\ziplfmin}{%.3f}", rng("LF")[1]),
  sprintf("\\newcommand{\\zippreddiff}{%.3f}",
          max(abs(wide$SCR_C - wide$pred2))),
  sprintf("\\newcommand{\\zippredonediff}{%.3f}",
          max(abs(wide$ONE_C - wide$pred1))),
  sprintf("\\newcommand{\\zippimin}{%.3f}", rng("SCR_LAM")[1]),
  sprintf("\\newcommand{\\zippimax}{%.3f}", rng("SCR_LAM")[2]),
  sprintf("\\newcommand{\\zipboxpmin}{%.3f}", rng("BOX_P")[1]),
  sprintf("\\newcommand{\\zipboxbmin}{%.3f}", rng("BOX_B")[1]),
  sprintf("\\newcommand{\\zipboxbmax}{%.3f}", rng("BOX_B")[2]),
  sprintf("\\newcommand{\\zipboxmmin}{%.3f}", rng("BOX_M")[1]),
  sprintf("\\newcommand{\\zipboxmmax}{%.3f}", rng("BOX_M")[2]),
  sprintf("\\newcommand{\\zipmbarmin}{%.2f}", min(wide$mbar)),
  sprintf("\\newcommand{\\zipmbarmax}{%.2f}", max(wide$mbar)))
writeLines(mac, file.path(res_dir, "zip_coverage_macros.tex"))

cat(tab, sep = "\n"); cat("\n"); cat(mac, sep = "\n")
