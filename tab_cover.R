## tab_cover.R -------------------------------------------------------------------
##
## Table 1: asymptotic lower and upper tail probabilities of the signed root
## with the constrained partial maximizer, from Corollary 1 with
## d_1 = d_2 = 1 and u_n = infinity, so that
##
##      R* = Z_psi + b min(Z_lam + h, 0),   b = kappa / (1 - kappa^2)^{1/2},
##
## with Z_psi and Z_lam (Z_lambda) independent standard normal; b is
## tilde-beta_n of (6). The table has kappa < 0,
## for which b < 0 and the shift b min(Z_lam + h, 0) is nonnegative.
## Conditionally on Z_lam = z the variable R* is normal with mean 
## b min(z + h, 0) and variance 1, which gives the integrals below. 
## Writes the table to ZIP_RESDIR/tab_cover.tex
##
## Usage, from zip_code/:   Rscript tab_cover.R

res_dir <- Sys.getenv("ZIP_RESDIR", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

Z_TWO  <- 1.96                               # two-sided nominal 0.95
KAPPAS <- c(-0.45, -0.70, -0.85, -0.95)      # score correlation
HS     <- c(0, 0.25, 0.5, 1, 1.5, 2, 3)      # standardized distance to the boundary

## P(R* > c) and P(R* < -c). On {Z_lam >= -h} the shift is zero, which gives
## the first term; the integral is over {Z_lam < -h}.
tail_upper <- function(kappa, h, c = Z_TWO) {
  b <- kappa / sqrt(1 - kappa^2)
  (1 - pnorm(c)) * pnorm(h) +
    integrate(function(z) (1 - pnorm(c - b * (z + h))) * dnorm(z), -Inf, -h,
              rel.tol = 1e-10)$value
}
tail_lower <- function(kappa, h, c = Z_TWO) {
  b <- kappa / sqrt(1 - kappa^2)
  pnorm(-c) * pnorm(h) +
    integrate(function(z) pnorm(-c - b * (z + h)) * dnorm(z), -Inf, -h,
              rel.tol = 1e-10)$value
}

lower <- outer(HS, KAPPAS, Vectorize(function(h, kappa) tail_lower(kappa, h)))
upper <- outer(HS, KAPPAS, Vectorize(function(h, kappa) tail_upper(kappa, h)))
size  <- lower[1, ] + upper[1, ]             # sum at h = 0, the largest over h

f3 <- function(x) sprintf("%.3f", x)
f2 <- function(x) sprintf("%.2f", x)

## Interleave the lower and upper columns, kappa by kappa.
cells <- matrix("", length(HS), 2 * length(KAPPAS))
cells[, seq(1, ncol(cells), 2)] <- f3(lower)
cells[, seq(2, ncol(cells), 2)] <- f3(upper)
rows <- sprintf("%s & %s \\\\", f2(HS), apply(cells, 1, paste, collapse = " & "))

head1 <- paste(" &", paste(sprintf("\\multicolumn{2}{c}{$\\kappa = %s$}", f2(KAPPAS)),
                           collapse = " & "), "\\\\")
rules <- paste(sprintf("\\cmidrule(lr){%d-%d}", 2 * seq_along(KAPPAS), 2 * seq_along(KAPPAS) + 1),
               collapse = " ")
head2 <- paste("$h$ &", paste(rep("Lower & Upper", length(KAPPAS)), collapse = " & "), "\\\\")
last  <- paste("Size &", paste(sprintf("\\multicolumn{2}{c}{%s}", f3(size)), collapse = " & "), "\\\\")

tab <- c(sprintf("\\begin{tabular}{l%s}", strrep("r", 2 * length(KAPPAS))),
         "\\toprule", head1, rules, head2, "\\midrule", rows,
         "\\midrule", last, "\\bottomrule", "\\end{tabular}")
writeLines(tab, file.path(res_dir, "tab_cover.tex"))
cat(tab, sep = "\n")
