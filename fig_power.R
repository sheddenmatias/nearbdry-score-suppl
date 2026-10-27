## Figure 1. Power of the upper one-sided test in the approximating experiment.
## Constrained profile score with least favorable critical value versus the
## efficient score with the normal critical value.  a < 0 throughout.
## Notation: a = b = kappa (1 - kappa^2)^{-1/2} of Section 3; Fstar(t, a, h) is
## the distribution function of R*_{b,h}; lf_cv(a) = c_alpha(kappa);
## pow_constr = P{R*_{n,delta} > c_alpha(kappa)}; pow_eff = 1 - Phi(z - delta).
Fstar <- function(t, a, h)
  pnorm(h) * pnorm(t) + integrate(function(m) pnorm(t + a * m) * dnorm(m + h),
                                  0, Inf, rel.tol = 1e-10, subdivisions = 2000L)$value
lf_cv <- function(a, alpha = 0.05)
  uniroot(function(c) 1 - Fstar(c, a, 0) - alpha, c(0, 20), tol = 1e-10)$root
pow_constr <- function(delta, a, h, cv) 1 - Fstar(cv - delta, a, h + a * delta)
pow_eff <- function(delta, alpha = 0.05) 1 - pnorm(qnorm(1 - alpha) - delta)

kap <- c(0.70, 0.85, 0.95)              # the |kappa| grid of Table 1, less 0.45
as <- -kap / sqrt(1 - kap^2); cvs <- sapply(as, lf_cv)
cat("kappa:", sprintf("%.2f", kap), " cv:", sprintf("%.2f", cvs), "\n")
hs <- seq(0, 3, by = 0.02)

res_dir <- Sys.getenv("ZIP_RESDIR", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)
pdf(file.path(res_dir, "fig_power.pdf"), width = 6.5, height = 3.1)
par(mfrow = c(1, 2), mar = c(4, 4, 2, 0.6), mgp = c(2.4, 0.7, 0), cex = 0.85)
ltys <- c(3, 2, 1)                       # kappa = 0.70 dotted, 0.85 dashed, 0.95 solid
for (k in 1:2) {
  delta <- c(1, 2)[k]
  plot(NA, xlim = c(0, 3), ylim = c(0, 1), xlab = "Standardized distance to the boundary",
       ylab = "Power", las = 1, xaxs = "i", yaxs = "i", bty = "l")
  abline(h = pow_eff(delta), col = "grey55", lwd = 1.5)
  for (i in seq_along(as)) {
    p <- sapply(hs, pow_constr, delta = delta, a = as[i], cv = cvs[i])
    lines(hs, p, lty = ltys[i], lwd = 1.5)
    ## direct labels: left end in panel (a); where the curves separate in (b)
    if (k == 1) {
      text(0.07, p[1], sprintf("%.2f", kap[i]), pos = 4, cex = 0.8, offset = 0.2)
    } else {
      xy <- list(c(2.05, 0.45, 0, 1), c(1.55, 0.84, 1, 1), c(1.5, 0.97, 0, 1))[[i]]  # kappa 0.70, 0.85, 0.95
      text(xy[1], xy[2], sprintf("%.2f", kap[i]), adj = xy[3:4], cex = 0.8)
    }
  }
  text(2.93, pow_eff(delta) + 0.025, "efficient score", adj = c(1, 0), cex = 0.8, col = "grey30")
  mtext(if (k == 1) expression("(a)  " * delta == 1) else expression("(b)  " * delta == 2),
        side = 3, line = 0.4, adj = 0, cex = 0.85)
}
dev.off()

## Checks for numbers quoted in the paper after Figure 1 (not in the figure).
## lower-tail loss at h = 0 (reject for small values, cv = -1.645), delta < 0
cat("\nlower one-sided, delta < 0, h = 0: constrained vs efficient\n")
for (i in seq_along(as)) for (d in c(-0.5, -1, -2))
  cat(sprintf(" kappa=%.2f delta=%.1f: %.3f vs %.3f\n", kap[i], d,
              Fstar(-1.645 - d, as[i], 0 + as[i] * d), pnorm(-1.645 - d)))

## Two-sided symmetric region |R| > c, with c set by the hypothesis at h = 0,
## which is where the two tails sum to the most (Table 1). Its power against
## delta = -1, where no boundary effect helps, versus the efficient score test.
cat(sprintf("\ntwo-sided symmetric, delta = -1: constrained vs efficient %.3f\n",
            1 - pnorm(1.96 + 1) + pnorm(-1.96 + 1)))
for (i in seq_along(as)) {
  cs <- uniroot(function(c) 1 - Fstar(c, as[i], 0) + Fstar(-c, as[i], 0) - 0.05,
                c(1, 20), tol = 1e-10)$root
  for (h in c(0, 1, 2))
    cat(sprintf(" kappa=%.2f h=%.0f: cv %.2f, power %.3f\n", kap[i], h, cs,
                1 - Fstar(cs + 1, as[i], h - as[i]) + Fstar(-cs + 1, as[i], h - as[i])))
}
