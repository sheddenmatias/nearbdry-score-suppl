## fig_zip_kde.R --------------------------------------------------------------
## Figure 2: four-panel figure for the zero-inflated Poisson example.
##
## DESIGN.  The two columns are two sample sizes, n = 25 and n = 200.  The
## parameter value (psi, lambda) differs between the columns and is chosen so
## that the two quantities that determine the distribution of R*_n,
##
##      kappa_n = kappa_n(theta)                  (score correlation)
##      h_n     = lambda * sqrt(n * I_ll(theta))  (standardized distance to boundary)
##
## are the same in both columns: kappa_n = -0.85 and h_n = 0.25.  With
## u_n = infinity, the approximating distribution is
##
##      R*_n = Z_psi + kappa (1 - kappa^2)^{-1/2} min(Z_lambda + h, 0),
##
## so it depends on the parameter and the sample size only through
## (kappa_n, h_n), and it is the same curve in both columns.  Anything that
## changes between the columns is therefore an effect of the sample size alone.
##
## The pairs (psi, lambda) are obtained by solving the 2x2 system
##      kappa_n(psi, lambda)       = KAPPA_TARGET
##      lambda sqrt(n I_ll)        = H_TARGET
## numerically; see solve_pair() below.
##
## ---------------------------------------------------------------------------

set.seed(405)
res_dir <- Sys.getenv("ZIP_RESDIR", "results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

KAPPA_TARGET <- -0.85 # score correlation kappa_n, common to both columns
H_TARGET     <- 0.25  # standardized distance to the boundary h_n, common to both
NSIM     <- 200000    # replicates per panel
NS       <- c(25, 200)

## ---- model quantities -------------------------------------------------------
## Per-observation information, with b = exp(-psi); entries are
## I_pp = I_{psi psi}, I_pl = I_{psi lambda}, I_ll = I_{lambda lambda}.
info_zip <- function(psi, lam) {
  b <- exp(-psi); p0 <- b + lam * (1 - b)
  c(I_pp = (1 - lam)^2 * b^2 / p0 + (1 - lam) / psi + lam - p0,
    I_pl = -b / p0,
    I_ll = (1 - b)^2 / (p0 * (1 - p0)))
}
## score correlation kappa_n(theta) = I_pl / (I_pp I_ll)^{1/2}
kappa_of <- function(psi, lam) {
  I <- info_zip(psi, lam)
  as.numeric(I["I_pl"] / sqrt(I["I_pp"] * I["I_ll"]))
}
h_of <- function(psi, lam, n) as.numeric(lam * sqrt(n * info_zip(psi, lam)["I_ll"]))

## Solve kappa_n(psi, lambda) = KAPPA_TARGET and h_n(psi, lambda) = H_TARGET for
## (psi, lambda).  Two equations, two unknowns; short Newton iteration.
solve_pair <- function(kappa_target, h_target, n, tol = 1e-12, maxit = 200) {
  x <- c(0.5, h_target / sqrt(n * 1.5))                 # (psi, lambda) start
  F <- function(x) c(kappa_of(x[1], x[2]) - kappa_target,
                     h_of(x[1], x[2], n) - h_target)
  for (it in seq_len(maxit)) {
    f <- F(x)
    if (max(abs(f)) < tol) break
    J <- matrix(0, 2, 2)                                 # numerical Jacobian
    for (j in 1:2) {
      e <- rep(0, 2); e[j] <- 1e-7 * max(1, abs(x[j]))
      J[, j] <- (F(x + e) - f) / e[j]
    }
    step <- solve(J, -f)
    while (any(x + step <= 0) || (x + step)[2] >= 1) step <- step / 2
    x <- x + step
  }
  x
}

pars <- lapply(NS, function(n) solve_pair(KAPPA_TARGET, H_TARGET, n))
names(pars) <- as.character(NS)
for (n in NS) {
  p <- pars[[as.character(n)]]
  cat(sprintf("n = %4d:  psi = %.4f  lambda = %.5f   kappa_n = %.4f   h_n = %.4f\n",
              n, p[1], p[2], kappa_of(p[1], p[2]), h_of(p[1], p[2], n)))
}

## ---- the approximating distribution --------------------------------------------------
## R*_n = Z_psi + beta min(Z_lambda + h, 0),  beta = kappa / (1 - kappa^2)^{1/2},
## with Z_psi, Z_lambda independent N(0,1).  Let V = -min(Z_lambda + h, 0) >= 0.
## P(V = 0) = Phi(h); on {V = m > 0} the density of V is phi(m + h) and R*_n is
## N(-beta m, 1).  Hence the density and distribution function below.
dRstar <- function(t, kappa, h, ngrid = 2000, mmax = 12) {
  beta <- kappa / sqrt(1 - kappa^2)
  m  <- seq(0, mmax, length.out = ngrid); w <- c(diff(m), 0)
  ph <- dnorm(m + h)
  pnorm(h) * dnorm(t) +
    vapply(t, function(x) sum(dnorm(x + beta * m) * ph * w), numeric(1))
}

## ---- simulation of the profile score statistic ------------------------------
## The statistic R-hat_n depends on the data only through n_0 = #{Y_i = 0} and
## sum(Y_i), and the partial maximizer is available in closed form:
##      lambda-hat(psi) = (n_0/n - b) / (1 - b),   b = exp(-psi),
## truncated to Lambda = [0, bar-lambda], bar-lambda = 0.999, for the
## constrained partial maximizer (17) (constrained = TRUE), and to the extended
## set tilde-Lambda(psi) = (-b/(1-b), 1) of (14) otherwise.  B = replicates.
sim_stat <- function(psi, lam, n, B, constrained, chunk = 5000) {
  b <- exp(-psi); out <- numeric(0)
  for (s in seq(1, B, by = chunk)) {
    m  <- min(chunk, B - s + 1)
    Y  <- matrix(rpois(m * n, psi), m, n)
    Y[matrix(runif(m * n) < lam, m, n)] <- 0        # zero inflation
    n0 <- rowSums(Y == 0)
    lh <- (n0 / n - b) / (1 - b)                       # hat-lambda_n(psi)
    lh <- if (constrained) pmin(pmax(lh, 0), 0.999)
          else             pmin(pmax(lh, -b / (1 - b) + 1e-9), 0.999)
    p0 <- b + lh * (1 - b)
    S_psi <- -n0 * (1 - lh) * b / p0 + (rowSums(Y) / psi - (n - n0))   # S_{n psi}, (A.3)
    iz <- t(vapply(lh, function(q) info_zip(psi, q), numeric(3)))
    out <- c(out, S_psi / sqrt(n * (iz[, 1] - iz[, 2]^2 / iz[, 3])))
  }
  out
}

## ---- plot -------------------------------------------------------------------
xg   <- seq(-4, 5.5, length.out = 600)
dstar <- dRstar(xg, KAPPA_TARGET, H_TARGET)         # same curve in every panel
dnrm  <- dnorm(xg)

col_kde  <- "#1F3B63"; fill_kde <- grDevices::adjustcolor(col_kde, 0.20)
col_star <- "#C0392B"; col_norm <- "grey35"

pdf(file.path(res_dir, "fig_zip_kde.pdf"), width = 7.4, height = 5.4, family = "Helvetica")
par(mfrow = c(2, 2), mar = c(3.3, 3.4, 2.4, 0.7), mgp = c(2.0, 0.55, 0),
    las = 1, cex.axis = 0.85, cex.lab = 0.95)

panels <- expand.grid(n = NS, constrained = c(TRUE, FALSE))
panels <- panels[order(!panels$constrained, panels$n), ]
labs <- c("(a)", "(b)", "(c)", "(d)")

for (i in seq_len(nrow(panels))) {
  n   <- panels$n[i]; con <- panels$constrained[i]
  p   <- pars[[as.character(n)]]
  Rn  <- sim_stat(p[1], p[2], n, NSIM, con)
  # floor the bandwidth at half the lattice spacing
  Ip    <- info_zip(p[1], p[2])
  eff   <- as.numeric(Ip["I_pp"] - Ip["I_pl"]^2 / Ip["I_ll"])   # efficient information
  delta <- 1 / (p[1] * sqrt(n * eff))     # lattice spacing of hat-R_n
  kd    <- density(Rn, bw = max(bw.nrd0(Rn), delta / 2), n = 512,
                   from = -4, to = 5.5)

  plot(NA, xlim = range(xg), ylim = c(0, 0.43), axes = FALSE, xlab = "", ylab = "")
  polygon(c(kd$x, rev(kd$x)), c(kd$y, rep(0, length(kd$y))), col = fill_kde,
          border = NA)
  lines(kd$x, kd$y, lwd = 2.6, col = col_kde)
  lines(xg, dnrm,  lwd = 1.8, lty = 2, col = col_norm)
  if (con) lines(xg, dstar, lwd = 2.2, col = col_star)
  axis(1, col = "grey40", col.axis = "grey20")
  axis(2, col = "grey40", col.axis = "grey20")
  title(xlab = expression(hat(R)[n]), ylab = "density", cex.lab = 1.05)
  title(main = sprintf("%s  %s,  n = %d", labs[i],
                       if (con) "constrained partial maximizer" else "extended-set maximizer", n),
        cex.main = 0.95, font.main = 1, adj = 0)

  if (i == 1)
    legend("topright", bty = "n", cex = 0.85, inset = c(0, 0.02),
           lwd = c(2.6, 2.2, 1.8), lty = c(1, 1, 2),
           col = c(col_kde, col_star, col_norm), text.col = "grey20",
           legend = expression(hat(R)[n], R[n]^symbol("*"), "N(0,1)"))
}
dev.off()
cat("wrote", file.path(res_dir, "fig_zip_kde.pdf"), "\n")
