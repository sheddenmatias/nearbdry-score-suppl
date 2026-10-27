## zip_coverage_engine.R -------------------------------------------------------
##
## The statistics and one cell of the coverage simulation for the zero-inflated
## Poisson example (Tables 2 and 3). Sourced by zip_coverage_sims.R.
## The model quantities come from zip_skewness.R
##
## b = exp(-psi) and p0 are as in Appendix A.3.

source("zip_skewness.R")   # zip_info(), model_constants(), h_of_n(), tstar_cdf()

Z_TWO  <- qnorm(0.975)
Z_ONE  <- qnorm(0.95)
Q_CHI  <- qchisq(0.95, 1)
LAM_MAX <- 0.999           # bar-lambda: Lambda(psi) = [0, bar-lambda] in (17)
A_MIN  <- -12              # below this the least favorable critical value is frozen
Z_BONF <- qnorm(1 - 0.05 / 4)   # Bonferroni critical value for two coordinates

## Least favorable critical value
## c_alpha(kappa) of Section 3, as a function of a = kappa (1 - kappa^2)^{-1/2}:
## the larger of z and the 0.95-quantile of R*_{a,0}. For a >= 0 the
## approximating variable is at most Z_psi, so c(a) = z. For a < 0, the
## quantile is tabulated on a grid and interpolated.
lf_critical <- local({
  grid <- seq(A_MIN, 0, by = 0.05)
  q <- vapply(grid, function(a)
    uniroot(function(t) tstar_cdf(t, a, 0) - 0.95, c(Z_ONE, 80), tol = 1e-9)$root,
    numeric(1))
  f <- splinefun(grid, q, method = "natural")
  function(a) ifelse(a >= 0, Z_ONE, pmax(Z_ONE, f(pmax(a, A_MIN))))
})

## Maximum-modulus critical value for two coordinates
## m(rho) = m_alpha(Omega) of Section 2.4, the 0.95-quantile of
## max(|Z_1|, |Z_2|) for a standard bivariate normal pair with correlation
## rho = the off-diagonal of Omega. Tabulated on a grid and interpolated.
maxmod_quantile <- local({
  F <- function(c, rho) {
    s <- sqrt(1 - rho^2)
    integrate(function(z) dnorm(z) * (pnorm((c - rho * z) / s) - pnorm((-c - rho * z) / s)),
              -c, c, rel.tol = 1e-10)$value
  }
  grid <- c(seq(0, 0.99, by = 0.01), 0.995, 0.999)
  q <- vapply(grid, function(r)
    uniroot(function(c) F(c, r) - 0.95, c(Z_TWO, Z_BONF + 1e-6), tol = 1e-9)$root,
    numeric(1))
  f <- splinefun(c(grid, 1), c(q, Z_TWO), method = "natural")
  function(rho) pmin(pmax(f(pmin(abs(rho), 1)), Z_TWO), Z_BONF)
})

## Root of the score for the mean at fixed zero-inflation
## The nuisance estimator of psi, a root over Psi(lambda), when lambda is the
## interest parameter (Table 3). Its set (0, Inf) is open, so no extension is needed.
psi_root <- function(lam, n0, sy, npos, iter = 60L) {
  ok <- npos > 0
  lo <- rep(1e-10, length(sy)); hi <- 2 * ifelse(ok, sy / pmax(npos, 1), 1) + 1
  g <- function(s) { b <- exp(-s); -n0 * (1 - lam) * b / (lam + (1 - lam) * b) + sy / s - npos }
  for (k in seq_len(iter)) {
    mid <- (lo + hi) / 2
    pos <- g(mid) > 0
    lo[pos] <- mid[pos]; hi[!pos] <- mid[!pos]
  }
  out <- (lo + hi) / 2; out[!ok] <- NA; out
}

## Zero-truncated Poisson root
## Solves psi / (1 - exp(-psi)) = ybar_pos for psi > 0 by bisection.
## The left side increases from 1 to infinity, so a root exists exactly when
## ybar_pos > 1, NA otherwise.
ztp_root <- function(ybar_pos, iter = 60L) {
  ok <- !is.na(ybar_pos) & ybar_pos > 1 + 1e-10
  y  <- ifelse(ok, ybar_pos, 2)
  lo <- rep(0, length(y)); hi <- y
  for (k in seq_len(iter)) {
    mid <- (lo + hi) / 2
    up  <- mid / (1 - exp(-mid)) > y
    hi[up] <- mid[up]; lo[!up] <- mid[!up]
  }
  out <- (lo + hi) / 2; out[!ok] <- NA; out
}

## One chunk of replications
## Returns a 0/1 matrix, one row per region and one column per replication,
## with 1 where the region contains the true mean, plus a row with the box
## critical value. NA where a statistic is undefined.
chunk_stats <- function(psi, lam, n, m) {
  Y <- matrix(rpois(m * n, psi), m, n)
  Y[matrix(runif(m * n) < lam, m, n)] <- 0
  n0 <- rowSums(Y == 0); sy <- rowSums(Y); npos <- n - n0
  b  <- exp(-psi); p0hat <- n0 / n            # p0hat = hat-p_{0,n}

  ## Profile maximizers at the hypothesized mean, in closed form:
  ## lh_u over the extended set tilde-Lambda(psi) of (13)-(14), lh_c over
  ## Lambda(psi) = [0, LAM_MAX] of (17). Off the event E_n = {0 < n0 < n} the
  ## first has no maximizer and is placed just inside tilde-Lambda(psi).
  lh_u <- pmin(pmax((p0hat - b) / (1 - b), -b / (1 - b) + 1e-9), LAM_MAX)
  lh_c <- pmin(pmax(lh_u, 0), LAM_MAX)

  ## S1 = S_{n psi}, S2 = S_{n lambda}, from (A.3)
  score <- function(lh) {
    p0 <- lh + (1 - lh) * b
    list(S1 = -n0 * (1 - lh) * b / p0 + sy / psi - npos,
         S2 =  n0 * (1 - b) / p0 - npos / (1 - lh))
  }
  ## eff = bar-I (per observation), K = beta_n = I_{psi lambda} / I_{lambda lambda}
  einfo <- function(s, l) {
    i <- zip_info(s, l)
    list(eff = i$i_pp - i$i_pl^2 / i$i_ll, K = i$i_pl / i$i_ll)
  }
  su <- score(lh_u); iu <- einfo(psi, lh_u)
  sc <- score(lh_c); ic <- einfo(psi, lh_c)

  scr_u <- su$S1^2 / (n * iu$eff)                       # extended set, Theorem 8
  r_c   <- sc$S1 / sqrt(n * ic$eff)                     # signed root hat-R_n, constrained
  scr_c <- r_c^2                                        # naive
  eff_c <- (sc$S1 - ic$K * sc$S2)^2 / (n * ic$eff)      # (hat-R^e_n)^2, Theorem 2
  a_hat <- model_constants(psi, lh_c)$a                 # tilde-beta_n(hat-theta_n)

  ## Constrained maximum likelihood estimator of (psi, lam), lam in [0, LAM_MAX].
  ## If there are no positive observations the Wald and likelihood ratio are NA.
  ybar_pos <- ifelse(npos > 0, sy / pmax(npos, 1), NA)
  psi_zt   <- ztp_root(ybar_pos)
  lam_zt   <- (p0hat - exp(-psi_zt)) / (1 - exp(-psi_zt))
  interior <- !is.na(lam_zt) & lam_zt >= 0
  psi_hat  <- ifelse(interior, psi_zt, sy / n)
  lam_hat  <- ifelse(interior, pmin(lam_zt, LAM_MAX), 0)

  ## Wald: squared error times the expected efficient information at the estimate
  wald <- (psi_hat - psi)^2 * n * einfo(psi_hat, lam_hat)$eff

  ## Likelihood ratio using constrained profile log-likelihood of the mean
  lp <- function(s) {
    e <- exp(-s); p0 <- pmax(p0hat, e)
    n0 * log(p0) + ifelse(npos > 0, npos * (log1p(-p0) - log1p(-e)), 0) +
      sy * log(s) - npos * s
  }
  lrt <- pmax(0, 2 * (lp(psi_hat) - lp(psi)))
  wald[n0 == n] <- NA; lrt[n0 == n] <- NA

  ## Box for both coordinates
  psi_t <- psi_root(lam, n0, sy, npos)
  bt <- exp(-psi_t); it <- zip_info(psi_t, lam)
  r_lam <- (n0 * (1 - bt) / (lam + (1 - lam) * bt) - npos / (1 - lam)) /   # hat-R_{n,2}(lambda)
             sqrt(n * (it$i_ll - it$i_pl^2 / it$i_pp))
  r_psi <- su$S1 / sqrt(n * iu$eff)                                       # hat-R_{n,1}(psi)
  ih <- zip_info(psi_hat, lam_hat)
  m_hat <- maxmod_quantile(abs(ih$i_pl) / sqrt(ih$i_pp * ih$i_ll))        # m_alpha(hat-Omega_n)
  mx <- pmax(abs(r_psi), abs(r_lam))

  rbind(SCR_U = scr_u <= Q_CHI,                # two-sided, nominal 0.95
        EFF_C = eff_c <= Q_CHI,
        SCR_C = scr_c <= Q_CHI,
        WLD   = wald  <= Q_CHI,
        LRT   = lrt   <= Q_CHI,
        ONE_C = r_c <= Z_ONE,                  # lower bound, nominal 0.95
        LF    = r_c <= lf_critical(a_hat),     # C^lf_n, Proposition 7
        SCR_LAM = r_lam^2 <= Q_CHI,            # marginal interval for lambda, extended-set mean
        BOX_P  = mx <= Z_TWO,                  # product of the two marginal intervals
        BOX_B  = mx <= Z_BONF,                 # Bonferroni box
        BOX_M  = mx <= m_hat,                  # box C^s_n of Proposition 4
        MM     = m_hat) * 1                    # its critical value, averaged in run_cell
}

## --- One cell ------------------------------------------------------------------
## n_sim replications in chunks. Returns one row per region with the coverage
## proportion and the number of replications that produced a value, together
## with the cell's parameters, a = tilde-beta_n(theta), h_n(theta), and the
## asymptotic coverages of the naive regions.
run_cell <- function(psi, lam, n, n_sim, chunk = 10000L) {
  t0 <- proc.time()[["elapsed"]]
  hits <- NULL; reps <- NULL
  for (s in seq(1L, n_sim, by = chunk)) {
    M <- chunk_stats(psi, lam, n, min(chunk, n_sim - s + 1L))
    h <- rowSums(M, na.rm = TRUE); r <- rowSums(!is.na(M))
    hits <- if (is.null(hits)) h else hits + h
    reps <- if (is.null(reps)) r else reps + r
  }
  prop <- hits / reps
  stats  <- c("SCR_U", "EFF_C", "SCR_C", "WLD", "LRT", "ONE_C", "LF",
              "SCR_LAM", "BOX_P", "BOX_B", "BOX_M")

  k <- model_constants(psi, lam)
  h_n <- h_of_n(n, lam, k$i_ll)
  cat(sprintf("[%s] psi = %g  lam = %g  n = %d:  a = %.3f  h = %.2f  %ds\n",
              format(Sys.time(), "%H:%M:%S"), psi, lam, n, k$a, h_n,
              round(proc.time()[["elapsed"]] - t0)))

  data.frame(psi = psi, lam = lam, n = n, a = k$a, h = h_n, mbar = prop[["MM"]],
             pred2 = tstar_cdf(Z_TWO, k$a, h_n) - tstar_cdf(-Z_TWO, k$a, h_n),
             pred1 = tstar_cdf(Z_ONE, k$a, h_n),
             stat = stats, prop = prop[stats],
             reps = reps[stats], row.names = NULL)
}
