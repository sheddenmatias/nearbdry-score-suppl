## ---------------------------------------------------------------------------
##  zip_skewness.R
##
##  Model quantities for the zero-inflated Poisson example (Section 4.1),
##  sourced by zip_coverage_engine.R.
## ---------------------------------------------------------------------------


## --- Model quantities ------------------------------------------------------

## Per-observation expected information for the zero-inflated Poisson law
##   P(Y = 0) = lam + (1 - lam) e^{-psi},  P(Y = y) = (1 - lam) e^{-psi} psi^y / y!
## Returns the three distinct entries of I(psi, lam), (A.4) in the paper:
##   i_pp = I_{psi psi},  i_pl = I_{psi lambda},  i_ll = I_{lambda lambda}.
## At lam = 0 this reduces to
##   i_pp = 1/psi,  i_pl = -1,  i_ll = e^{psi} - 1.
zip_info <- function(psi, lam) {
  el <- exp(-psi)                  # b in the paper
  cc <- 1 - el
  p0 <- lam + (1 - lam) * el       # p_0 in the paper
  list(
    i_pp = (1 - lam)^2 * el^2 / p0 + (1 - lam) / psi + lam - p0,
    i_pl = -el * ((1 - lam) * cc / p0 + 1),
    i_ll = cc^2 / p0 + cc / (1 - lam)
  )
}

## a = tilde-beta_n(theta) of (6), = kappa_n / (1 - kappa_n^2)^{1/2}.  The
## factors of n cancel, so this does not depend on the sample size.
model_constants <- function(psi, lam) {
  i <- zip_info(psi, lam)
  eff <- i$i_pp - i$i_pl^2 / i$i_ll   # per-observation efficient information bar-I
  list(
    a    = i$i_pl / sqrt(eff * i$i_ll),
    i_ll = i$i_ll
  )
}

## Standardized distance from the true nuisance value to the boundary,
##   h_n(theta) = lambda {n I_{lambda lambda}(theta)}^{1/2}.
h_of_n <- function(n, lam, i_ll) lam * sqrt(n * i_ll)


## --- Distribution function of the approximating variable -------------------

## T* = R*_n = Z_psi + a min(Z_lambda + h, 0) of Section 3 (R*_{b,h} with b = a).
## P(T* <= t) = Phi(h) Phi(t) + int_0^inf Phi(t + a m) phi(m + h) dm.
tstar_cdf <- function(t, a, h) {
  integrand <- function(m) pnorm(t + a * m) * dnorm(m + h)
  pnorm(h) * pnorm(t) +
    integrate(integrand, 0, Inf, rel.tol = 1e-12, subdivisions = 2000L)$value
}
