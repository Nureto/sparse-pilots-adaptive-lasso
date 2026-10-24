## =====================================================================
##  numerical_checks.R : numerical checks of Propositions 1, 3 and 4
##  (Appendix, "Numerical checks")
##  Usage:  Rscript R/numerical_checks.R     (run from PilotWeights_Q)
##  Needs:  R/common.R (for the weights as implemented) and, for check 2,
##          results/E1/E1_R200.rds. Expected run time 5 to 20 minutes.
##  Writes: results/numerical_checks.txt
##
##  Checks 4 and 5 use a plain coordinate descent solver written below, so
##  that they do not depend on glmnet.
## =====================================================================
source("R/common.R")                  # adaptive_weights(), cap_inert_predicted()
log_file <- "results/numerical_checks.txt"
sink(log_file, split = TRUE)
cat("Numerical checks of the propositions,", format(Sys.time()), "\n")
cat(R.version.string, "\n\n")

## ---------------------------------------------------------------------
##  Check 1. Proposition 1(a): the cap is inert exactly when
##  |Z| >= p - floor(0.95 (p - 1)) (type 7 quantiles), for every |Z|
## ---------------------------------------------------------------------
set.seed(20261021)
mism <- 0L; tried <- 0L
for (p in c(20L, 50L, 100L, 200L, 500L, 1000L)) {
  for (nz in 1:(p - 1L)) {
    b <- c(numeric(nz), runif(p - nz, 0.05, 3))
    act <- adaptive_weights(b, gamma = 2, delta = 1e-6, cap = 0.95)$cap_active
    tried <- tried + 1L
    if (act == cap_inert_predicted(nz, p)) mism <- mism + 1L
  }
}
cat(sprintf("Check 1, Prop 1(a): cap status disagrees with the condition in %d of %d configurations\n",
            mism, tried))

## ---------------------------------------------------------------------
##  Check 2. Proposition 1(c) on the 8,000 weight vectors of E1:
##  min_Z w / max_off w = 1 + min_off |b|^gamma / delta (log10 scale)
## ---------------------------------------------------------------------
e1f <- "results/E1/E1_R200.rds"
if (file.exists(e1f)) {
  e1 <- readRDS(e1f)
  ok <- !e1$cap_active & is.finite(e1$log10_ratio)
  pred <- log10(1 + e1$minabs_off^2 / e1$delta)
  cat(sprintf("Check 2, Prop 1(c) on E1: largest absolute error %.1e over %d weight vectors (cap active in %d)\n",
              max(abs(e1$log10_ratio[ok] - pred[ok])), sum(ok), sum(e1$cap_active)))
} else cat("Check 2 skipped: ", e1f, " not found\n")

## ---------------------------------------------------------------------
##  Check 3. Proposition 1(b) and (c) on 3,000 random configurations with
##  the cap inert: w_Z = p / Z_eff and w_k = kappa_k p / Z_eff
## ---------------------------------------------------------------------
set.seed(20261022)
worst_b <- 0; worst_c <- 0; capped <- 0L; N3 <- 3000L
for (t in seq_len(N3)) {
  p  <- sample(20:1499, 1L)
  z0 <- p - floor(0.95 * (p - 1))               # smallest |Z| with the cap inert
  nz <- if (z0 <= p - 1) sample(z0:(p - 1), 1L) else p - 1L
  dl <- 10^runif(1, -9, -1)
  off <- runif(p - nz, 1e-5, 3) * 10^runif(p - nz, -4, 0)
  b  <- c(numeric(nz), off)
  aw <- adaptive_weights(b, gamma = 2, delta = dl, cap = 0.95)
  if (aw$cap_active) { capped <- capped + 1L; next }
  w <- aw$w
  kappa <- dl / (off^2 + dl)
  zeff <- nz + sum(kappa)
  eb <- max(abs(w[seq_len(nz)] - p / zeff) / (p / zeff),
            abs(w[-seq_len(nz)] - kappa * p / zeff) / w[-seq_len(nz)])
  ec <- max(abs((w[1] / w[-seq_len(nz)]) - (1 + off^2 / dl)) / (1 + off^2 / dl))
  worst_b <- max(worst_b, eb); worst_c <- max(worst_c, ec)
}
cat(sprintf("Check 3, Prop 1(b),(c): largest relative error %.1e for (b) and %.1e for (c) over %d configurations (cap active in %d)\n",
            worst_b, worst_c, N3 - capped, capped))

## ---------------------------------------------------------------------
##  Coordinate descent for
##    (1/2n)||y - X b||^2 + lam * sum_j w_j alpha |b_j|
##                        + lam (1 - alpha)/2 * sum_j v_j b_j^2,
##  v = w (convention W) or v = 1 (convention L); columns of X have
##  mean 0 and mean square 1, y has mean 0.
## ---------------------------------------------------------------------
cd_en <- function(X, y, w, lam, a, conv = c("W", "L"), b0 = NULL,
                  tol = 1e-12, maxit = 20000L) {
  conv <- match.arg(conv)
  n <- nrow(X); p <- ncol(X)
  b <- if (is.null(b0)) numeric(p) else b0
  r <- as.vector(y - X %*% b)
  thr   <- lam * a * w
  ridge <- lam * (1 - a) * (if (conv == "W") w else rep(1, p))
  for (it in seq_len(maxit)) {
    dm <- 0
    for (j in seq_len(p)) {
      xj <- X[, j]
      z  <- sum(xj * r) / n + b[j]
      bn <- sign(z) * max(abs(z) - thr[j], 0) / (1 + ridge[j])
      if (bn != b[j]) { r <- r - xj * (bn - b[j]); dm <- max(dm, abs(bn - b[j])); b[j] <- bn }
    }
    if (dm < tol) break
  }
  attr(b, "iter") <- it
  b
}
std_cols <- function(X) { X <- sweep(X, 2, colMeans(X)); sweep(X, 2, sqrt(colMeans(X^2)), "/") }
unit_w   <- function(b, dl) {                     # weights as implemented
  aw <- adaptive_weights(b, gamma = 2, delta = dl, cap = 0.95)
  if (aw$cap_active) warning("cap active: Proposition 1(b) does not apply")
  aw$w
}

## ---------------------------------------------------------------------
##  Check 4. Proposition 3 at a fixed lambda: under both conventions an
##  inactive pilot zero satisfies |x_j' r| / n <= lam alpha p / Z_eff,
##  30 random problems, four values of delta, tolerance 1e-12
## ---------------------------------------------------------------------
set.seed(20261023)
N4 <- 30L; n <- 60L; p <- 120L; lam <- 0.15; a <- 0.2
viol <- c(W = 0L, L = 0L); solves <- 0L; maxiter <- 0L; Ssz <- integer(N4)
for (t in seq_len(N4)) {
  X <- matrix(rnorm(n * p), n, p) + 0.5 * rnorm(n)        # common factor
  X <- std_cols(X)
  beta <- numeric(p); beta[1:8] <- seq(3, 1.5, length.out = 8)
  y <- as.vector(X %*% beta) + rt(n, df = 3); y <- y - mean(y); y <- y / sqrt(mean(y^2))
  bp <- cd_en(X, y, rep(1, p), 0.08, 1, "W")              # sparse lasso pilot
  S <- which(bp != 0); Z <- which(bp == 0); Ssz[t] <- length(S)
  for (cv in c("W", "L")) {
    bprev <- NULL
    for (dl in c(1e-4, 1e-6, 1e-8, 1e-10)) {
      w <- unit_w(bp, dl)
      b <- cd_en(X, y, w, lam, a, cv, bprev); bprev <- b
      solves <- solves + 1L; maxiter <- max(maxiter, attr(b, "iter"))
      zeff <- length(Z) + sum(dl / (bp[S]^2 + dl))
      thr  <- lam * a * p / zeff                          # the threshold of Prop 3
      inact <- Z[b[Z] == 0]
      g <- abs(crossprod(X[, inact, drop = FALSE], y - X %*% b)) / n
      if (any(g > thr * (1 + 1e-6))) viol[cv] <- viol[cv] + 1L
    }
  }
}
cat(sprintf("Check 4, Prop 3: %d problems (mean |S| = %.1f), %d solves; threshold violated in %d (W) and %d (L); most sweeps used %d\n",
            N4, mean(Ssz), solves, viol["W"], viol["L"], maxiter))

## ---------------------------------------------------------------------
##  Check 5. Proposition 4 on glmnet's default path: 40 random problems,
##  n = 60, p = 100, 20 lambda values from lambda_max down to
##  1e-2 lambda_max with warm starts, convention W, alpha = 0.2
## ---------------------------------------------------------------------
set.seed(20261024)
N5 <- 40L; n <- 60L; p <- 100L; a <- 0.2; eps <- 1e-2; dl <- 1e-6
held <- 0L; entered_when_held <- 0L; entered_any <- 0L
for (t in seq_len(N5)) {
  X <- std_cols(matrix(rnorm(n * p), n, p))
  beta <- numeric(p); beta[1:8] <- runif(8, 1, 3)
  y <- as.vector(X %*% beta) + rnorm(n); y <- y - mean(y)
  bp <- cd_en(X, y, rep(1, p), 0.2, 1, "W", tol = 1e-11)  # sparse lasso pilot
  S <- which(bp != 0); Z <- which(bp == 0)
  w <- unit_w(bp, dl)
  g <- abs(crossprod(X, y)) / n
  lmax <- max(g / (a * w)); lmin <- eps * lmax
  lhs <- eps * max(g[S] * (1 + bp[S]^2 / dl))            # condition (5)
  rhs <- sqrt(sum(y^2)) / sqrt(n)
  b <- NULL; ent <- FALSE
  for (lm in exp(seq(log(lmax), log(lmin), length.out = 20L))) {
    b <- cd_en(X, y, w, lm, a, "W", b, tol = 1e-11)
    if (any(b[Z] != 0)) ent <- TRUE
  }
  entered_any <- entered_any + ent
  if (lhs > rhs) { held <- held + 1L; entered_when_held <- entered_when_held + ent }
}
cat(sprintf("Check 5, Prop 4: condition (5) held in %d of %d problems; a pilot zero entered at some point of the path in %d of those (in %d problems overall)\n",
            held, N5, entered_when_held, entered_any))
sink()
cat("\nSaved", log_file, "\n")
