## =====================================================================
##  common.R : shared code for the pilot-weights study
##  Sourced by E1 to E4. Requires: glmnet (and hdi for E4).
## =====================================================================
suppressPackageStartupMessages(library(glmnet))

## ---------------------------------------------------------------------
##  Data generation
## ---------------------------------------------------------------------
gen_x <- function(n, p, corr = c("equi", "ar1", "block"),
                  rho = 0.5, block = 10L, rho_block = 0.7) {
  corr <- match.arg(corr)
  Z <- matrix(rnorm(n * p), n, p)
  if (corr == "equi") {
    f <- rnorm(n)
    X <- sqrt(rho) * f + sqrt(1 - rho) * Z          # f recycles down columns
  } else if (corr == "ar1") {
    X <- Z
    if (p > 1L) for (j in 2:p)
      X[, j] <- rho * X[, j - 1L] + sqrt(1 - rho^2) * Z[, j]
  } else {
    X <- Z
    nb <- ceiling(p / block)
    for (b in seq_len(nb)) {
      idx <- ((b - 1L) * block + 1L):min(b * block, p)
      f <- rnorm(n)
      X[, idx] <- sqrt(rho_block) * f + sqrt(1 - rho_block) * Z[, idx]
    }
  }
  X
}

## s active coefficients, spread evenly, decreasing linearly from 3 to 1.5
true_beta <- function(p, s = 10L) {
  b <- numeric(p)
  idx <- unique(round(seq(1, p, length.out = s)))
  b[idx] <- seq(3, 1.5, length.out = length(idx))
  b
}

## scen: "normal", "t3", or "t3_y10" (t3 errors, 10% of responses x 10)
gen_data <- function(n, p, corr, scen, s = 10L, beta = NULL) {
  X <- gen_x(n, p, corr)
  if (is.null(beta)) beta <- true_beta(p, s)
  e <- if (scen == "normal") rnorm(n) else rt(n, df = 3)
  y <- as.vector(X %*% beta) + e
  bad <- integer(0)
  if (scen == "t3_y10") {
    bad <- sample.int(n, max(1L, round(0.10 * n)))
    y[bad] <- 10 * y[bad]
  }
  list(x = X, y = y, beta = beta, bad = bad)
}

## ---------------------------------------------------------------------
##  Pilots
## ---------------------------------------------------------------------
safe_coef <- function(fit_obj, s, p) {
  if (is.null(fit_obj)) return(rep(0, p))
  out <- tryCatch(as.vector(stats::coef(fit_obj, s = s))[-1],
                  error = function(e) rep(NA_real_, p))
  if (length(out) != p || any(!is.finite(out))) rep(0, p) else out
}

pilot_lasso <- function(x, y) {
  f <- tryCatch(cv.glmnet(x, y, alpha = 1, nfolds = 5L), error = function(e) NULL)
  safe_coef(f, "lambda.min", ncol(x))
}

pilot_ridge <- function(x, y) {
  f <- tryCatch(cv.glmnet(x, y, alpha = 0, nfolds = 5L), error = function(e) NULL)
  safe_coef(f, "lambda.min", ncol(x))
}

## Huber-IRLS elastic-net proxy, copied from raenrmm::huber_pilot.
## alpha = 0.2 gives the sparse pilot used by the estimator; alpha = 0 gives a
## dense robust (ridge) pilot.
pilot_huber <- function(x, y, alpha = 0.2, passes = 3L, cc = 1.345) {
  p <- ncol(x)
  f0 <- tryCatch(cv.glmnet(x, y, alpha = alpha, nfolds = 3L, nlambda = 20L),
                 error = function(e) NULL)
  if (is.null(f0)) return(rep(0, p))
  b <- safe_coef(f0, "lambda.1se", p)
  lam <- f0$lambda.1se
  for (h in seq_len(passes)) {
    r <- as.vector(y - x %*% b)
    s <- stats::mad(r)
    if (!is.finite(s) || s < 1e-10) break
    u <- r / s
    rw <- pmin(1, cc / abs(u)); rw[!is.finite(rw)] <- 1
    srw <- sqrt(rw)
    fg <- tryCatch(glmnet(x * srw, y * srw, alpha = alpha, lambda = lam,
                          standardize = FALSE), error = function(e) NULL)
    if (is.null(fg)) break
    bn <- as.vector(stats::coef(fg))[-1]
    if (length(bn) != p || any(!is.finite(bn))) break
    if (sqrt(sum((bn - b)^2)) < 1e-3) { b <- bn; break }
    b <- bn
  }
  b
}

run_pilot <- function(name, x, y) {
  switch(name,
         lasso       = pilot_lasso(x, y),
         ridge       = pilot_ridge(x, y),
         huber       = pilot_huber(x, y, alpha = 0.2),
         huber_dense = pilot_huber(x, y, alpha = 0),
         stop("unknown pilot: ", name))
}

## ---------------------------------------------------------------------
##  Adaptive weights as implemented: raw -> optional 95% cap -> unit mean
## ---------------------------------------------------------------------
adaptive_weights <- function(b, gamma = 2, delta = 1e-6, cap = 0.95) {
  w <- 1 / (abs(b)^gamma + delta)
  cap_active <- FALSE
  if (!is.null(cap)) {
    q <- unname(stats::quantile(w, cap))          # type 7, R default
    cap_active <- any(w > q)
    w <- pmin(w, q)
  }
  list(w = w / mean(w), cap_active = cap_active)
}

## Corrected Proposition 1(a): with type-7 quantiles the cap is inert
## exactly when |Z| >= p - floor(0.95 (p - 1)).
cap_inert_predicted <- function(nz, p, cap = 0.95)
  nz >= p - floor(cap * (p - 1))

weight_diag <- function(b, w, cap_active, cap = 0.95) {
  p <- length(b)
  Z <- which(b == 0)
  off <- setdiff(seq_len(p), Z)
  nz <- length(Z)
  data.frame(
    Zsize        = nz,
    Zfrac        = nz / p,
    cap_active   = cap_active,
    cap_inert_pred = if (is.null(cap)) NA else cap_inert_predicted(nz, p, cap),
    maxwZ        = if (nz) max(w[Z]) else NA_real_,
    bound        = if (nz) p / nz else NA_real_,
    minwZ        = if (nz) min(w[Z]) else NA_real_,
    maxw_off     = if (length(off)) max(w[off]) else NA_real_,
    minw_off     = if (length(off)) min(w[off]) else NA_real_,
    log10_ratio  = if (nz && length(off)) log10(min(w[Z]) / max(w[off])) else NA_real_,
    minabs_off   = if (length(off)) min(abs(b[off])) else NA_real_
  )
}

## ---------------------------------------------------------------------
##  Fitting (squared loss, glmnet). Returns c(intercept, beta).
## ---------------------------------------------------------------------
lambda_fixed_nonadaptive <- function(x, y, alpha) {
  f <- cv.glmnet(x, y, alpha = alpha, nfolds = 3L, nlambda = 20L)
  f$lambda.1se
}

## glmnet needs >= 2 columns; pad with a constant column (coefficient 0).
.pad <- function(xk, pf) {
  if (ncol(xk) >= 2L) return(list(x = xk, pf = pf, padded = FALSE))
  list(x = cbind(xk, 0), pf = c(pf, 1), padded = TRUE)
}

## glmnet convention: penalty.factor multiplies the WHOLE elastic-net penalty,
##   lambda * sum_j w_j [ alpha |b_j| + (1 - alpha)/2 b_j^2 ].
fit_adaptive <- function(x, y, alpha, w, keep = NULL,
                         lambda_mode = c("fixed", "retuned"), lam_fixed = NULL) {
  lambda_mode <- match.arg(lambda_mode)
  p <- ncol(x)
  if (is.null(keep)) keep <- seq_len(p)
  out <- numeric(p + 1L)
  if (!length(keep)) { out[1] <- mean(y); return(out) }
  pd <- .pad(x[, keep, drop = FALSE], w[keep])
  lam_used <- lam_fixed; at_edge <- NA
  if (lambda_mode == "fixed") {
    f <- glmnet(pd$x, y, alpha = alpha, lambda = lam_fixed, penalty.factor = pd$pf)
    cf <- as.vector(stats::coef(f))
  } else {
    f <- cv.glmnet(pd$x, y, alpha = alpha, penalty.factor = pd$pf,
                   nfolds = 3L, nlambda = 20L)
    cf <- as.vector(stats::coef(f, s = "lambda.1se"))
    lam_used <- f$lambda.1se
    at_edge <- if (f$lambda.1se >= max(f$lambda)) "upper" else
               if (f$lambda.1se <= min(f$lambda)) "lower" else "interior"
  }
  if (pd$padded) cf <- cf[-length(cf)]
  out[1] <- cf[1]
  out[1L + keep] <- cf[-1]
  attr(out, "lam_used") <- lam_used; attr(out, "lam_at_edge") <- at_edge
  out
}

## ---------------------------------------------------------------------
##  L1-only convention (Zou and Zhang 2009; RAENR-MM): the weights multiply
##  the L1 term only,
##    (1/2n)||y - X b||^2 + lambda alpha sum_j w_j |b_j| + lambda (1-alpha)/2 ||b||^2,
##  on standardised columns, so lambda is on the same scale as glmnet's.
## ---------------------------------------------------------------------
.std <- function(x, y) {
  n <- nrow(x); xm <- colMeans(x); xc <- sweep(x, 2, xm)
  xs <- sqrt(colSums(xc^2) / n); xs[xs < 1e-12] <- 1
  list(x = sweep(xc, 2, xs, "/"), y = y - mean(y), xm = xm, xs = xs, ym = mean(y))
}

## Solved exactly by glmnet's compiled code via data augmentation (Zou and
## Hastie 2005, Lemma 1): appending sqrt(n lam (1 - alpha)) I_p to the
## standardised design and zeros to the centred response turns the problem
## into a weighted lasso with lambda' = lam alpha n / (n + p). Verified to
## 1e-15 against a direct coordinate-descent solution (15 October 2026).
l1only_at <- function(s, w, lam, alpha) {
  n <- nrow(s$x); p <- ncol(s$x)
  if (p == 1L) {                                    # closed form, one column
    z <- sum(s$x * s$y) / n
    return(sign(z) * max(abs(z) - lam * alpha * w, 0) / (1 + lam * (1 - alpha)))
  }
  if (alpha < 1) {
    xa <- rbind(s$x, diag(sqrt(n * lam * (1 - alpha)), p))
    ya <- c(s$y, numeric(p)); lamp <- lam * alpha * n / (n + p)
  } else { xa <- s$x; ya <- s$y; lamp <- lam }
  f <- glmnet(xa, ya, alpha = 1, lambda = lamp, penalty.factor = w,
              intercept = FALSE, standardize = FALSE, thresh = 1e-9)
  as.vector(f$beta)
}

fit_l1only <- function(x, y, alpha, w, keep = NULL,
                       lambda_mode = c("fixed", "retuned"), lam_fixed = NULL,
                       nfolds = 3L, nlam = 20L) {
  lambda_mode <- match.arg(lambda_mode)
  p <- ncol(x)
  if (is.null(keep)) keep <- seq_len(p)
  out <- numeric(p + 1L)
  if (!length(keep)) { out[1] <- mean(y); return(out) }
  xk <- x[, keep, drop = FALSE]
  wk <- w[keep]; wk <- wk / mean(wk)              # unit mean, as glmnet does
  path_fit <- function(xtr, ytr, lams) {
    s <- .std(xtr, ytr); res <- vector("list", length(lams))
    for (i in seq_along(lams)) {
      b <- l1only_at(s, wk, lams[i], alpha)
      bo <- b / s$xs
      res[[i]] <- c(s$ym - sum(s$xm * bo), bo)
    }
    res
  }
  if (lambda_mode == "fixed") {
    cf <- path_fit(xk, y, lam_fixed)[[1]]; lam_used <- lam_fixed; at_edge <- NA
  } else {
    ## grid wide enough to cover both the weighted and the unweighted scale
    s <- .std(xk, y); g <- abs(colSums(s$x * s$y)) / nrow(xk)
    ## lower end follows glmnet's own default: 1e-2 of the unweighted
    ## lambda_max when n < p, 1e-4 otherwise
    ratio <- if (nrow(xk) < ncol(xk)) 1e-2 else 1e-4
    lmax <- max(g / (alpha * wk)); lmin <- ratio * max(g) / alpha
    if (lmin >= lmax) lmin <- ratio * lmax
    lams <- exp(seq(log(lmax), log(lmin), length.out = nlam))
    fold <- sample(rep(seq_len(nfolds), length.out = nrow(xk)))
    err <- matrix(NA_real_, nfolds, nlam)
    for (k in seq_len(nfolds)) {
      tr <- fold != k
      pth <- path_fit(xk[tr, , drop = FALSE], y[tr], lams)
      for (i in seq_len(nlam))
        err[k, i] <- mean((y[!tr] - pth[[i]][1] -
                           as.vector(xk[!tr, , drop = FALSE] %*% pth[[i]][-1]))^2)
    }
    cvm <- colMeans(err); cvsd <- apply(err, 2, stats::sd) / sqrt(nfolds)
    imin <- which.min(cvm)
    i1 <- min(which(cvm <= cvm[imin] + cvsd[imin]))   # largest lambda within 1 SE
    lam_used <- lams[i1]
    at_edge <- if (i1 == 1L) "upper" else if (i1 == nlam) "lower" else "interior"
    cf <- path_fit(xk, y, lams[seq_len(i1)])[[i1]]
  }
  out[1] <- cf[1]; out[1L + keep] <- cf[-1]
  attr(out, "lam_used") <- lam_used; attr(out, "lam_at_edge") <- at_edge
  out
}

## ---------------------------------------------------------------------
##  Metrics
## ---------------------------------------------------------------------
sel_metrics <- function(coefs, beta, b_pilot, x_test = NULL, y_test = NULL) {
  b0 <- coefs[1]; bhat <- coefs[-1]
  p  <- length(beta)
  S  <- which(beta != 0)
  sel <- which(abs(bhat) > 1e-8)
  supp_pilot <- which(b_pilot != 0)
  Z  <- setdiff(seq_len(p), supp_pilot)
  tp <- as.numeric(length(intersect(sel, S)))
  fp <- as.numeric(length(setdiff(sel, S)))
  fn <- as.numeric(length(setdiff(S, sel)))
  tn <- as.numeric(p) - tp - fp - fn
  den <- sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
  un <- union(sel, supp_pilot)
  data.frame(
    nsel     = length(sel),
    tp = tp, fp = fp,
    tpr      = if (length(S)) tp / length(S) else NA_real_,
    fpr      = if (p > length(S)) fp / (p - length(S)) else NA_real_,
    prec     = if (length(sel)) tp / length(sel) else NA_real_,
    f1       = if (tp > 0) 2 * tp / (2 * tp + fp + fn) else 0,
    mcc      = if (den > 0) (tp * tn - fp * fn) / den else 0,
    pilot_nsel = length(supp_pilot),
    pilot_tpr  = if (length(S)) length(intersect(supp_pilot, S)) / length(S) else NA_real_,
    jac_pilot  = if (length(un)) length(intersect(sel, supp_pilot)) / length(un) else 1,
    readmit    = length(intersect(sel, Z)),
    dropped    = length(setdiff(supp_pilot, sel)),
    nsel_minus_pmZ = length(sel) - length(supp_pilot),
    coef_err   = sum((bhat - beta)^2),
    rmse_test  = if (is.null(x_test)) NA_real_ else
                   sqrt(mean((y_test - b0 - as.vector(x_test %*% bhat))^2))
  )
}

jaccard <- function(a, b) {
  u <- union(a, b); if (!length(u)) return(1)
  length(intersect(a, b)) / length(u)
}

mean_pairwise_jaccard <- function(sets) {
  k <- length(sets); if (k < 2L) return(NA_real_)
  s <- 0; m <- 0
  for (i in 1:(k - 1L)) for (j in (i + 1L):k) { s <- s + jaccard(sets[[i]], sets[[j]]); m <- m + 1 }
  s / m
}

## command-line helper: arg i or default
## In RStudio, set e.g.  ARGS <- c("1", "4")  in the console before source().
arg_or <- function(i, default) {
  a <- if (exists("ARGS", envir = globalenv())) get("ARGS", envir = globalenv())
       else commandArgs(trailingOnly = TRUE)
  if (length(a) >= i && nzchar(a[i])) a[i] else default
}

dir.create("results", showWarnings = FALSE)
