## =====================================================================
##  E3: the robust estimator (RAENR-MM as implemented)      (prediction P6)
##  Pilot {Huber EN sparse, Huber ridge dense} x lambda {fixed, re-tuned}
##
##  Usage:  Rscript R/E3_robust.R <shard> <nshards> [R]
##  Smoke:  Rscript R/E3_robust.R 1 1 2
## =====================================================================
source("R/common.R")

SHARD  <- as.integer(arg_or(1, "1"))
NSHARD <- as.integer(arg_or(2, "1"))
R      <- as.integer(arg_or(3, "100"))
SEED   <- 20261017L
n <- 100L; n_test <- 500L
alpha <- 0.2; gamma <- 2; delta <- 1e-6; cc <- 4.685

## ---- core functions copied from the raenrmm package (as implemented) ----
rho_biweight <- function(u, cc = 4.685) ifelse(abs(u) <= cc, 1 - (1 - (u / cc)^2)^3, 1)
w_biweight   <- function(u, cc = 4.685) ifelse(abs(u) <= cc, (1 - (u / cc)^2)^2, 0)
mscale <- function(u, cc = 4.685, b = 0.5, tol = 1e-8, max_it = 100L) {
  u <- u[is.finite(u)]; if (!length(u)) return(NA_real_)
  s <- stats::mad(u); if (!is.finite(s) || s < 1e-10) return(max(s, 1e-10))
  for (it in seq_len(max_it)) {
    s_new <- s * sqrt(mean(rho_biweight(u / s, cc)) / b)
    if (!is.finite(s_new) || s_new < 1e-10) break
    if (abs(s_new - s) < tol * s) { s <- s_new; break }
    s <- s_new
  }
  s
}
soft_threshold <- function(z, g) sign(z) * pmax(abs(z) - g, 0)
cd_sweep <- function(xw, yw, beta, wp, lam, alpha, d) {
  p <- length(beta); r <- as.vector(yw - xw %*% beta)
  lam_a <- lam * alpha; den_v <- d + lam * (1 - alpha)
  for (j in seq_len(p)) {
    xj <- xw[, j]; r <- r + xj * beta[j]
    den <- den_v[j]
    bnew <- if (abs(den) < 1e-12) 0 else soft_threshold(sum(xj * r) / den, lam_a * wp[j] / den)
    beta[j] <- bnew; r <- r - xj * bnew
  }
  beta
}

## Cross-validation of lambda for the estimator's OWN penalty (weights on the
## L1 term only), scored by the M-scale of held-out residuals as in
## raenrmm::cv_raenrmm. Row weights are the biweight weights of the initial
## residuals; the path is warm-started from large to small lambda.
cv_lambda_robust <- function(x, y, wp, beta0, sigma, nfolds = 3L, nlam = 35L) {
  rw <- w_biweight(as.vector(y - x %*% beta0) / sigma, cc)
  srw <- sqrt(rw); xw <- x * srw; yw <- y * srw; S <- sum(rw)
  g <- abs(colSums(xw * yw)) / S
  ratio <- if (nrow(x) < ncol(x)) 1e-2 else 1e-4   # glmnet's default rule
  lmax <- max(g / (alpha * wp)); lmin <- ratio * max(g) / alpha
  if (lmin >= lmax) lmin <- ratio * lmax
  lams <- exp(seq(log(lmax), log(lmin), length.out = nlam))
  fold <- sample(rep(seq_len(nfolds), length.out = nrow(x)))
  sc <- matrix(NA_real_, nfolds, nlam)
  for (k in seq_len(nfolds)) {
    tr <- fold != k
    xt <- xw[tr, , drop = FALSE]; yt <- yw[tr]; dt <- colSums(xt^2); St <- sum(rw[tr])
    b <- numeric(ncol(x))
    for (i in seq_len(nlam)) {
      for (cd in seq_len(50L)) {
        b_old <- b
        b <- cd_sweep(xt, yt, b, wp, lams[i] * St, alpha, dt)
        if (sqrt(sum((b - b_old)^2)) < 1e-4) break
      }
      sc[k, i] <- mscale(as.vector(y[!tr] - x[!tr, , drop = FALSE] %*% b), cc)
    }
  }
  cvm <- colMeans(sc, na.rm = TRUE)
  i <- which.min(cvm)
  list(lambda = lams[i],
       at_edge = if (i == 1L) "upper" else if (i == nlam) "lower" else "interior")
}

## RAENR-MM with the adaptive weights wp supplied. lambda: "fixed" = the
## implementation's non-adaptive choice; "retuned" = cv_lambda_robust above.
fit_raenrmm <- function(x, y, wp, lambda_mode = c("fixed", "retuned"),
                        max_outer = 5L, max_inner = 50L) {
  lambda_mode <- match.arg(lambda_mode)
  n <- nrow(x); p <- ncol(x)
  pil <- cv.glmnet(x, y, alpha = alpha, nfolds = 3L, nlambda = 20L)
  lambda_use <- pil$lambda.1se
  beta <- safe_coef(pil, "lambda.1se", p)
  r_pil <- as.vector(y - x %*% beta)
  df_pil <- max(1L, sum(abs(beta) > 1e-8))
  sigma <- mscale(r_pil, cc) * sqrt(n / max(n - df_pil, 1L))
  if (!is.finite(sigma) || sigma < 1e-10) sigma <- stats::mad(r_pil)
  frozen <- FALSE; at_edge <- NA
  if (lambda_mode == "retuned") {
    cvr <- cv_lambda_robust(x, y, wp, beta, sigma)
    lambda_use <- cvr$lambda; at_edge <- cvr$at_edge; frozen <- TRUE
  }
  for (it in seq_len(max_outer)) {
    r <- as.vector(y - x %*% beta)
    rw <- w_biweight(r / sigma, cc)
    if (sum(rw) < 1e-6) break
    srw <- sqrt(rw); xw <- x * srw; yw <- y * srw; d_w <- colSums(xw^2)
    keep <- rw > 1e-8
    if (!frozen && sum(keep) >= 12L) {        # as implemented: non-adaptive
      pk <- tryCatch(cv.glmnet(x[keep, , drop = FALSE], y[keep], alpha = alpha,
                               nfolds = 3L, nlambda = 20L), error = function(e) NULL)
      if (!is.null(pk) && is.finite(pk$lambda.1se)) { lambda_use <- pk$lambda.1se; frozen <- TRUE }
    }
    beta_pre <- beta
    for (cd in seq_len(max_inner)) {
      b_old <- beta
      beta <- cd_sweep(xw, yw, beta, wp, lambda_use * sum(rw), alpha, d_w)
      if (sqrt(sum((beta - b_old)^2)) < 1e-4) break
    }
    if (sqrt(sum((beta - beta_pre)^2)) < 1e-3) break
  }
  out <- c(0, beta)                  # intercept = FALSE, as implemented
  attr(out, "lam_used") <- lambda_use; attr(out, "lam_at_edge") <- at_edge
  out
}

## ---- design -------------------------------------------------------------
cells <- expand.grid(p = c(50L, 100L, 200L), contam = c(0, 0.10),
                     stringsAsFactors = FALSE)
pilots <- c("huber", "huber_dense")
lmodes <- c("fixed", "retuned")

odir <- file.path("results", paste0("E3v4_R", R))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
mine <- which((seq_len(nrow(cells)) - 1L) %% NSHARD == (SHARD - 1L))

for (ci in mine) {
  ofile <- file.path(odir, sprintf("cell_%03d.rds", ci))
  if (file.exists(ofile)) { cat("cell", ci, "exists, skipped\n"); next }
  cp <- cells[ci, ]; t0 <- Sys.time()
  rows <- list(); k <- 0L; nfail <- 0L
  for (r in seq_len(R)) {
    set.seed(SEED + 100000L * ci + r)
    d <- gen_data(n, cp$p, "equi", "t3")
    if (cp$contam > 0) {                       # rowwise: whole rows x 10
      bad <- sample.int(n, round(cp$contam * n))
      d$x[bad, ] <- 10 * d$x[bad, ]; d$y[bad] <- 10 * d$y[bad]
    }
    te <- gen_data(n_test, cp$p, "equi", "t3", beta = d$beta)
    for (pil in pilots) {
      b <- tryCatch(run_pilot(pil, d$x, d$y), error = function(e) rep(0, cp$p))
      aw <- adaptive_weights(b, gamma, delta, cap = 0.95)
      wd <- weight_diag(b, aw$w, aw$cap_active)
      for (lm_ in lmodes) {
        cf <- tryCatch(fit_raenrmm(d$x, d$y, aw$w, lm_), error = function(e) NULL)
        if (is.null(cf)) { nfail <- nfail + 1L; next }
        k <- k + 1L
        rows[[k]] <- cbind(
          data.frame(cell = ci, p = cp$p, contam = cp$contam, rep = r,
                     pilot = pil, lambda_mode = lm_,
                     lam_used = attr(cf, "lam_used"), lam_at_edge = attr(cf, "lam_at_edge")),
          wd[, c("Zsize", "Zfrac", "cap_active", "maxwZ", "bound", "log10_ratio")],
          sel_metrics(cf, d$beta, b, te$x, te$y))
      }
    }
  }
  res <- if (length(rows)) do.call(rbind, rows) else data.frame()
  saveRDS(res, ofile)
  cat(sprintf("E3 cell %d  p=%d contam=%.2f  %d rows, %d failed  %.1f min\n", ci, cp$p,
              cp$contam, nrow(res), nfail,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("Shard", SHARD, "finished.\n")
