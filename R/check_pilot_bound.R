## =====================================================================
##  check_pilot_bound.R : can an adaptive lasso or elastic net fit add a
##  variable that its pilot left out?
##
##  check_pilot_bound(x, y, pilot) computes, for a pilot estimate and a data
##  set, the adaptive weights as software builds them (1/(|b|^gamma + delta),
##  95% cap, rescaled to mean one), the effective zero set, and the
##  condition of Proposition 4 of the paper. When the condition holds and
##  lambda is chosen anywhere on glmnet's default path, no variable with a
##  zero pilot coefficient can be selected: the adaptive selection is a
##  subset of the pilot's.
##
##  Arguments
##    x, y      design matrix (n x p) and response, as passed to glmnet
##    pilot     pilot coefficients (length p, intercept excluded), e.g.
##              coef(cv.glmnet(x, y), s = "lambda.min")[-1]
##    gamma, delta, cap   weight construction (defaults as in the paper)
##    lambda.min.ratio    lower end of the path as a fraction of lambda_max;
##              NULL uses glmnet's default (0.01 if n < p, 1e-4 otherwise)
##
##  The columns of x are standardised with divisor n and y is centred, as
##  glmnet does internally (standardize = TRUE, intercept = TRUE).
## =====================================================================
check_pilot_bound <- function(x, y, pilot, gamma = 2, delta = 1e-6, cap = 0.95,
                              lambda.min.ratio = NULL) {
  x <- as.matrix(x); y <- as.numeric(y); b <- as.numeric(pilot)
  n <- nrow(x); p <- ncol(x)
  if (length(b) != p) stop("pilot must have one coefficient per column of x")
  if (is.null(lambda.min.ratio)) lambda.min.ratio <- if (n < p) 1e-2 else 1e-4

  xc <- sweep(x, 2, colMeans(x))
  sx <- sqrt(colMeans(xc^2)); sx[sx == 0] <- 1
  xs <- sweep(xc, 2, sx, "/")
  yc <- y - mean(y)

  ## weights as implemented: raw, capped at the 95% quantile (type 7), mean one
  raw <- 1 / (abs(b)^gamma + delta)
  q <- unname(stats::quantile(raw, cap, type = 7))
  w <- pmin(raw, q); w <- w / mean(w)

  S <- which(b != 0); Z <- which(b == 0)
  kappa <- delta / (abs(b[S])^gamma + delta)
  zeff <- length(Z) + sum(kappa)
  cap_inert <- length(Z) >= p - floor(cap * (p - 1))
  eff_zero <- sum(abs(b[S])^gamma < delta)

  g <- abs(crossprod(xs, yc)) / n                       # |x_k' y| / n
  if (length(S) && length(Z)) {
    lhs <- lambda.min.ratio * max(g[S] * (1 + abs(b[S])^gamma / delta))
    rhs <- sqrt(sum(yc^2)) / sqrt(n)
    holds <- lhs > rhs
  } else {
    lhs <- NA_real_; rhs <- NA_real_; holds <- NA
  }
  out <- list(n = n, p = p, pilot_size = length(S), zeros = length(Z),
              Z_eff = zeff, cap_inert = cap_inert,
              weight_on_zeros = if (length(Z)) w[Z[1]] else NA_real_,
              effective_zeros = eff_zero,
              lambda.min.ratio = lambda.min.ratio, lhs = lhs, rhs = rhs,
              condition_holds = holds, weights = w)
  class(out) <- "pilot_bound"
  out
}

print.pilot_bound <- function(x, ...) {
  cat("Pilot: ", x$pilot_size, " nonzero and ", x$zeros, " zero coefficients (n = ",
      x$n, ", p = ", x$p, ")\n", sep = "")
  if (is.na(x$condition_holds)) {
    if (x$pilot_size == 0)
      cat("The pilot is empty: every weight equals one and the adaptive fit is\n",
          "the ordinary (nonadaptive) elastic net (Lemma 1).\n", sep = "")
    else
      cat("The pilot has no zero coefficients: there is nothing to check.\n")
    return(invisible(x))
  }
  cat(sprintf("Cap inert: %s; common weight on pilot zeros: %.4g (p / Z_eff, Z_eff = %.2f)\n",
              x$cap_inert, x$weight_on_zeros, x$Z_eff))
  cat(sprintf("Pilot coefficients acting as zeros (|b|^gamma < delta): %d\n", x$effective_zeros))
  cat(sprintf("Condition of Proposition 4: %.4g > %.4g ? %s\n", x$lhs, x$rhs,
              if (x$condition_holds) "yes" else "no"))
  if (x$condition_holds)
    cat("With lambda anywhere on glmnet's default path, no pilot zero can be\n",
        "selected: the adaptive selection is a subset of the pilot's.\n", sep = "")
  else
    cat("The condition fails: pilot zeros may enter near the bottom of the path.\n")
  invisible(x)
}
