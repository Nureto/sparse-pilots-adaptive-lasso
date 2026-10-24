## =====================================================================
##  check_ystd.R  (17 October 2026; a few seconds)
##  Which problem does glmnet solve for the Gaussian family when alpha < 1?
##  Its help page says it standardises y internally (1/n variance) and then
##  unstandardises the coefficients. If so, the ridge term it applies is
##  lambda * (1 - alpha) / sd_y rather than lambda * (1 - alpha).
##  This script checks the KKT conditions of both candidate problems on the
##  active set. The one with residual near 0 (about 1e-7) is the one solved.
##  Usage: source("R/check_ystd.R")
## =====================================================================
library(glmnet)
set.seed(1)
n <- 100; p <- 50; alpha <- 0.2
x <- matrix(rnorm(n * p), n, p)
x <- scale(x, center = TRUE, scale = FALSE)
x <- sweep(x, 2, sqrt(colMeans(x^2)), "/")          # 1/n unit variance
y <- as.vector(x[, 1:5] %*% c(3, 2.5, 2, 1.5, 1) + rnorm(n) * 2)
y <- y - mean(y)
sdy <- sqrt(mean(y^2))                               # about 4 to 5 here
w <- runif(p, 0.5, 1.5); w <- w * p / sum(w)         # weights summing to p
lam <- 0.1

fit <- suppressWarnings(glmnet(x, y, alpha = alpha, lambda = lam,
              penalty.factor = w, standardize = FALSE, intercept = FALSE,
              thresh = 1e-14))
## The two candidate ridge terms differ by a factor of about 4 here, so the
## check discriminates clearly even at glmnet's default convergence threshold.
b <- as.vector(coef(fit))[-1]
act <- which(b != 0)
g <- as.vector(crossprod(x, y - x %*% b)) / n

kkt <- function(cfac)                                # ridge multiplier c
  max(abs(g[act] - lam * alpha * w[act] * sign(b[act]) -
          lam * (1 - alpha) * cfac * w[act] * b[act]))

cat(sprintf("sd(y) = %.3f, active variables = %d\n", sdy, length(act)))
cat(sprintf("KKT residual, ridge = lambda(1-alpha)        : %.2e\n", kkt(1)))
cat(sprintf("KKT residual, ridge = lambda(1-alpha)/sd(y)  : %.2e\n", kkt(1 / sdy)))
cat("glmnet version:", as.character(packageVersion("glmnet")), "\n")
