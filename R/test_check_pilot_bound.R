## =====================================================================
##  test_check_pilot_bound.R : checks check_pilot_bound() against glmnet
##  Usage: source("R/test_check_pilot_bound.R")   (run from PilotWeights_Q)
##  About 1 to 3 minutes. For 50 random data sets it fits the lasso pilot,
##  applies check_pilot_bound(), then fits the adaptive elastic net on
##  glmnet's whole default path and counts any pilot zero that enters.
##  Expected last line: "pilot zeros entered when the condition held: 0".
## =====================================================================
suppressPackageStartupMessages(library(glmnet))
source("R/check_pilot_bound.R")

set.seed(20261023)
held <- 0L; entered_when_held <- 0L; entered_any <- 0L; N <- 50L
for (t in seq_len(N)) {
  n <- 100L; p <- 200L
  x <- matrix(rnorm(n * p), n, p) + 0.5 * rnorm(n)          # equicorrelated
  beta <- numeric(p); beta[round(seq(1, p, length.out = 10))] <- seq(3, 1.5, length.out = 10)
  y <- as.vector(x %*% beta) + rt(n, df = 3)
  pilot <- as.vector(coef(cv.glmnet(x, y, nfolds = 5), s = "lambda.min"))[-1]
  chk <- check_pilot_bound(x, y, pilot)
  if (is.na(chk$condition_holds)) next
  fit <- glmnet(x, y, alpha = 0.2, penalty.factor = chk$weights)   # default path
  Z <- which(pilot == 0)
  ent <- any(as.matrix(fit$beta)[Z, , drop = FALSE] != 0)
  entered_any <- entered_any + ent
  if (chk$condition_holds) { held <- held + 1L; entered_when_held <- entered_when_held + ent }
}
cat("Example output for the last data set:\n"); print(chk)
cat(sprintf("\ncondition held in %d of %d data sets; a pilot zero entered somewhere on the path in %d data sets overall\n",
            held, N, entered_any))
cat(sprintf("pilot zeros entered when the condition held: %d\n", entered_when_held))
