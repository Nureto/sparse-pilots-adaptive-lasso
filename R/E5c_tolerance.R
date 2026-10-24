## =====================================================================
##  E5c: numerical tolerance check for E5 Part B (post hoc, 17 October 2026)
##
##  At alpha = 1 the two penalty conventions are the same problem, so E5 Part B
##  should show them agreeing. With glmnet's default convergence threshold they
##  agreed in only 72% of fixed-lambda fits (coefficient differences up to
##  0.04), because the pilot support is almost unpenalised and the problem is
##  poorly conditioned. This script refits the fixed-lambda part of E5 Part B
##  on the SAME data, pilots and lambda (same seeds) with a much tighter
##  threshold, and certifies every solution by its own optimality (KKT)
##  conditions. Fixed lambda only: no cross-validation of the adaptive fit,
##  so it is fast.
##
##  Usage:  Rscript R/E5c_tolerance.R [R]        (default R = 100)
##  Smoke:  ARGS <- c("2"); source("R/E5c_tolerance.R")
## =====================================================================
source("R/common.R")

R     <- as.integer(arg_or(1, "100"))
SEED  <- 20261017L                      # same as E5, so the data are identical
n <- 100L; gamma <- 2; delta <- 1e-6; K <- 3L
TH <- 1e-14                             # convergence threshold for glmnet

gfit <- function(x, y, ...) {           # glmnet with a tight threshold, any version
  f <- tryCatch(glmnet(x, y, ..., control = list(thresh = TH)),
                error = function(e) NULL)
  if (is.null(f)) f <- suppressWarnings(glmnet(x, y, ..., thresh = TH))
  f
}

## KKT residual on standardised columns; ridge multiplier cvec_j per coordinate
kkt_resid <- function(s, bs, w, lam, alpha, cvec) {
  g <- as.vector(crossprod(s$x, s$y - s$x %*% bs)) / nrow(s$x)
  act <- abs(bs) > 0
  r_act <- if (any(act)) max(abs(g[act] - lam * alpha * w[act] * sign(bs[act]) -
                                 cvec[act] * bs[act])) else 0
  r_in  <- if (any(!act)) max(pmax(abs(g[!act]) - lam * alpha * w[!act], 0)) else 0
  max(r_act, r_in)
}

cells <- expand.grid(p = c(100L, 500L), corr = c("equi", "ar1"),
                     scen = c("normal", "t3_y10"), stringsAsFactors = FALSE)
odir <- file.path("results", paste0("E5c_R", R))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

for (ci in seq_len(nrow(cells))) {
  ofile <- file.path(odir, sprintf("cell_%03d.rds", ci))
  if (file.exists(ofile)) { cat("cell", ci, "exists, skipped\n"); next }
  cp <- cells[ci, ]; t0 <- Sys.time(); rows <- list(); k <- 0L
  for (r in seq_len(R)) {
    set.seed(SEED + 100000L * ci + r)               # same stream as E5
    d <- gen_data(n, cp$p, cp$corr, cp$scen)
    foldid <- sample(rep(seq_len(K), length.out = n))
    yc <- d$y - mean(d$y); ys <- yc / sqrt(mean(yc^2))
    s <- .std(d$x, ys)                              # s$y is ys centred (already)
    for (pil in c("lasso", "huber")) {
      b <- tryCatch(run_pilot(pil, d$x, d$y), error = function(e) rep(0, cp$p))
      supp <- which(b != 0)
      w <- adaptive_weights(b, gamma, delta, cap = 0.95)$w   # mean 1
      for (a in c(0.2, 1)) {
        lam <- cv.glmnet(d$x, ys, alpha = a, foldid = foldid, nlambda = 20L)$lambda.1se
        ## convention W: glmnet, penalty.factor = w (y already standardised)
        fW <- gfit(d$x, ys, alpha = a, lambda = c(lam * 10, lam), penalty.factor = w)
        bW <- as.vector(fW$beta[, 2]) * s$xs        # to the standardised scale
        ## convention L: exact augmentation (Zou and Hastie 2005, Lemma 1)
        if (a < 1) {
          xa <- rbind(s$x, diag(sqrt(n * lam * (1 - a)), cp$p))
          ya <- c(s$y, numeric(cp$p)); lamp <- lam * a * n / (n + cp$p)
        } else { xa <- s$x; ya <- s$y; lamp <- lam }
        fL <- gfit(xa, ya, alpha = 1, lambda = c(lamp * 10, lamp), penalty.factor = w,
                   intercept = FALSE, standardize = FALSE)
        bL <- as.vector(fL$beta[, 2])
        selW <- which(abs(bW) > 1e-8); selL <- which(abs(bL) > 1e-8)
        k <- k + 1L
        rows[[k]] <- data.frame(
          cell = ci, p = cp$p, corr = cp$corr, scen = cp$scen, rep = r,
          pilot = pil, alpha = a, pilot_nsel = length(supp),
          kkt_W = kkt_resid(s, bW, w, lam, a, lam * (1 - a) * w),
          kkt_L = kkt_resid(s, bL, w, lam, a, rep(lam * (1 - a), cp$p)),
          coef_diff = max(abs(bW - bL)),
          agree = setequal(selW, selL),
          nsel_W = length(selW), nsel_L = length(selL),
          readmit_W = length(setdiff(selW, supp)), readmit_L = length(setdiff(selL, supp)),
          dropped_W = length(setdiff(supp, selW)), dropped_L = length(setdiff(supp, selL)))
      }
    }
  }
  res <- do.call(rbind, rows); saveRDS(res, ofile)
  cat(sprintf("E5c cell %d  p=%d %-4s %-6s  agree(alpha=1) %.3f  max KKT W %.1e  L %.1e  %.1f min\n",
              ci, cp$p, cp$corr, cp$scen, mean(res$agree[res$alpha == 1]),
              max(res$kkt_W), max(res$kkt_L),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("E5c finished.\n")
