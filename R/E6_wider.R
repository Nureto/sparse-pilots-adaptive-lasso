## =====================================================================
##  E6: wider design and direct checks of the theory   (added 17 October 2026)
##
##  Three things in one run, on new replications:
##   1. Wider design: n in {50, 100, 200}, strong or weak signal (true
##      coefficients 3 -> 1.5, or half that), p in {200, 1000},
##      equicorrelated or AR(1) predictors, t3 errors. 24 cells.
##      The n = 100, strong-signal cells repeat E2 conditions on new data.
##   2. Assumption 2 checked directly: the limit problems F_W^0 and F_L^0 of
##      Proposition 2 are solved at the same fixed lambda, and we record
##      whether every pilot variable is nonzero in the limit solution.
##   3. Which pilot variables a fixed lambda drops, with their pilot
##      coefficients, so the effective-zero-set explanation is tested.
##   Also: condition (5) of Proposition 4 is computed for every fit.
##
##  The response is standardised (1/n variance) so that glmnet solves the
##  stated convention-W problem exactly (see check_ystd.R).
##
##  Usage:  Rscript R/E6_wider.R <shard> <nshards> [R]
##  Smoke:  ARGS <- c("1", "1", "2"); source("R/E6_wider.R")
## =====================================================================
source("R/common.R")
## Never stop a glmnet path early. By default glmnet stops when 99.9% of the
## deviance is explained, which happens when the pilot support is nearly
## unpenalised and has about n or more variables (first run, 17 October:
## 292 fits failed this way, mostly at n = 50). Reset at the end of the script.
glmnet.control(fdev = 0, devmax = 1)

SHARD  <- as.integer(arg_or(1, "1"))
NSHARD <- as.integer(arg_or(2, "1"))
R      <- as.integer(arg_or(3, "100"))
SEED   <- 20261018L
gamma <- 2; delta <- 1e-6; alpha <- 0.2; K <- 3L
deltas <- c(1e-6, 1e-8, 1e-10, 1e-12)
dlab   <- c("1e-6", "1e-8", "1e-10", "1e-12")       # labels used in column names
TH <- 1e-12

gfit <- function(x, y, ...) {           # glmnet with a tight threshold
  f <- tryCatch(glmnet(x, y, ..., control = list(thresh = TH)),
                error = function(e) NULL)
  if (is.null(f)) f <- suppressWarnings(glmnet(x, y, ..., thresh = TH))
  f
}
## If the tight threshold is not reached within glmnet's iteration limit (this
## happens when the pilot support is nearly unpenalised and has about n or
## more variables), glmnet returns only the first lambda. Then refit at
## glmnet's default threshold and flag the fit as "loose" (v3, 17 October).
gpath <- function(x, y, ...) {
  f <- gfit(x, y, ...); loose <- FALSE
  if (ncol(f$beta) < 2L) { f <- suppressWarnings(glmnet(x, y, ...)); loose <- TRUE }
  if (ncol(f$beta) < 2L) stop("glmnet path stopped before the target lambda")
  b <- as.vector(f$beta[, 2]); attr(b, "loose") <- loose; b
}
safe <- function(expr) tryCatch(expr, error = function(e) NULL)
isloose <- function(b) isTRUE(attr(b, "loose"))
## convention W on standardised data, fixed lambda (path of two values)
fitW <- function(s, pf, lam)
  gpath(s$x, s$y, alpha = alpha, lambda = c(10 * lam, lam),
        penalty.factor = pf, standardize = FALSE, intercept = FALSE)
## convention L by exact augmentation (Zou and Hastie 2005, Lemma 1)
fitL <- function(s, pf, lam) {
  n <- nrow(s$x); p <- ncol(s$x)
  xa <- rbind(s$x, diag(sqrt(n * lam * (1 - alpha)), p))
  ya <- c(s$y, numeric(p)); lamp <- lam * alpha * n / (n + p)
  gpath(xa, ya, alpha = 1, lambda = c(10 * lamp, lamp),
        penalty.factor = pf, standardize = FALSE, intercept = FALSE)
}
nz <- function(b) which(abs(b) > 1e-8)

cells <- expand.grid(n = c(50L, 100L, 200L), p = c(200L, 1000L),
                     signal = c("strong", "weak"), corr = c("equi", "ar1"),
                     stringsAsFactors = FALSE)
odir <- file.path("results", paste0("E6v2_R", R))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
mine <- which((seq_len(nrow(cells)) - 1L) %% NSHARD == (SHARD - 1L))

for (ci in mine) {
  ofile <- file.path(odir, sprintf("cell_%03d.rds", ci))
  if (file.exists(ofile)) { cat("cell", ci, "exists, skipped\n"); next }
  cp <- cells[ci, ]; n <- cp$n; p <- cp$p; t0 <- Sys.time()
  rows <- list(); k <- 0L; n_empty <- c(lasso = 0L, huber = 0L); nfail <- 0L
  for (r in seq_len(R)) {
    set.seed(SEED + 100000L * ci + r)
    beta <- true_beta(p) * (if (cp$signal == "weak") 0.5 else 1)
    d <- gen_data(n, p, cp$corr, "t3", beta = beta)
    truth <- which(d$beta != 0)
    foldid <- sample(rep(seq_len(K), length.out = n))
    yc <- d$y - mean(d$y); ys <- yc / sqrt(mean(yc^2))
    s <- .std(d$x, ys)
    lam <- cv.glmnet(d$x, ys, alpha = alpha, foldid = foldid,
                     nlambda = 20L)$lambda.1se          # as implemented
    g <- abs(colSums(s$x * s$y)) / n
    for (pil in c("lasso", "huber")) {
      res <- tryCatch({
        b <- run_pilot(pil, d$x, ys)
        S <- which(b != 0); Z <- setdiff(seq_len(p), S)
        if (!length(S)) { n_empty[pil] <- n_empty[pil] + 1L; NULL } else {
          met <- function(sel, pre) {
            if (is.null(sel)) {
              o <- as.list(rep(NA_real_, 6))
              names(o) <- paste0(pre, c("_nsel", "_tpr", "_prec", "_jac", "_drop", "_add"))
              return(o)
            }
            tp <- length(intersect(sel, truth))
            o <- list(length(sel), tp / length(truth),
                      if (length(sel)) tp / length(sel) else NA_real_,
                      jaccard(sel, S), length(setdiff(S, sel)), length(setdiff(sel, S)))
            names(o) <- paste0(pre, c("_nsel", "_tpr", "_prec", "_jac", "_drop", "_add"))
            o
          }
          aw <- adaptive_weights(b, gamma, delta, cap = 0.95); w <- aw$w
          kap <- delta / (abs(b)^gamma + delta)
          eps <- if (n < p) 1e-2 else 1e-4
          Cl <- eps * max(g[S] * (1 + abs(b[S])^gamma / delta))
          Cr <- sqrt(mean(s$y^2))
          L <- list(cell = ci, n = n, p = p, signal = cp$signal, corr = cp$corr,
                    rep = r, pilot = pil, lam = lam, nS = length(S),
                    pilot_tpr = length(intersect(S, truth)) / length(truth),
                    cap_active = aw$cap_active, Zeff = sum(kap),
                    n_effzero = sum(kap[S] >= 0.5), condC = Cl > Cr, condC_ratio = Cl / Cr)
          ## retuned lambda, convention W, glmnet default path (as in E2)
          fr <- cv.glmnet(s$x, s$y, alpha = alpha, penalty.factor = w, foldid = foldid,
                          nlambda = 20L, standardize = FALSE, intercept = FALSE)
          selR <- nz(as.vector(coef(fr, s = "lambda.1se"))[-1])
          L <- c(L, met(selR, "Wret"),
                 list(Wret_edge = if (fr$lambda.1se >= max(fr$lambda)) "upper" else
                      if (fr$lambda.1se <= min(fr$lambda)) "lower" else "interior"))
          ## fixed lambda, convention W, at four values of delta
          bW <- lapply(deltas, function(dl)
            safe(fitW(s, adaptive_weights(b, gamma, dl, cap = 0.95)$w, lam)))
          selW <- if (is.null(bW[[1]])) NULL else nz(bW[[1]])   # delta = 1e-6
          L <- c(L, met(selW, "Wfix"))
          ## fixed lambda, convention L, delta = 1e-6
          bL <- safe(fitL(s, w, lam)); selL <- if (is.null(bL)) NULL else nz(bL)
          L <- c(L, met(selL, "Lfix"))
          L$fix_failed <- sum(vapply(c(bW, list(bL)), is.null, logical(1)))
          L$fix_loose  <- sum(vapply(c(bW, list(bL)), isloose, logical(1)))
          ## limit problems (delta -> 0): pilot support unpenalised in the l1 part
          L[c("rankS_full", "A2_W", "A2_L", "minS_W", "minS_L", "W0_nsel", "L0_nsel",
              "W0_add", "L0_add", "W0_margin", "W_same_1e-6", "W_same_1e-8",
              "W_same_1e-10", "W_same_1e-12", "L_same_1e-6", "L0_loose", "W0_loose")] <- NA
          if (length(Z)) {
            pf0 <- numeric(p); pf0[Z] <- p / length(Z)
            L$rankS_full <- length(S) < n && qr(s$x[, S, drop = FALSE])$rank == length(S)
            b0L <- safe(fitL(s, pf0, lam))
            if (!is.null(b0L)) {
              s0L <- nz(b0L); L$L0_loose <- isloose(b0L)
              L$A2_L <- all(abs(b0L[S]) > 1e-8); L$minS_L <- min(abs(b0L[S]))
              L$L0_nsel <- length(s0L); L$L0_add <- length(setdiff(s0L, S))
              if (!is.null(selL)) L[["L_same_1e-6"]] <- setequal(selL, s0L)
            }
            b0W <- if (isTRUE(L$rankS_full)) safe(fitW(s, pf0, lam)) else NULL
            if (!is.null(b0W)) {
              s0W <- nz(b0W); L$W0_loose <- isloose(b0W)
              L$A2_W <- all(abs(b0W[S]) > 1e-8); L$minS_W <- min(abs(b0W[S]))
              L$W0_nsel <- length(s0W); L$W0_add <- length(setdiff(s0W, S))
              r0 <- as.vector(s$y - s$x %*% b0W)
              thr <- lam * alpha * p / length(Z)
              inZ <- Z[abs(b0W[Z]) <= 1e-8]
              if (length(inZ))
                L$W0_margin <- min(thr - abs(colSums(s$x[, inZ, drop = FALSE] * r0)) / n) / thr
              for (i in seq_along(deltas)) if (!is.null(bW[[i]]))
                L[[paste0("W_same_", dlab[i])]] <- setequal(nz(bW[[i]]), s0W)
            }
          }
          ## the pilot variables dropped at a fixed lambda (W and L), and why
          for (cv_ in c("W", "L")) {
            sl <- if (cv_ == "W") selW else selL
            dr <- if (is.null(sl)) integer(0) else setdiff(S, sl)
            L[[paste0(cv_, "drop_effzero")]]   <- sum(kap[dr] >= 0.5)      # |b|^2 <= delta
            L[[paste0(cv_, "drop_near")]]      <- sum(kap[dr] >= 0.01)     # |b|^2 <= 99 delta
            L[[paste0(cv_, "drop_maxabs")]]    <- if (length(dr)) max(abs(b[dr])) else NA_real_
            L[[paste0(cv_, "drop_medabs")]]    <- if (length(dr)) stats::median(abs(b[dr])) else NA_real_
          }
          L$pilot_medabs <- stats::median(abs(b[S]))
          as.data.frame(L, check.names = FALSE, stringsAsFactors = FALSE)
        }
      }, error = function(e) { message("cell ", ci, " rep ", r, " ", pil, ": ",
                                       conditionMessage(e)); "fail" })
      if (identical(res, "fail")) nfail <- nfail + 1L
      else if (is.data.frame(res)) { k <- k + 1L; rows[[k]] <- res }
    }
  }
  df <- if (k) do.call(rbind, rows) else NULL
  saveRDS(list(rows = df, n_empty = n_empty, n_failed = nfail), ofile)
  cat(sprintf("E6 cell %2d  n=%3d p=%4d %-6s %-4s  %d rows, empty %d/%d, failed %d  %.1f min\n",
              ci, n, p, cp$signal, cp$corr, k, n_empty["lasso"], n_empty["huber"], nfail,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
glmnet.control(factory = TRUE)
cat("Shard", SHARD, "finished.\n")
