## =====================================================================
##  E7: why does the robust estimator add so many pilot zeros? (17 Oct 2026)
##
##  E3 found that RAENR-MM, at a fixed lambda, added about 29 pilot zeros,
##  whereas E5c found only about 1 for squared loss once the response was
##  standardised. Hypothesis, recorded before this run: the implementation
##  chooses lambda with glmnet, which standardises y internally, and then
##  applies it to its own objective on the original response, so its ridge
##  term is about sd(y) times stronger than the one lambda was chosen for.
##  Test: on the SAME data and the SAME pilot, run the estimator exactly as
##  implemented (a) on y as observed and (b) on y standardised (1/n variance).
##  If (b) adds far fewer pilot zeros, the hypothesis is supported.
##  The estimator's code is taken unchanged from R/E3_robust.R.
##
##  Usage:  Rscript R/E7_robust_scaled.R [R]        (default R = 100)
##  Smoke:  ARGS <- c("2"); source("R/E7_robust_scaled.R")
## =====================================================================
R_E7 <- as.integer(if (exists("ARGS", envir = globalenv()))
  get("ARGS", envir = globalenv())[1] else commandArgs(trailingOnly = TRUE)[1])
if (is.na(R_E7)) R_E7 <- 100L

## load the estimator's functions from E3 without running E3's design
src <- readLines("R/E3_robust.R")
cut <- grep("^## ---- design", src)[1]
ARGS_saved <- if (exists("ARGS", envir = globalenv())) get("ARGS", envir = globalenv()) else NULL
assign("ARGS", c("1", "1", "1"), envir = globalenv())
eval(parse(text = src[seq_len(cut - 1L)]))
if (is.null(ARGS_saved)) rm("ARGS", envir = globalenv()) else assign("ARGS", ARGS_saved, envir = globalenv())

SEED7 <- 20261019L
cells <- expand.grid(p = c(50L, 100L, 200L), contam = c(0, 0.10),
                     stringsAsFactors = FALSE)
odir <- file.path("results", paste0("E7v2_R", R_E7))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)

for (ci in seq_len(nrow(cells))) {
  ofile <- file.path(odir, sprintf("cell_%03d.rds", ci))
  if (file.exists(ofile)) { cat("cell", ci, "exists, skipped\n"); next }
  cp <- cells[ci, ]; t0 <- Sys.time(); rows <- list(); k <- 0L; nfail <- 0L
  for (r in seq_len(R_E7)) {
    set.seed(SEED7 + 100000L * ci + r)
    d <- gen_data(n, cp$p, "equi", "t3")
    if (cp$contam > 0) {                       # rowwise, as in E3
      bad <- sample.int(n, round(cp$contam * n))
      d$x[bad, ] <- 10 * d$x[bad, ]; d$y[bad] <- 10 * d$y[bad]
    }
    sy <- sqrt(mean((d$y - mean(d$y))^2))
    b <- tryCatch(run_pilot("huber", d$x, d$y), error = function(e) rep(0, cp$p))
    S <- which(b != 0); truth <- which(d$beta != 0)
    ## v2 (17 October): a third version keeps the weights of the observed
    ## scale, so that only the scaling of the objective changes. In the
    ## second version the pilot is rescaled too, which also moves pilot
    ## coefficients into the effective zero set (|b|^2 <= delta).
    for (ver in c("observed", "standardised", "standardised_sameweights")) {
      yv <- if (ver == "observed") d$y else (d$y - mean(d$y)) / sy
      bv <- if (ver == "standardised") b / sy else b   # pilot on the y scale used
      w <- adaptive_weights(bv, gamma, delta, cap = 0.95)$w
      set.seed(SEED7 + 100000L * ci + r + 50000L)     # same CV folds in both versions
      cf <- tryCatch(fit_raenrmm(d$x, yv, w, "fixed"), error = function(e) NULL)
      if (is.null(cf)) { nfail <- nfail + 1L; next }
      sel <- which(abs(cf[-1]) > 1e-8); tp <- length(intersect(sel, truth))
      k <- k + 1L
      rows[[k]] <- data.frame(cell = ci, p = cp$p, contam = cp$contam, rep = r,
        version = ver, sd_y = sy, lam_used = attr(cf, "lam_used"),
        pilot_nsel = length(S), nsel = length(sel), tpr = tp / length(truth),
        prec = if (length(sel)) tp / length(sel) else NA_real_,
        jac_pilot = jaccard(sel, S), dropped = length(setdiff(S, sel)),
        added = length(setdiff(sel, S)))
    }
  }
  res <- do.call(rbind, rows); saveRDS(res, ofile)
  a <- tapply(res$added, res$version, mean)
  cat(sprintf("E7 cell %d  p=%d contam=%.2f  added: observed %.1f, standardised %.1f, same weights %.1f  failed %d  %.1f min\n",
              ci, cp$p, cp$contam, a["observed"], a["standardised"], a["standardised_sameweights"], nfail,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("E7 finished.\n")
