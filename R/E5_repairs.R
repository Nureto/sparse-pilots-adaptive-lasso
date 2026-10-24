## =====================================================================
##  E5: two POST HOC repairs (17 October 2026), labelled as such in the paper
##
##  A. P4 repaired. Finite-delta weights against the delta = 0 theory, both
##     re-tuned, now with the SAME cross-validation folds and with lambda grids
##     matched so that the retained coefficients face identical penalties. Any
##     remaining disagreement is due to the zero-set variables alone.
##  B. Penalty conventions compared on a COMMON lambda grid, the same folds and
##     the same 1-SE rule, and also at a common FIXED lambda. y is standardised
##     first (1/n variance): glmnet standardises y internally for the Gaussian
##     family, which divides its effective ridge term by sd(y), so without this
##     step the two conventions would also differ in ridge strength (this
##     confounded the E2 comparison). At alpha = 1 the two conventions are the
##     same problem, so they must agree: that is the built-in check.
##
##  Usage:  Rscript R/E5_repairs.R <shard> <nshards> [R]
##  Smoke:  ARGS <- c("1", "1", "2"); source("R/E5_repairs.R")
## =====================================================================
source("R/common.R")

SHARD  <- as.integer(arg_or(1, "1"))
NSHARD <- as.integer(arg_or(2, "1"))
R      <- as.integer(arg_or(3, "100"))
SEED   <- 20261017L
n <- 100L; gamma <- 2; delta <- 1e-6; NLAM <- 20L; K <- 3L

## ---- one CV engine shared by every solver ---------------------------------
## solver(x, y, lams) returns a list of c(intercept, beta) for each lambda.
cv_common <- function(x, y, solver, lams, foldid) {
  err <- matrix(NA_real_, K, length(lams))
  for (k in seq_len(K)) {
    tr <- foldid != k
    pth <- solver(x[tr, , drop = FALSE], y[tr], lams)
    for (i in seq_along(lams))
      err[k, i] <- mean((y[!tr] - pth[[i]][1] -
                         as.vector(x[!tr, , drop = FALSE] %*% pth[[i]][-1]))^2)
  }
  cvm <- colMeans(err); cvsd <- apply(err, 2, stats::sd) / sqrt(K)
  imin <- which.min(cvm)
  i1 <- min(which(cvm <= cvm[imin] + cvsd[imin]))          # 1-SE rule
  list(coef = solver(x, y, lams[seq_len(i1)])[[i1]], i = i1,
       edge = if (i1 == 1L) "upper" else if (i1 == length(lams)) "lower" else "interior")
}

## glmnet convention on a given grid (weights multiply the whole penalty)
solver_glmnet <- function(w, alpha, keep = NULL) function(x, y, lams) {
  p <- ncol(x); if (is.null(keep)) keep <- seq_len(p)
  out0 <- numeric(p + 1L)
  if (!length(keep)) { out0[1] <- mean(y); return(rep(list(out0), length(lams))) }
  pd <- .pad(x[, keep, drop = FALSE], w[keep])
  f <- glmnet(pd$x, y, alpha = alpha, lambda = lams, penalty.factor = pd$pf)
  cf <- as.matrix(stats::coef(f))
  if (pd$padded) cf <- cf[-nrow(cf), , drop = FALSE]
  if (ncol(cf) < length(lams))                       # glmnet may stop early
    cf <- cbind(cf, cf[, rep(ncol(cf), length(lams) - ncol(cf))])
  lapply(seq_along(lams), function(i) { o <- out0; o[1] <- cf[1, i]; o[1L + keep] <- cf[-1, i]; o })
}

## l1-only convention on a given grid (weights multiply the l1 term only)
solver_l1 <- function(w, alpha) function(x, y, lams) {
  wk <- w / mean(w); s <- .std(x, y)
  lapply(lams, function(l) {
    b <- l1only_at(s, wk, l, alpha) / s$xs
    c(s$ym - sum(s$xm * b), b)
  })
}

make_grid <- function(x, y, w, alpha) {
  s <- .std(x, y); g <- abs(colSums(s$x * s$y)) / nrow(x)
  ratio <- if (nrow(x) < ncol(x)) 1e-2 else 1e-4
  lmax <- max(g / (alpha * w)); lmin <- ratio * max(g) / alpha
  if (lmin >= lmax) lmin <- ratio * lmax
  exp(seq(log(lmax), log(lmin), length.out = NLAM))
}

## ---- design -------------------------------------------------------------
cells <- expand.grid(p = c(100L, 500L), corr = c("equi", "ar1"),
                     scen = c("normal", "t3_y10"), stringsAsFactors = FALSE)
odir <- file.path("results", paste0("E5v2_R", R))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
mine <- which((seq_len(nrow(cells)) - 1L) %% NSHARD == (SHARD - 1L))

for (ci in mine) {
  ofile <- file.path(odir, sprintf("cell_%03d.rds", ci))
  if (file.exists(ofile)) { cat("cell", ci, "exists, skipped\n"); next }
  cp <- cells[ci, ]; t0 <- Sys.time(); rows <- list(); k <- 0L
  nfail <- 0L; nempty <- 0L
  for (r in seq_len(R)) {
    set.seed(SEED + 100000L * ci + r)
    d <- gen_data(n, cp$p, cp$corr, cp$scen)
    foldid <- sample(rep(seq_len(K), length.out = n))       # shared by all fits
    for (pil in c("lasso", "huber")) {
      b <- tryCatch(run_pilot(pil, d$x, d$y), error = function(e) rep(0, cp$p))
      supp <- which(b != 0)
      w <- adaptive_weights(b, gamma, delta, cap = 0.95)$w

      ## ---- A: finite delta vs theory, alpha = 0.2, glmnet convention ----
      res <- if (length(supp) < 1L) { nempty <- nempty + 1L; NA } else tryCatch({
        wth <- numeric(cp$p); wth[supp] <- 1 / abs(b[supp])^gamma
        wth[supp] <- wth[supp] / mean(wth[supp])
        grid_th <- make_grid(d$x[, supp, drop = FALSE], d$y, wth[supp], 0.2)
        cfac <- stats::median(w[supp] / wth[supp])          # common rescaling
        f_fin <- cv_common(d$x, d$y, solver_glmnet(w, 0.2), grid_th / cfac, foldid)
        f_th  <- cv_common(d$x, d$y, solver_glmnet(wth, 0.2, keep = supp), grid_th, foldid)
        s_fin <- which(abs(f_fin$coef[-1]) > 1e-8); s_th <- which(abs(f_th$coef[-1]) > 1e-8)
        data.frame(part = "A", cell = ci, p = cp$p, corr = cp$corr, scen = cp$scen,
                   rep = r, pilot = pil, alpha = 0.2, pen_conv = "glmnet",
                   pilot_nsel = length(supp),
                   nsel_fin = length(s_fin), nsel_th = length(s_th),
                   same = setequal(s_fin, s_th),
                   readmit_fin = length(setdiff(s_fin, supp)),
                   i_fin = f_fin$i, i_th = f_th$i)
      }, error = function(e) NULL)
      if (is.null(res)) nfail <- nfail + 1L
      else if (is.data.frame(res)) { k <- k + 1L; rows[[k]] <- res }

      ## ---- B: conventions on standardised y, alpha 0.2 and 1 ----
      yc <- d$y - mean(d$y); ys <- yc / sqrt(mean(yc^2))
      for (a in c(0.2, 1)) {
        grid <- make_grid(d$x, ys, w, a)
        fits <- tryCatch({
          lam_fix <- cv.glmnet(d$x, ys, alpha = a, foldid = foldid,
                               nlambda = 20L)$lambda.1se       # non-adaptive, as implemented
          fx <- list(coef = NULL, i = NA_integer_, edge = NA_character_)
          list(retuned = list(
                 glmnet = cv_common(d$x, ys, solver_glmnet(w, a), grid, foldid),
                 l1only = cv_common(d$x, ys, solver_l1(w, a), grid, foldid)),
               fixed = list(
                 glmnet = modifyList(fx, list(coef = solver_glmnet(w, a)(d$x, ys, c(lam_fix * 10, lam_fix))[[2]])),
                 l1only = modifyList(fx, list(coef = solver_l1(w, a)(d$x, ys, lam_fix)[[1]]))))
        }, error = function(e) NULL)
        if (is.null(fits)) { nfail <- nfail + 1L; next }
        for (lm_ in names(fits)) {
          ff <- fits[[lm_]]
          sel <- lapply(ff, function(f) which(abs(f$coef[-1]) > 1e-8))
          cdiff <- max(abs(ff$glmnet$coef[-1] - ff$l1only$coef[-1]))
          for (cv_ in names(ff)) {
            k <- k + 1L
            rows[[k]] <- cbind(
              data.frame(part = "B", cell = ci, p = cp$p, corr = cp$corr, scen = cp$scen,
                         rep = r, pilot = pil, alpha = a, lambda_mode = lm_, pen_conv = cv_,
                         agree = setequal(sel$glmnet, sel$l1only), coef_diff = cdiff,
                         lam_index = ff[[cv_]]$i, lam_edge = ff[[cv_]]$edge),
              sel_metrics(ff[[cv_]]$coef, d$beta / sqrt(mean(yc^2)), b))
          }
        }
      }
    }
  }
  ## parts A and B have different columns, so they are saved as two tables
  isA <- vapply(rows, function(z) z$part[1] == "A", logical(1))
  saveRDS(list(A = if (any(isA)) do.call(rbind, rows[isA]) else NULL,
               B = if (any(!isA)) do.call(rbind, rows[!isA]) else NULL,
               n_empty_pilot = nempty, n_failed = nfail), ofile)
  cat(sprintf("E5 cell %d  p=%d %-4s %-6s  %d records, %d empty pilots, %d failed  %.1f min\n",
              ci, cp$p, cp$corr, cp$scen, length(rows), nempty, nfail,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("Shard", SHARD, "finished.\n")
