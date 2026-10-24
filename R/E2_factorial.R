## =====================================================================
##  E2: pilot x weights x lambda factorial, squared loss  (predictions P2-P5)
##
##  Usage (from the PilotWeights_Q folder):
##    Rscript R/E2_factorial.R <shard> <nshards> [R] [design]
##      design = "main"  (alpha = 0.2, full grid)            default
##             = "lasso" (alpha = 1, AR(1), normal, p in {100, 500, 1000})
##  Examples:
##    Rscript R/E2_factorial.R 1 1 2 main      # smoke test, about 5 minutes
##    Rscript R/E2_factorial.R 1 4             # shard 1 of 4, R = 200
##
##  Each design cell is saved to its own file as soon as it finishes, and
##  finished cells are skipped, so an interrupted run can simply be restarted.
## =====================================================================
source("R/common.R")

SHARD  <- as.integer(arg_or(1, "1"))
NSHARD <- as.integer(arg_or(2, "1"))
R      <- as.integer(arg_or(3, "200"))
DESIGN <- arg_or(4, "main")
SEED   <- 20261016L
n      <- 100L
n_test <- 500L
gamma  <- 2
delta  <- 1e-6

if (DESIGN == "main") {
  alpha <- 0.2
  cells <- expand.grid(p = c(50L, 100L, 200L, 500L, 1000L),
                       corr = c("equi", "ar1", "block"),
                       scen = c("normal", "t3", "t3_y10"),
                       stringsAsFactors = FALSE)
} else if (DESIGN == "lasso") {
  alpha <- 1
  cells <- expand.grid(p = c(100L, 500L, 1000L), corr = "ar1", scen = "normal",
                       stringsAsFactors = FALSE)
} else stop("design must be 'main' or 'lasso'")

nz <- function(v) if (is.null(v)) NA else v
pilots   <- c("lasso", "huber", "ridge")
wschemes <- c("cap", "nocap", "theory")
lmodes   <- c("fixed", "retuned")

odir <- file.path("results", paste0("E2v4_", DESIGN, "_R", R))
dir.create(odir, showWarnings = FALSE, recursive = TRUE)
mine <- which((seq_len(nrow(cells)) - 1L) %% NSHARD == (SHARD - 1L))
cat(sprintf("E2 design=%s  shard %d/%d  cells %s  R=%d\n", DESIGN, SHARD, NSHARD,
            paste(mine, collapse = ","), R))

for (ci in mine) {
  ofile <- file.path(odir, sprintf("cell_%03d.rds", ci))
  if (file.exists(ofile)) { cat("cell", ci, "exists, skipped\n"); next }
  cp <- cells[ci, ]
  t0 <- Sys.time()
  rows <- list(); k <- 0L; fails <- list()
  for (r in seq_len(R)) {
    set.seed(SEED + 100000L * ci + r)
    d  <- gen_data(n, cp$p, cp$corr, cp$scen)
    te <- gen_data(n_test, cp$p, cp$corr,
                   if (cp$scen == "normal") "normal" else "t3", beta = d$beta)
    lam_fixed <- tryCatch(lambda_fixed_nonadaptive(d$x, d$y, alpha),
                          error = function(e) NA_real_)
    for (pil in pilots) {
      b <- tryCatch(run_pilot(pil, d$x, d$y), error = function(e) rep(0, cp$p))
      supp <- which(b != 0)
      for (ws in wschemes) {
        if (ws == "theory") {
          w <- numeric(cp$p); w[supp] <- 1 / abs(b[supp])^gamma
          keep <- supp; cap_act <- NA
          wd <- weight_diag(b, w, FALSE, cap = NULL)
          ## theory: pilot zeros carry an infinite weight (excluded outright)
          if (wd$Zsize > 0) { wd$maxwZ <- Inf; wd$minwZ <- Inf; wd$log10_ratio <- Inf }
        } else {
          aw <- adaptive_weights(b, gamma, delta, cap = if (ws == "cap") 0.95 else NULL)
          w <- aw$w; keep <- NULL; cap_act <- aw$cap_active
          wd <- weight_diag(b, w, aw$cap_active, cap = if (ws == "cap") 0.95 else NULL)
        }
        ## penalty convention: glmnet weights the whole penalty; "l1only"
        ## weights the L1 term only (Zou and Zhang 2009; RAENR-MM). The second
        ## is fitted for the as-implemented weights ("cap") only.
        convs <- if (ws == "cap") c("glmnet", "l1only") else "glmnet"
        for (cv_ in convs) for (lm_ in lmodes) {
          err <- ""
          cf <- tryCatch(
            if (cv_ == "glmnet") fit_adaptive(d$x, d$y, alpha, w, keep, lm_, lam_fixed)
            else fit_l1only(d$x, d$y, alpha, w, keep, lm_, lam_fixed),
            error = function(e) { err <<- conditionMessage(e); NULL })
          met <- if (is.null(cf)) NULL else
            sel_metrics(cf, d$beta, b, te$x, te$y)
          if (is.null(met)) {
            fails[[length(fails) + 1L]] <- data.frame(cell = ci, rep = r, pilot = pil,
              wscheme = ws, pen_conv = cv_, lambda_mode = lm_, error = err)
            next
          }
          k <- k + 1L
          rows[[k]] <- cbind(
            data.frame(design = DESIGN, cell = ci, p = cp$p, corr = cp$corr,
                       scen = cp$scen, rep = r, alpha = alpha, pilot = pil,
                       wscheme = ws, pen_conv = cv_, lambda_mode = lm_,
                       lam_fixed = lam_fixed,
                       lam_used = nz(attr(cf, "lam_used")),
                       lam_at_edge = nz(attr(cf, "lam_at_edge"))),
            wd[, c("Zsize", "Zfrac", "cap_active", "cap_inert_pred",
                   "maxwZ", "bound", "log10_ratio")],
            met)
        }
      }
    }
  }
  res <- if (length(rows)) do.call(rbind, rows) else data.frame()
  attr(res, "failures") <- if (length(fails)) do.call(rbind, fails) else NULL
  saveRDS(res, ofile)
  cat(sprintf("cell %3d  p=%4d %-5s %-6s  %d rows, %d failed fits  %.1f min\n", ci, cp$p, cp$corr,
              cp$scen, nrow(res), length(fails),
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
cat("Shard", SHARD, "finished.\n")
