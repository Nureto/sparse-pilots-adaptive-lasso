## =====================================================================
##  E1: weight structure and delta invariance           (predictions P1a-c)
##  Usage:  Rscript R/E1_weights.R [R]          (run from PilotWeights_Q)
##  Smoke:  Rscript R/E1_weights.R 2
## =====================================================================
source("R/common.R")

R      <- as.integer(arg_or(1, "200"))
SEED   <- 20261015L
n      <- 100L
p_grid <- c(50L, 100L, 200L, 500L, 1000L)
deltas <- c(1e-8, 1e-6, 1e-4, 1e-2)
pilots <- c("lasso", "huber")

dir.create("results/E1", showWarnings = FALSE, recursive = TRUE)
out_file <- sprintf("results/E1/E1_R%d.rds", R)
rows <- list(); k <- 0L
t0 <- Sys.time()
for (pi_ in seq_along(p_grid)) {
  p <- p_grid[pi_]
  for (r in seq_len(R)) {
    set.seed(SEED + 100000L * pi_ + r)
    d <- gen_data(n, p, corr = "equi", scen = "t3")
    for (pil in pilots) {
      b <- run_pilot(pil, d$x, d$y)
      for (dl in deltas) {
        aw <- adaptive_weights(b, gamma = 2, delta = dl, cap = 0.95)
        k <- k + 1L
        rows[[k]] <- cbind(data.frame(p = p, rep = r, pilot = pil, delta = dl),
                           weight_diag(b, aw$w, aw$cap_active))
      }
    }
  }
  cat(sprintf("E1  p = %4d  done  (%.1f min elapsed)\n", p,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res <- do.call(rbind, rows)
saveRDS(res, out_file)
cat("Saved", out_file, "with", nrow(res), "rows\n")

## quick check of P1a
chk <- with(res, cap_active == !cap_inert_pred)
cat(sprintf("P1a: cap status matches corrected Prop 1(a) in %d of %d fits (%.2f%%)\n",
            sum(chk), length(chk), 100 * mean(chk)))
