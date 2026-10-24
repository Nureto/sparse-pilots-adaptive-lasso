## =====================================================================
##  combine.R : collect E1-E3 results and write summary tables
##  Usage:  Rscript R/combine.R           (after the runs have finished)
##  Output: results/summary/*.csv  -- zip the whole results folder and send it.
## =====================================================================
dir.create("results/summary", showWarnings = FALSE, recursive = TRUE)
mcse <- function(v) { v <- v[is.finite(v)]; if (length(v) < 2) NA else sd(v) / sqrt(length(v)) }
agg <- function(df, by, vars) {
  sp <- split(df, df[by], drop = TRUE)
  do.call(rbind, lapply(sp, function(g) {
    out <- g[1, by, drop = FALSE]; out$n_fits <- nrow(g)
    for (v in vars) {
      x <- g[[v]]
      out[[paste0(v, "_mean")]] <- mean(x[is.finite(x)])
      out[[paste0(v, "_mcse")]] <- mcse(x)
    }
    out
  }))
}
## among files tagged _R<n> or _B<n>, take the one with the largest n
largest <- function(f) {
  if (!length(f)) return(f)
  n <- as.numeric(sub(".*_[RB](\\d+)\\.rds$", "\\1", basename(f)))
  f[which.max(n)]
}
read_dir <- function(d) {
  f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(f)) return(NULL)
  do.call(rbind, lapply(f, readRDS))
}

## ---- E1 ------------------------------------------------------------------
e1f <- list.files("results/E1", pattern = "^E1_R\\d+\\.rds$", full.names = TRUE)
if (length(e1f)) {
  e1 <- readRDS(largest(e1f))
  e1$P1a_match <- e1$cap_active == !e1$cap_inert_pred
  e1$ratio_to_bound <- e1$maxwZ / e1$bound
  ## A2 identity (cap inert): min_Z w / max_off w = 1 + min|b_off|^2 / delta
  if ("minabs_off" %in% names(e1)) {
    ok <- !e1$cap_active & is.finite(e1$log10_ratio)
    pred <- log10(1 + e1$minabs_off^2 / e1$delta)
    cat(sprintf("A2 identity: max abs error %.2e over %d fits\n",
                max(abs(e1$log10_ratio[ok] - pred[ok])), sum(ok)))
  }
  t1 <- agg(e1, c("pilot", "p", "delta"),
            c("Zfrac", "P1a_match", "cap_active", "ratio_to_bound", "log10_ratio"))
  write.csv(t1, "results/summary/E1_weights.csv", row.names = FALSE)
  cat(sprintf("E1: P1a holds in %.2f%% of %d fits\n", 100 * mean(e1$P1a_match), nrow(e1)))
}

## ---- E2 ------------------------------------------------------------------
vars2 <- c("nsel", "tpr", "fpr", "prec", "f1", "mcc", "jac_pilot", "readmit",
           "dropped", "nsel_minus_pmZ", "coef_err", "rmse_test", "Zfrac")
for (d in list.dirs("results", recursive = FALSE)) {
  if (!grepl("^E2", basename(d))) next
  e2 <- read_dir(d); if (is.null(e2)) next
  if (!"pen_conv" %in% names(e2)) e2$pen_conv <- "glmnet"   # first smoke runs
  tag <- basename(d)
  write.csv(agg(e2, c("pen_conv", "pilot", "wscheme", "lambda_mode", "p"), vars2),
            sprintf("results/summary/%s_by_p.csv", tag), row.names = FALSE)
  write.csv(agg(e2, c("pen_conv", "pilot", "wscheme", "lambda_mode", "corr", "scen"), vars2),
            sprintf("results/summary/%s_by_design.csv", tag), row.names = FALSE)
  ## P2: n_selected against p - |Z| for sparse pilots
  reg <- list(); k <- 0L
  for (cv_ in unique(e2$pen_conv)) for (pil in c("lasso", "huber")) for (lm_ in c("fixed", "retuned")) {
    g <- e2[e2$pen_conv == cv_ & e2$pilot == pil & e2$wscheme == "cap" & e2$lambda_mode == lm_, ]
    if (nrow(g) < 10) next
    pmz <- g$pilot_nsel
    fit <- lm(g$nsel ~ pmz); ci <- confint(fit)
    k <- k + 1L
    reg[[k]] <- data.frame(pen_conv = cv_, pilot = pil, lambda_mode = lm_, n = nrow(g),
      intercept = coef(fit)[1], slope = coef(fit)[2],
      slope_lo = ci[2, 1], slope_hi = ci[2, 2], r2 = summary(fit)$r.squared)
  }
  if (k) write.csv(do.call(rbind, reg), sprintf("results/summary/%s_P2_regression.csv", tag),
                   row.names = FALSE)
  ## P4: finite delta vs theory, same pilot and lambda mode, same replication
  key <- c("cell", "rep", "pilot", "lambda_mode")
  a <- e2[e2$wscheme == "cap" & e2$pen_conv == "glmnet", c(key, "nsel", "tp", "fp")]
  b <- e2[e2$wscheme == "theory", c(key, "nsel", "tp", "fp")]
  m <- merge(a, b, by = key, suffixes = c("_cap", "_theory"))
  m$same <- m$nsel_cap == m$nsel_theory & m$tp_cap == m$tp_theory & m$fp_cap == m$fp_theory
  write.csv(agg(m, c("pilot", "lambda_mode"), "same"),
            sprintf("results/summary/%s_P4_finite_vs_theory.csv", tag), row.names = FALSE)
  fl <- lapply(list.files(d, "\\.rds$", full.names = TRUE), function(f) attr(readRDS(f), "failures"))
  fl <- do.call(rbind, fl)
  cat(sprintf("%s: %d fits, %d failed fits\n", tag, nrow(e2), if (is.null(fl)) 0 else nrow(fl)))
}

## ---- E3 ------------------------------------------------------------------
for (d in list.dirs("results", recursive = FALSE)) {
  if (!grepl("^E3", basename(d))) next
  e3 <- read_dir(d); if (is.null(e3)) next
  write.csv(agg(e3, c("pilot", "lambda_mode", "p", "contam"), vars2),
            sprintf("results/summary/%s.csv", basename(d)), row.names = FALSE)
  cat(sprintf("%s: %d fits\n", basename(d), nrow(e3)))
}

## ---- E4 ------------------------------------------------------------------
if (file.exists("results/E4/E4_full.rds")) {
  write.csv(readRDS("results/E4/E4_full.rds")$table, "results/summary/E4_full.csv", row.names = FALSE)
  sf <- list.files("results/E4", "^E4_stability", full.names = TRUE)
  if (length(sf)) write.csv(readRDS(largest(sf)), "results/summary/E4_stability.csv", row.names = FALSE)
}
cat("Summaries written to results/summary/\n")
