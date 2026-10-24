## =====================================================================
##  E4: real data, riboflavin (n = 71, p = 4088)             (prediction P7)
##  Usage:  Rscript R/E4_riboflavin.R [B]        B = bootstrap resamples
##  Smoke:  Rscript R/E4_riboflavin.R 3
##  Data:   riboflavin from hdi 0.1-10 (see load_riboflavin below);
##          the hdi package itself does not need to be installed.
## =====================================================================
source("R/common.R")
## The riboflavin data (n = 71, p = 4088) ship with the hdi package, version
## 0.1-10. Ready-built hdi is only published for R 4.5 and later, so rather
## than installing the package this script takes the data file straight out
## of hdi's source archive, once, and keeps a copy in data/riboflavin.RData.
## Works on any R version; needs internet only the first time.
load_riboflavin <- function(cache = "data/riboflavin.RData") {
  if (requireNamespace("hdi", quietly = TRUE)) {
    e <- new.env(); utils::data("riboflavin", package = "hdi", envir = e)
    return(e$riboflavin)
  }
  if (!file.exists(cache)) {
    dir.create(dirname(cache), showWarnings = FALSE, recursive = TRUE)
    urls <- c("https://cloud.r-project.org/src/contrib/hdi_0.1-10.tar.gz",
              "https://cloud.r-project.org/src/contrib/Archive/hdi/hdi_0.1-10.tar.gz",
              "https://meierluk.r-universe.dev/src/contrib/hdi_0.1-10.tar.gz")
    tgz <- tempfile(fileext = ".tar.gz"); got <- FALSE
    for (u in urls) {
      message("Downloading the riboflavin data from ", u)
      ok <- tryCatch(utils::download.file(u, tgz, mode = "wb", quiet = TRUE) == 0,
                     error = function(e) FALSE, warning = function(w) FALSE)
      if (ok && file.exists(tgz) && file.size(tgz) > 1e5) { got <- TRUE; break }
    }
    if (!got) stop("Could not download hdi_0.1-10.tar.gz. Check the internet ",
                   "connection, or install hdi 0.1-10 from CRAN.")
    member <- grep("^hdi/data/riboflavin\\.(RData|rda)$", utils::untar(tgz, list = TRUE),
                   value = TRUE)
    if (!length(member)) stop("riboflavin data file not found inside the archive.")
    td <- tempfile(); utils::untar(tgz, files = member[1], exdir = td)
    file.copy(file.path(td, member[1]), cache, overwrite = TRUE)
    message("Saved a copy of the data in ", cache)
  }
  e <- new.env(); load(cache, envir = e)
  if (!exists("riboflavin", envir = e)) stop(cache, " does not contain 'riboflavin'.")
  e$riboflavin
}

B     <- as.integer(arg_or(1, "100"))
SEED  <- 20261018L
alpha <- 0.2; gamma <- 2; delta <- 1e-6

riboflavin <- load_riboflavin()
x <- as.matrix(riboflavin$x); y <- as.numeric(riboflavin$y)
storage.mode(x) <- "double"
cat(sprintf("riboflavin: n = %d, p = %d\n", nrow(x), ncol(x)))

dir.create("results/E4", showWarnings = FALSE, recursive = TRUE)
pilots <- c("lasso", "huber", "ridge")
lmodes <- c("fixed", "retuned")

one_pass <- function(x, y) {
  p <- ncol(x)
  lam_fixed <- lambda_fixed_nonadaptive(x, y, alpha)
  out <- list(); k <- 0L
  for (pil in pilots) {
    b <- run_pilot(pil, x, y)
    aw <- adaptive_weights(b, gamma, delta, cap = 0.95)
    for (lm_ in lmodes) {
      cf <- tryCatch(fit_adaptive(x, y, alpha, aw$w, NULL, lm_, lam_fixed),
                     error = function(e) NULL)
      if (is.null(cf)) next
      sel <- which(abs(cf[-1]) > 1e-8); supp <- which(b != 0)
      k <- k + 1L
      out[[k]] <- list(pilot = pil, lambda_mode = lm_, sel = sel, supp = supp,
                       Zsize = p - length(supp), cap_active = aw$cap_active)
    }
  }
  out
}

## ---- full data -----------------------------------------------------------
set.seed(SEED)
full <- one_pass(x, y)
tab <- do.call(rbind, lapply(full, function(z) data.frame(
  pilot = z$pilot, lambda_mode = z$lambda_mode, nsel = length(z$sel),
  pilot_nsel = length(z$supp), Zsize = z$Zsize, p_minus_Z = ncol(x) - z$Zsize,
  jac_pilot = jaccard(z$sel, z$supp), readmit = length(setdiff(z$sel, z$supp)),
  dropped = length(setdiff(z$supp, z$sel)), cap_active = z$cap_active)))
print(tab, row.names = FALSE)
saveRDS(list(table = tab, sets = full), "results/E4/E4_full.rds")

## ---- bootstrap stability -------------------------------------------------
boot <- vector("list", B)
t0 <- Sys.time()
for (bb in seq_len(B)) {
  set.seed(SEED + bb)
  idx <- sample.int(nrow(x), replace = TRUE)
  boot[[bb]] <- tryCatch(one_pass(x[idx, ], y[idx]), error = function(e) NULL)
  if (bb %% 10 == 0) cat(sprintf("bootstrap %d/%d  (%.1f min)\n", bb, B,
                                 as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
saveRDS(boot, sprintf("results/E4/E4_boot_B%d.rds", B))

ok <- Filter(Negate(is.null), boot)
stab <- list(); k <- 0L
for (pil in pilots) for (lm_ in lmodes) {
  pick <- function(run) Filter(function(z) z$pilot == pil && z$lambda_mode == lm_, run)
  got <- lapply(ok, pick); got <- Filter(length, got)
  if (length(got) < 2L) next
  sel_sets  <- lapply(got, function(g) g[[1]]$sel)
  supp_sets <- lapply(got, function(g) g[[1]]$supp)
  k <- k + 1L
  stab[[k]] <- data.frame(pilot = pil, lambda_mode = lm_, B_used = length(got),
    stab_selected = mean_pairwise_jaccard(sel_sets),
    stab_pilot    = mean_pairwise_jaccard(supp_sets),
    mean_jac_sel_vs_pilot = mean(mapply(jaccard, sel_sets, supp_sets)))
}
stab <- do.call(rbind, stab)
print(stab, row.names = FALSE)
saveRDS(stab, sprintf("results/E4/E4_stability_B%d.rds", B))
cat("E4 finished.\n")
