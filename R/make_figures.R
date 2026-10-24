## =====================================================================
##  make_figures.R : Figures 1 and 2 of the paper, from the saved results
##  Usage:  Rscript R/make_figures.R        (run from PilotWeights_Q)
##  Needs:  results/E2v4_main_R200 and results/E5c_R100. Base R only.
##  Writes: figures/Figure1.{pdf,eps,tif,png} and figures/Figure2.{...},
##          and, in the authors' working folder, copies them into the
##          manuscript folders.
## =====================================================================

read_dir <- function(d) {
  f <- list.files(d, pattern = "\\.rds$", full.names = TRUE)
  if (!length(f)) stop("no results found in ", d)
  do.call(rbind, lapply(f, readRDS))
}
E2 <- read_dir("results/E2v4_main_R200")
C  <- read_dir("results/E5c_R100")
cat(sprintf("E2: %d fits; E5c: %d fits\n", nrow(E2), nrow(C)))

dir.create("figures", showWarnings = FALSE)
col_lasso <- "#1f5aa6"; col_huber <- "#e07b00"
col_W     <- "#3a7d44"; col_L     <- "#7b4fa0"
alpha_col <- function(col, a) grDevices::adjustcolor(col, alpha.f = a)

## draw once on every device: vector PDF and EPS (fonts embedded by cairo),
## and 600 dpi TIFF and 200 dpi PNG for previews
save_all <- function(name, width, height, draw) {
  base <- file.path("figures", name)
  grDevices::cairo_pdf(paste0(base, ".pdf"), width, height, family = "sans")
  draw(); grDevices::dev.off()
  grDevices::cairo_ps(paste0(base, ".eps"), width, height, family = "sans",
                      onefile = FALSE, fallback_resolution = 600)
  draw(); grDevices::dev.off()
  grDevices::tiff(paste0(base, ".tif"), width, height, units = "in", res = 600,
                  compression = "lzw", type = "cairo", family = "sans")
  draw(); grDevices::dev.off()
  grDevices::png(paste0(base, ".png"), width, height, units = "in", res = 200,
                 type = "cairo", family = "sans")
  draw(); grDevices::dev.off()
  cat("written", base, "(.pdf .eps .tif .png)\n")
}

## ---------------------------------------------------------------------
##  Figure 1. E2, convention W, nonempty pilots: |S-hat| against |S|,
##  3,000 fits per pilot and panel, jittered
## ---------------------------------------------------------------------
d1 <- E2[E2$wscheme == "cap" & E2$pilot_nsel > 0 & E2$pen_conv == "glmnet", ]
set.seed(1)
pts <- list()
for (lm in c("fixed", "retuned")) for (pil in c("lasso", "huber")) {
  z <- d1[d1$lambda_mode == lm & d1$pilot == pil, ]
  z <- z[sample.int(nrow(z), 3000L), ]
  pts[[paste(lm, pil)]] <- data.frame(
    x = z$pilot_nsel + runif(nrow(z), -0.3, 0.3),
    y = z$nsel       + runif(nrow(z), -0.3, 0.3))
}
draw_fig1 <- function() {
  op <- par(mfrow = c(1, 2), mar = c(3.1, 3.3, 1.6, 0.6), mgp = c(2, 0.6, 0),
            cex = 0.75, las = 1, tcl = -0.3)
  on.exit(par(op))
  titles <- list(fixed   = expression("(a) Fixed " * lambda),
                 retuned = expression("(b) Retuned " * lambda))
  for (lm in c("fixed", "retuned")) {
    plot(NA, xlim = c(0, 100), ylim = c(0, 100), xaxs = "i", yaxs = "i",
         xlab = expression("Pilot support size |" * italic(S) * "|"),
         ylab = if (lm == "fixed") expression("Selected set size |" * hat(italic(S)) * "|") else "",
         yaxt = if (lm == "fixed") "s" else "n")
    if (lm != "fixed") axis(2, labels = FALSE)
    title(main = titles[[lm]], font.main = 1, cex.main = 1.05, adj = 0.5)
    points(pts[[paste(lm, "lasso")]], pch = 16, cex = 0.25, col = alpha_col(col_lasso, 0.3))
    points(pts[[paste(lm, "huber")]], pch = 17, cex = 0.25, col = alpha_col(col_huber, 0.3))
    abline(0, 1, lty = 2, lwd = 0.8)
    abline(h = 10, lty = 3, lwd = 0.7, col = "grey45")
    if (lm == "fixed")
      legend("topleft", c("Lasso pilot", "Huber pilot"), pch = c(16, 17),
             col = c(col_lasso, col_huber), bty = "n", pt.cex = 1)
  }
}
save_all("Figure1", 6.3, 2.9, draw_fig1)

## ---------------------------------------------------------------------
##  Figure 2. Pilot zeros selected at a fixed lambda, conventions W and L,
##  nonempty pilots: (a) E2 in the cells shared with E5; (b) E5c
## ---------------------------------------------------------------------
se <- function(v) sd(v) / sqrt(length(v))
E2s <- E2[E2$p %in% c(100, 500) & E2$corr %in% c("equi", "ar1") &
          E2$scen %in% c("normal", "t3_y10") & E2$wscheme == "cap" &
          E2$lambda_mode == "fixed" & E2$pilot_nsel > 0, ]
Cs <- C[C$alpha == 0.2 & C$pilot_nsel > 0, ]
tab <- NULL
for (pil in c("lasso", "huber")) for (cv in c("W", "L")) {
  conv_name <- if (cv == "W") "glmnet" else "l1only"
  za <- E2s$readmit[E2s$pilot == pil & E2s$pen_conv == conv_name]
  zb <- Cs[[if (cv == "W") "readmit_W" else "readmit_L"]][Cs$pilot == pil]
  tab <- rbind(tab,
    data.frame(panel = "a", pilot = pil, conv = cv, n = length(za), mean = mean(za), se = se(za)),
    data.frame(panel = "b", pilot = pil, conv = cv, n = length(zb), mean = mean(zb), se = se(zb)))
}
tab <- tab[order(tab$panel, tab$pilot != "lasso", tab$conv != "W"), ]
cat("\nFigure 2 values (mean and Monte Carlo standard error):\n")
print(format(tab, digits = 3), row.names = FALSE)
write.csv(tab, "figures/Figure2_values.csv", row.names = FALSE)

draw_fig2 <- function() {
  op <- par(mfrow = c(1, 2), mar = c(2.4, 3.6, 1.6, 0.6), mgp = c(2.4, 0.6, 0),
            cex = 0.75, las = 1, tcl = -0.3)
  on.exit(par(op))
  titles <- c(a = "(a) E2, response as observed", b = "(b) E5c, standardised response")
  for (pn in c("a", "b")) {
    t <- tab[tab$panel == pn, ]
    m  <- rbind(W = t$mean[t$conv == "W"], L = t$mean[t$conv == "L"])
    s  <- rbind(W = t$se[t$conv == "W"],   L = t$se[t$conv == "L"])
    top <- max(m + 1.96 * s) * 1.08
    x <- barplot(m, beside = TRUE, ylim = c(0, top), col = c(col_W, col_L),
                 border = "white", names.arg = c("Lasso pilot", "Huber pilot"),
                 ylab = "Pilot zeros selected (mean)", space = c(0.15, 0.8))
    ## hatching on the convention L bars
    rect(x[2, ] - 0.5, 0, x[2, ] + 0.5, m[2, ], density = 18, angle = 45,
         col = "white", border = NA, lwd = 0.8)
    arrows(x, m - 1.96 * s, x, m + 1.96 * s, angle = 90, code = 3,
           length = 0.03, lwd = 0.8)
    title(main = titles[[pn]], font.main = 1, cex.main = 1.05)
    box(bty = "l")
    if (pn == "a")
      legend("topright", c("Convention W", "Convention L"), fill = c(col_W, col_L),
             border = "white", bty = "n")
  }
}
save_all("Figure2", 6.3, 2.7, draw_fig2)

## ---------------------------------------------------------------------
##  Copy into the manuscript folders (authors' working copy only)
## ---------------------------------------------------------------------
if (dir.exists("manuscript_csam/figures")) {
for (nm in c("Figure1", "Figure2")) {
  for (ext in c("pdf", "eps", "tif", "png"))
    file.copy(file.path("figures", paste0(nm, ".", ext)),
              file.path("manuscript_csam/figures", paste0(nm, ".", ext)), overwrite = TRUE)
  for (ext in c("pdf", "eps"))
    file.copy(file.path("figures", paste0(nm, ".", ext)),
              file.path("CSAM_submission/figures", paste0(nm, ".", ext)), overwrite = TRUE)
}
for (ext in c("pdf", "png")) {
  file.copy(file.path("figures", paste0("Figure1.", ext)),
            file.path("manuscript/figures", paste0("fig1_sel_vs_pilot.", ext)), overwrite = TRUE)
  file.copy(file.path("figures", paste0("Figure2.", ext)),
            file.path("manuscript/figures", paste0("fig2_conventions.", ext)), overwrite = TRUE)
}
cat("\nFigures copied to manuscript_csam/figures, CSAM_submission/figures and manuscript/figures\n")
}
