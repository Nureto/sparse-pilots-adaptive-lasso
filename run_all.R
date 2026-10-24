## =====================================================================
##  run_all.R : reproduce every result, table, figure and check in
##  "Sparse pilots bound adaptive lasso selection: exact weights and
##   software conventions"
##
##  Set the working directory to this folder (the one that contains R/,
##  results/ and data/), then run this file. Requires R >= 4.3 and the
##  glmnet package (version 5.1 was used): install.packages("glmnet").
##
##  The saved results in results/ are those reported in the paper, so the
##  figures and the numerical checks (last block) can be reproduced in a few
##  minutes without rerunning the simulations. Rerunning all experiments
##  overwrites results/ and takes about a day on one core; E2 dominates.
##  Each experiment script can also be split into shards and run in
##  parallel with Rscript: see the usage line at the top of each script.
## =====================================================================

RERUN_SIMULATIONS <- FALSE      # set to TRUE to rerun E1 to E7

if (RERUN_SIMULATIONS) {
  ARGS <- c("200");                    source("R/E1_weights.R")       # E1, weights
  ARGS <- c("1", "1", "200", "main");  source("R/E2_factorial.R")     # E2, main factorial
  ARGS <- c("1", "1", "200", "lasso"); source("R/E2_factorial.R")     # E2, alpha = 1
  ARGS <- c("1", "1", "100");          source("R/E3_robust.R")        # E3, robust estimator
  ARGS <- c("100");                    source("R/E4_riboflavin.R")    # E4, riboflavin data
  ARGS <- c("1", "1", "100");          source("R/E5_repairs.R")       # E5, post hoc
  ARGS <- c("100");                    source("R/E5c_tolerance.R")    # E5c, post hoc
  ARGS <- c("1", "1", "100");          source("R/E6_wider.R")         # E6, post hoc
  ARGS <- c("100");                    source("R/E7_robust_scaled.R") # E7, post hoc
  rm(ARGS)
}

source("R/check_ystd.R")        # glmnet's internal scaling of y (Section 2.5)
source("R/combine.R")           # summary tables in results/summary/
source("R/make_figures.R")      # Figures 1 and 2 in figures/
source("R/numerical_checks.R")  # Appendix checks, results/numerical_checks.txt
