# Sparse pilots bound adaptive lasso selection

R code and results for the paper

> Alabi, N. O. and Oyeyemi, G. M. *Sparse Pilots Bound Adaptive Lasso Selection: Exact Weights and Software Conventions* (submitted).

When the pilot estimate of an adaptive lasso or adaptive elastic net has many exact zeros, and the penalty level is chosen by cross validation on glmnet's default path, the adaptive step cannot bring back a variable the pilot left out. The function `check_pilot_bound()` in `R/check_pilot_bound.R` tells you, for your own data and pilot, whether this applies.

## Quick start: check your own analysis

```r
install.packages("glmnet")
source("R/check_pilot_bound.R")
library(glmnet)
pilot <- as.vector(coef(cv.glmnet(x, y), s = "lambda.min"))[-1]
check_pilot_bound(x, y, pilot)
```

## Reproducing the paper

Everything is written in R. The only package needed is glmnet (version 5.1 was used, under R 4.3).

## How to reproduce

1. Open R and set this folder (the one that contains `R/`, `results/` and `data/`) as the working directory.
2. Run `source("run_all.R")`.

By default `run_all.R` uses the saved results in `results/`, which are the ones reported in the paper. In a few minutes it regenerates:

- the summary tables (`results/summary/`);
- Figures 1 and 2 (`figures/`);
- the numerical checks of the Appendix (`results/numerical_checks.txt`).

To rerun the simulations as well, set `RERUN_SIMULATIONS <- TRUE` at the top of `run_all.R`. This overwrites `results/` and takes about a day on one core. Each experiment script can also be split into shards and run in parallel with `Rscript`; the usage line is at the top of each script. Every script sets its seeds, so a rerun with the same versions of R and glmnet should reproduce the saved results.

## Scripts (`R/`)

| Script | Experiment in the paper | Output |
|---|---|---|
| `common.R` | Shared code: data generation, pilots, the adaptive weights as implemented (cap and rescaling), and the solvers for conventions W and L | none |
| `E1_weights.R` | E1, weight structure (P1) | `results/E1/` |
| `E2_factorial.R` | E2, main factorial (P2 to P5); argument `main`, or `lasso` for alpha = 1 | `results/E2v4_main_R200/`, `results/E2v4_lasso_R200/` |
| `E3_robust.R` | E3, robust adaptive elastic net (P6) | `results/E3v4_R100/` |
| `E4_riboflavin.R` | E4, riboflavin data (P7) | `results/E4/` |
| `E5_repairs.R` | E5, post hoc, Parts A and B | `results/E5v2_R100/` |
| `E5c_tolerance.R` | E5c, post hoc, tight tolerance check of E5 Part B | `results/E5c_R100/` |
| `E6_wider.R` | E6, post hoc, wider design and direct checks of the theory | `results/E6v2_R100/` |
| `E7_robust_scaled.R` | E7, post hoc, the robust estimator on the observed and the standardised response | `results/E7v2_R100/` |
| `check_ystd.R` | Check of glmnet's internal standardisation of y (Section 2.5) | printed |
| `combine.R` | Summary tables | `results/summary/` |
| `make_figures.R` | Figures 1 and 2 | `figures/` |
| `numerical_checks.R` | Numerical checks of Propositions 1, 3 and 4 (Appendix) | `results/numerical_checks.txt` |
| `check_pilot_bound.R` | `check_pilot_bound(x, y, pilot)`: for your own data and pilot, the implemented weights, the effective zero set and the condition of Proposition 4 (whether any variable the pilot left out can be selected) | printed |
| `test_check_pilot_bound.R` | Checks `check_pilot_bound()` against glmnet on 50 simulated data sets | printed |

### Checking your own analysis

```r
source("R/check_pilot_bound.R")
library(glmnet)
pilot <- as.vector(coef(cv.glmnet(x, y), s = "lambda.min"))[-1]
check_pilot_bound(x, y, pilot)
```

## Data

`data/riboflavin.RData` holds the riboflavin data (n = 71, p = 4088), taken unchanged from the R package hdi, version 0.1-10, on CRAN. If the file is missing, `E4_riboflavin.R` loads the data from an installed hdi package, or downloads the hdi source archive and extracts the data from it.

## Study plan

Supplementary File S1 reproduces the study plan, with its predictions P1 to P7 and amendments A1 to A6. Experiments E5 to E7 were designed after the planned runs and are labelled post hoc in the paper.

## Licence

The code is released under the GNU General Public License, version 3 (see `LICENSE`). The riboflavin data in `data/` come unchanged from the R package hdi (version 0.1-10, CRAN), which is distributed under the GPL.

## Contact

Nurudeen Olawale Alabi, Department of Mathematics and Statistics, Federal Polytechnic Ilaro, Nigeria (nurudeen.alabi@federalpolyilaro.edu.ng).
