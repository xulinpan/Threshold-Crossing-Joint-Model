## 00_setup.R --------------------------------------------------------------
## Global configuration, packages, palettes, and paths for the Bayesian
## hierarchical threshold-crossing model pipeline.
## -------------------------------------------------------------------------

## ---- packages ----
## cmdstanr is preferred; rstan is a fallback (see 02_fit_models.R).
suppressWarnings(suppressMessages({
  have <- function(p) requireNamespace(p, quietly = TRUE)
  need <- c("posterior", "RColorBrewer", "ggsci")
  for (p in need) if (!have(p)) message("NOTE: package '", p, "' not installed.")
  if (have("cmdstanr")) library(cmdstanr)
  if (have("posterior")) library(posterior)
  library(RColorBrewer)
}))

## ---- reproducibility ----
set.seed(20260709)
options(mc.cores = max(1L, parallel::detectCores() - 1L))

## ---- project paths (edit ROOT if running elsewhere) ----
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
get_script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))
  getwd()
}
CODE_DIR <- get_script_dir()                          # .../04_Code/R
ROOT     <- normalizePath(file.path(CODE_DIR, ".."), mustWork = FALSE)  # .../04_Code
if (!dir.exists(file.path(ROOT, "stan"))) ROOT <- ".."  # interactive fallback
## All outputs are self-contained under this model folder to avoid touching
## the manuscript's real figures/tables. Copy selected figures out manually.
DIR_STAN  <- file.path(ROOT, "stan")
DIR_OUT   <- file.path(ROOT, "outputs");           dir.create(DIR_OUT, showWarnings = FALSE, recursive = TRUE)
DIR_FIG   <- file.path(DIR_OUT, "figures");         dir.create(DIR_FIG, showWarnings = FALSE, recursive = TRUE)
DIR_TAB   <- file.path(DIR_OUT, "tables");          dir.create(DIR_TAB, showWarnings = FALSE, recursive = TRUE)

## ---- model constants ----
C_F   <- -5.0     # assay floor
C_D   <- -4.5     # DMR threshold
GRID  <- 60L      # grid points for soft-min running minimum
KAPPA <- 50       # soft-min sharpness

## ---- ground-truth parameters (posterior means of the fitted CML model) ----
TRUTH <- list(
  beta0 = -2.074, beta1 = -3.561, beta2 = 0.502, beta_bm = 0.607,
  sigma_y = 1.813, tau0 = 0.968, tau1 = 2.203, sigma_thr = 0.20,
  gamma0 = -4.095, gamma1 = -0.228, gamma2 = -1.831, alpha = -1.275
)

## ---- publication palette (RColorBrewer, colour-blind aware) ----
PAL_SET1 <- ggsci::pal_nejm("default")(8)   # ggsci NEJM (categorical)
PAL_DARK <- ggsci::pal_lancet("lanonc")(8)  # ggsci Lancet (categorical)
PAL_BLUE <- brewer.pal(9, "Blues")    # sequential
COL <- list(
  floor   = PAL_DARK[1],
  exact   = PAL_DARK[2],
  thresh  = PAL_DARK[3],
  interval= PAL_DARK[4],
  obs     = PAL_SET1[2],
  cens    = PAL_SET1[1],
  dmr     = PAL_SET1[1],
  cmr     = PAL_SET1[3],
  trend   = PAL_SET1[5]
)

## ---- helper: open a 600-dpi PNG device (ragg if available, else grDevices) ----
open_png <- function(file, width = 7, height = 5, res = 600) {
  if (requireNamespace("ragg", quietly = TRUE)) {
    ragg::agg_png(file, width = width, height = height, units = "in", res = res)
  } else {
    grDevices::png(file, width = width, height = height, units = "in", res = res, type = "cairo")
  }
}

message("Setup complete. Stan dir: ", DIR_STAN, " | figures -> ", DIR_FIG)
