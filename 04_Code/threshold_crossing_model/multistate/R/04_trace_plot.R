## 04_trace_plot.R ---------------------------------------------------------
## MCMC trace plots (600 dpi) for the spline multi-state fit.
## Uses, in order of preference:
##   (1) the `fit` cmdstanr object still in your session (no refit needed), or
##   (2) outputs/spline_real_draws.rds saved by 03_fit_spline.R.
## Renders with bayesplot if available; otherwise a base-R fallback.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else getwd() })
suppressMessages(library(posterior))
OUT <- file.path(this_dir, "..", "outputs")
FIG <- file.path(OUT, "figures"); dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
pars <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","sigma_zeta","sigma_thr")

## ---- obtain draws ----
get_draws <- function() {
  if (exists("fit", envir = .GlobalEnv)) {
    message("Using `fit` from the current session.")
    return(as_draws_df(get("fit", .GlobalEnv)$draws(variables = pars)))
  }
  f <- file.path(OUT, "spline_real_draws.rds")
  if (file.exists(f)) { message("Using saved draws: ", f)
    return(subset_draws(readRDS(f), variable = pars)) }
  stop("No draws found. Re-run 03_fit_spline.R (it now saves spline_real_draws.rds), ",
       "or keep the `fit` object in your session before sourcing this script.")
}
dr <- get_draws()

out_png <- file.path(FIG, "figR_18_trace_spline.png")

## ---- preferred: bayesplot ----
if (requireNamespace("bayesplot", quietly = TRUE) &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  suppressMessages({ library(bayesplot); library(ggplot2) })
  bayesplot::color_scheme_set("mix-blue-red")   # distinct per-chain colors
  arr <- as_draws_array(dr)
  p <- mcmc_trace(arr, pars = pars, facet_args = list(ncol = 3)) +
       ggtitle("MCMC trace plots: spline multi-state model (real CML cohort)")
  ggsave(out_png, p, width = 11, height = 8, dpi = 600)
  message("Wrote ", out_png)
} else {
  ## ---- base-R fallback (RColorBrewer palette) ----
  suppressMessages(if (requireNamespace("RColorBrewer", quietly = TRUE)) library(RColorBrewer))
  d  <- as.data.frame(dr)
  ch <- sort(unique(d$.chain)); nch <- length(ch)
  cols <- if (requireNamespace("RColorBrewer", quietly = TRUE))
            ggsci::pal_lancet("lanonc")(max(3, nch))[seq_len(nch)]   # ggsci Lancet
          else seq_len(nch) + 1
  png(out_png, width = 11, height = 8, units = "in", res = 600, type = "cairo")
  op <- par(mfrow = c(3, 3), mar = c(3.2, 3.6, 2, 1), mgp = c(2, 0.6, 0))
  for (p_ in pars) {
    rng <- range(d[[p_]])
    plot(NA, xlim = c(1, max(d$.iteration)), ylim = rng, xlab = "iteration",
         ylab = p_, main = p_, las = 1)
    for (k in seq_along(ch)) {
      dk <- d[d$.chain == ch[k], ]
      lines(dk$.iteration, dk[[p_]], col = adjustcolor(cols[k], 0.75), lwd = 0.5)
    }
  }
  par(op); dev.off()
  message("Wrote ", out_png, " (base-R fallback)")
}
