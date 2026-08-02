## =====================================================================
## 17_loo_diagnostics.R  —  make the M3 vs M3s LOO comparison defensible
##
## Run in the SAME session as 16_fit_unified_model.R (uses m3, m3s).
##
## Two problems with the first LOO run:
##   (1) r_eff was not supplied, so "MCSE of elpd_loo is NA" and the SEs
##       assume independent draws -- they are optimistic under MCMC
##       autocorrelation.
##   (2) 10/495 observations had Pareto k > 0.7, so their pointwise LOO
##       terms are unreliable.
##
## This script fixes (1), characterises (2) -- in particular whether the
## bad points are the LEFT-CENSORED floor observations, which would be the
## expected and reportable explanation -- and tests whether the conclusion
## survives dropping them.
## =====================================================================

suppressPackageStartupMessages({library(loo); library(posterior)})

stopifnot(exists("m3"), exists("m3s"))   # from 16_fit_unified_model.R

ll  <- function(f) f$draws("log_lik")
## r_eff accounts for MCMC autocorrelation in the pointwise terms
reff <- function(f) loo::relative_eff(exp(ll(f)))

l3  <- loo::loo(ll(m3),  r_eff = reff(m3))
l3s <- loo::loo(ll(m3s), r_eff = reff(m3s))

cat("\n=== LOO with r_eff (M3 = no spline, M3s = spline) ===\n")
print(l3s); print(l3)
cat("\n=== comparison ===\n")
print(loo::loo_compare(list(M3 = l3, M3s = l3s)))

## ---- are the high-k points the censored floor observations? ----------
k3s   <- l3s$diagnostics$pareto_k
bad   <- which(k3s > 0.7)
isf   <- as.integer(long$log_mrd <= FLOOR)      # from the fitting script
cat(sprintf("\nhigh-k points: %d of %d\n", length(bad), length(k3s)))
cat(sprintf("  of these, %d are FLOOR (left-censored) observations (%.0f%%)\n",
            sum(isf[bad]), 100*mean(isf[bad])))
cat(sprintf("  baseline floor rate in the data: %.0f%%\n", 100*mean(isf)))
print(data.frame(obs = bad, k = round(k3s[bad], 2),
                 floor = isf[bad],
                 log_mrd = long$log_mrd[bad],
                 t_months = round(long$t_months[bad], 1)))

## ---- does the conclusion survive dropping the unreliable points? -----
e3  <- l3$pointwise[, "elpd_loo"]
e3s <- l3s$pointwise[, "elpd_loo"]
keep <- setdiff(seq_along(e3), union(bad, which(l3$diagnostics$pareto_k > 0.7)))
d_all  <- sum(e3s) - sum(e3)
d_keep <- sum(e3s[keep]) - sum(e3[keep])
se_keep <- sqrt(length(keep) * var(e3s[keep] - e3[keep]))
cat(sprintf("\nelpd_diff (M3s - M3): all points %.1f | excluding high-k %.1f (SE %.1f)\n",
            d_all, d_keep, se_keep))
cat("If the sign and rough magnitude are unchanged, the spline conclusion\n",
    "is robust to the unreliable points.\n")

## ---- gold standard if you want it: exact refit of the bad points -----
## Requires model methods compiled:
##   mod <- cmdstan_model(stan_file, force_recompile = TRUE,
##                        compile_model_methods = TRUE)
## then:  l3s_mm <- m3s$loo(variables = "log_lik", moment_match = TRUE)
