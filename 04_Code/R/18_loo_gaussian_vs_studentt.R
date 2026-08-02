## =============================================================================
## 18_loo_gaussian_vs_studentt.R
##
## Approximate leave-one-out comparison of the Gaussian and Student-t
## observation models for the primary (independent random-effects) joint model.
##
## WHY. The Student-t refit drives nu to its lower bound (2.20), which says the
## likelihood prefers heavy tails. That is not the same as the t model being
## better: heavy tails buy pointwise fit by letting the model follow individual
## observations, and the penalty for that flexibility may cancel the gain. A
## rough WAIC calculation on these draws gave an elpd difference of about -4
## with a standard error near 10 -- i.e. no material difference -- but with
## p_waic near 119 on 495 observations, WAIC is not trustworthy here. PSIS-LOO
## with Pareto-k diagnostics is the right tool, and this script runs it.
##
## The manuscript states the qualitative conclusion (the two specifications are
## not distinguishable) without quoting a figure. INSERT THE NUMBERS BELOW into
## Section "Sensitivity to the censoring point and to serial correlation" once
## this has been run, and check the Pareto-k column before doing so.
##
## Only the LONGITUDINAL contributions are compared: log_lik_y is the
## observation-level log-likelihood, and the event sub-model is identical in the
## two fits, so this isolates the observation model.
##
## RUN from the project root:
##   Rscript 04_Code\R\18_loo_gaussian_vs_studentt.R
## REQUIREMENTS: loo, posterior, cmdstanr (or just the CSVs and loo/posterior).
## =============================================================================
suppressMessages({ library(posterior); library(cmdstanr); library(loo) })

root <- normalizePath(".", winslash = "/")
draws_dir <- file.path(root, "08_Model", "stan_joint_interval_dmr_independent_draws")
if (!dir.exists(draws_dir)) stop("Run from the project root: ", draws_dir, " not found.")

## NB. The Gaussian model's file prefix is a strict prefix of the Student-t
## model's, so the Student-t files must be excluded BEFORE selecting the most
## recent run -- otherwise the newest "Gaussian" match is a Student-t file and
## the subsequent filter leaves nothing.
pick <- function(prefix, exclude_studentt) {
  f <- list.files(draws_dir, pattern = paste0("^", prefix, ".*\\.csv$"), full.names = TRUE)
  if (exclude_studentt) f <- f[!grepl("studentt", basename(f))]
  if (!length(f)) stop("No draws found for ", prefix)
  tag <- sub("^.*-([0-9]{12})-[0-9]+-[a-f0-9]+\\.csv$", "\\1", basename(f))
  f[tag == max(tag)]
}
f_gauss <- pick("glw_joint_interval_dmr_independent",           TRUE)
f_t     <- pick("glw_joint_interval_dmr_independent_studentt",  FALSE)
message("Gaussian chains: ", length(f_gauss), " | Student-t chains: ", length(f_t))
stopifnot(length(f_gauss) > 0, length(f_t) > 0)

## read_cmdstan_csv lives in cmdstanr, not posterior. post_warmup_draws is a
## draws_array with dimensions iteration x chain x variable, which is exactly
## what relative_eff() wants; loo() itself takes the matrix form.
ll_array <- function(files)
  cmdstanr::read_cmdstan_csv(files, variables = "log_lik_y")$post_warmup_draws

AG <- ll_array(f_gauss); AT <- ll_array(f_t)
LG <- posterior::as_draws_matrix(AG)
LT <- posterior::as_draws_matrix(AT)
message("log-lik matrices: Gaussian ", nrow(LG), "x", ncol(LG),
        " | Student-t ", nrow(LT), "x", ncol(LT))

## relative effective sample sizes improve the PSIS tail estimate
rG <- loo::relative_eff(exp(as.array(AG)))
rT <- loo::relative_eff(exp(as.array(AT)))

loo_G <- loo::loo(LG, r_eff = rG)
loo_T <- loo::loo(LT, r_eff = rT)

cat("\n===== Gaussian =====\n");  print(loo_G)
cat("\n===== Student-t =====\n"); print(loo_T)

cat("\n===== comparison (positive elpd_diff favours the FIRST row) =====\n")
cmp <- loo::loo_compare(list(Gaussian = loo_G, StudentT = loo_T))
print(cmp)

## Pareto-k: values above 0.7 mean the PSIS approximation is unreliable for
## those observations and the comparison should be treated with caution.
kG <- loo::pareto_k_values(loo_G); kT <- loo::pareto_k_values(loo_T)
cat(sprintf("\nPareto-k > 0.7:  Gaussian %d/%d,  Student-t %d/%d\n",
            sum(kG > 0.7), length(kG), sum(kT > 0.7), length(kT)))
cat(sprintf("p_loo:           Gaussian %.1f,  Student-t %.1f  (n = %d observations)\n",
            loo_G$estimates["p_loo","Estimate"], loo_T$estimates["p_loo","Estimate"], ncol(LG)))

out <- file.path(root, "08_Model", "loo_gaussian_vs_studentt.csv")
write.csv(data.frame(
  model = c("Gaussian","StudentT"),
  elpd_loo = c(loo_G$estimates["elpd_loo","Estimate"], loo_T$estimates["elpd_loo","Estimate"]),
  se_elpd  = c(loo_G$estimates["elpd_loo","SE"],       loo_T$estimates["elpd_loo","SE"]),
  p_loo    = c(loo_G$estimates["p_loo","Estimate"],    loo_T$estimates["p_loo","Estimate"]),
  k_gt_0.7 = c(sum(kG > 0.7), sum(kT > 0.7))), out, row.names = FALSE)
message("\nWrote ", out)
message("Insert elpd_diff and its SE into Section sec:sens-extra, ",
        "and report the Pareto-k counts if any exceed 0.7.")
