
# Build Stan data lists for Bayesian joint longitudinal--interval survival model
library(readr)
library(dplyr)

build_stan_data <- function(folder = "cml_paper_datasets", prefix = "real") {
  long <- read_csv(file.path(folder, paste0(prefix, "_longitudinal_analysis.csv")), show_col_types = FALSE)
  interval <- read_csv(file.path(folder, paste0(prefix, "_interval_survival_analysis.csv")), show_col_types = FALSE)
  patient <- read_csv(file.path(folder, paste0(prefix, "_patient_level_analysis.csv")), show_col_types = FALSE)

  # Primary molecular-only analysis: K = 0
  stan_data <- list(
    N_obs = nrow(long),
    N_pat = length(unique(long$patient_num)),
    N_int = nrow(interval),
    id_obs = as.integer(long$patient_num),
    id_int = as.integer(interval$patient_num),
    y = as.numeric(scale(long$log_mrd)),
    t_obs = as.numeric(scale(long$t_months)),
    log_bcrabl = as.numeric(scale(long$log_bcrabl)),
    log_abl = as.numeric(scale(long$log_abl)),
    sample_bm = as.numeric(long$sample_bm),
    t_start = as.numeric(scale(interval$t_start)),
    t_end = as.numeric(scale(interval$t_end)),
    log_gap = as.numeric(scale(interval$log_gap)),
    event_interval = as.integer(interval$event_interval),
    K = 0,
    Z = matrix(0, nrow = length(unique(long$patient_num)), ncol = 0),
    is_event_patient = as.integer(patient$dmr_event[match(sort(unique(long$patient_num)), patient$patient_num)])
  )
  stan_data
}

# Example:
# stan_data_real <- build_stan_data("/mnt/data/cml_paper_datasets", "real")
# stan_data_sim  <- build_stan_data("/mnt/data/cml_paper_datasets", "simulated")
