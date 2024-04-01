library(chkptstanr)
library(cmdstanr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

source("functions.R")

args = commandArgs(trailingOnly = TRUE)

set_cmdstan_path("cmdstan-2.34.1")

stan_code = read_lines("solanum-aa4.stan") |> 
  paste(collapse = "\n")

chkpt_path = paste0("chkpt_folder_aa4_", args[1])

stan_rh_curves = read_rds("stan_rh_curves.rds")
init = read_rds("init.rds")

fit_m = chkpt_stan1(
  model_code = stan_code,
  data = stan_rh_curves,
  iter_warmup = 4e3,
  iter_sampling = 4e3,
  iter_typical = 150,
  thin = 4e0,
  iter_per_chkpt = 1e1,
  exit_after = NULL,
  chkpt_progress = TRUE,
  path = chkpt_path,
  init = init,
  max_treedepth = 12L
)

if (length(list.files(paste0(chkpt_path, "/cp_samples")) > 0)) {
  draws = combine_chkpt_draws1(path = chkpt_path)
  write_rds(draws, paste0("draws_aa4_", args[1], ".rds"))
}

tar(paste0(chkpt_path, ".tar"), chkpt_path)
