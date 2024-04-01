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

s = "chkpt_folder_aa5"
chkpt_path = paste0("checkpoints/", s)

# Only untar if path doesn't exist (first run)
if (!dir.exists(chkpt_path)) {
  untar(paste0(s, ".tar"))
  options(cmdstanr_force_recompile = TRUE)
  file.remove(paste0(s, ".tar"))
}

fit_m = chkpt_stan1(
  chkpt_path,
  iter_warmup = 40, #4e3,
  iter_sampling = 40, #4e3,
  thin = 1, #4e0,
  iter_per_chkpt = 10, #2e0,
  chkpt_progress = TRUE,
  max_treedepth = 10L #12L
)

draws = combine_chkpt_draws1(path = chkpt_path)

write_rds(draws, paste0("draws_aa5_", args[1], ".rds"))
tar(paste0(s, "_", args[1], ".tar"), chkpt_path)
