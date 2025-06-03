# Fit model to actual data
# This version runs locally; I have other code for running on HPC
source("r/header.R")

args = commandArgs(trailingOnly = TRUE)
# args = list("aa1", "0")

m = args[1]
chain = args[2]

chkpt_path = glue("checkpoints/chkpt_folder_{m}")

# Only untar if path doesn't exist
if (!dir.exists(chkpt_path)) {
  untar(paste0(chkpt_path, ".tar"))
  options(cmdstanr_force_recompile = TRUE)
  # file.remove(paste0(chkpt_path, ".tar"))
}

chkpt_path1 = glue("checkpoints/chkpt_folder_{m}_{chain}")

if (!dir.exists(chkpt_path1)) {
  R.utils::copyDirectory(from = chkpt_path, to = chkpt_path1)
}

fit_m = chkpt_stan1(
  chkpt_path1,
  iter_warmup = 4e3,
  iter_sampling = 4e3,
  thin = 4e0,
  iter_per_chkpt = 2e2,
  chkpt_progress = TRUE,
  max_treedepth = 12L
)

draws = combine_chkpt_draws1(path = chkpt_path1)

write_rds(draws, glue("objects/draws_{m}_{chain}.rds"))
tar(paste0(chkpt_path1, ".tar"), chkpt_path1)
