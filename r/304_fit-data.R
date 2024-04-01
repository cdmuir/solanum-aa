# Fit model to actual data
# This version runs locally; I have other code for running on HPC
source("r/header.R")

m = "aa2"

chkpt_path = glue("checkpoints/chkpt_folder_{m}")

# Only untar if path doesn't exist
if (!dir.exists(chkpt_path)) {
  untar(paste0(chkpt_path, ".tar"))
  file.remove(paste0(chkpt_path, ".tar"))
}

dir.create("checkpoints/another_name")
untar(paste0(chkpt_path, ".tar"), "checkpoints/another_name")

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

write_rds(draws, glue("objects/draws_{m}.rds"))
tar(paste0(chkpt_path, ".tar"), chkpt_path)
