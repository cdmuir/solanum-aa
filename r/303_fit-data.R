# Fit model to actual data
# This version runs locally; I have other code for running on HPC
source("r/header.R")

m = "aa1"

stan_code = read_lines(glue("stan/solanum-{m}.stan")) |> 
  paste(collapse = "\n")

chkpt_path = glue("checkpoints/chkpt_folder_{m}")

stan_rh_curves = read_rds("data/stan_rh_curves.rds")
init = read_rds("objects/init.rds")

fit_m = chkpt_stan1(
  model_code = stan_code,
  data = stan_rh_curves,
  iter_warmup = 4e3,
  iter_sampling = 4e3,
  thin = 4e0,
  iter_per_chkpt = 2e2,
  chkpt_progress = TRUE,
  path = chkpt_path,
  init = init,
  max_treedepth = 12L
)

draws = combine_chkpt_draws1(path = chkpt_path)

write_rds(draws, glue("objects/draws_{m}.rds"))
tar(paste0(chkpt_path, ".tar"), chkpt_path)
