# Initialize fit on HTC
source("r/header.R")

m = "aa1"

stan_code = read_lines(glue("stan/solanum-{m}.stan")) |> 
  paste(collapse = "\n")

chkpt_path = glue("checkpoints/chkpt_folder_{m}")

stan_rh_curves = read_rds("data/stan_rh_curves.rds")
init = read_rds("objects/init.rds")

fit_m = inititalize_stan(
  model_code = stan_code,
  data = stan_rh_curves,
  iter_typical = 10, #150,
  path = chkpt_path,
  init = init,
  max_treedepth = 10 #12L
)

tar(paste0(chkpt_path, ".tar"), chkpt_path)
