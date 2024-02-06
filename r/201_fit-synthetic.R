# Fit Stan model to synthetic data
source("r/header.R")

solanum_aa = cmdstan_model("stan/solanum-aa.stan", dir = "stan/bin")

library(furrr)
plan(multisession(workers = 3))
list.files("synthetic-data", pattern = "stan_sim[0-9]{4}.rds", full.names = TRUE) |>
  future_walk(\(.x) {
    n = str_extract(.x, "[0-9]{4}")
    stan_sim = read_rds(.x)
    
    fit_sim = solanum_aa$sample(
      data = stan_sim,
      chains = 1L,
      parallel_chains = 1L,
      init = list(
        list(c_a = stan_sim$CO2_s, log_gsw = log(stan_sim$g_sw))
      ),
      seed = 20240203
    )
    
    fit_sim$save_object(glue("objects/fit_sim{n}.rds"))
    
  })
