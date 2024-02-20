# Fit Stan model to synthetic data
source("r/header.R")

solanum_aa = cmdstan_model("stan/solanum-aa.stan", dir = "stan/bin")

library(furrr)
plan(multisession(workers = 3))
list.files("synthetic-data", pattern = "stan_sim[0-9]{4}.rds", full.names = TRUE) |>
  future_walk(\(.x) {
    n = str_extract(.x, "[0-9]{4}")
    stan_sim = read_rds(.x)
    
    # This could be useful if I want to pre-scale variables rather than calculate in Stan
    # stan_sim$sd_A = sd(stan_sim$A)
    # stan_sim$mean_A = mean(stan_sim$A)
    # stan_sim$sd_log_gsw = sd(log(stan_sim$g_sw))
    # stan_sim$mean_log_gsw = mean(log(stan_sim$g_sw))
    stan_sim$n_comp = 2 # add this to previous script
    fit_sim = solanum_aa$sample(
      data = stan_sim,
      chains = 1L,
      parallel_chains = 1L,
      # init = list(
      #   list(
      #     c_a = rep(stan_sim$CO2_s, stan_sim$n_comp),
      #     log_gsw = rep(log(stan_sim$g_sw), stan_sim$n_comp),
      #     sd_A = rep(sd(stan_sim$A), stan_sim$n_comp), 
      #     mean_A = rep(mean(stan_sim$A), stan_sim$n_comp)
      #   )
      # ),
      seed = 20240203
    )
    
    fit_sim$save_object(glue("objects/fit_sim{n}.rds"))
    
  })
