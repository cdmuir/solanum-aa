# Fit Stan model to synthetic data
source("r/header.R")

solanum_aa = cmdstan_model("stan/solanum-aa2.stan", dir = "stan/bin")

library(furrr)
plan(multisession(workers = 3))

list.files("synthetic-data", pattern = "stan_sim[0-9]{4}.rds", full.names = TRUE) |>
  future_walk(\(.x) {
    n = str_extract(.x, "[0-9]{4}")
    stan_sim = read_rds(.x)
    
    ad = c(stan_sim$n, stan_sim$n_comp)
    
    fit_sim = solanum_aa$sample(
      data = stan_sim,
      chains = 1L,
      parallel_chains = 1L,
      init = list(list(
        c_a = array(rep(stan_sim$CO2_s, stan_sim$n_comp), dim = ad),
        log_gsw = array(rep(log(stan_sim$g_sw), stan_sim$n_comp), dim = ad),
        sd_A = rep(sd(stan_sim$A), stan_sim$n_comp),
        mean_A = rep(mean(stan_sim$A), stan_sim$n_comp)
      )),
      seed = 20240203
    )
    
    fit_sim$save_object(glue("objects/fit_sim{n}.rds"))
    
  })
