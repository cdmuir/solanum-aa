library(cmdstanr)
library(glue)
library(readr)
library(stringr)

set_cmdstan_path("cmdstan-2.34.1")

solanum_aa = cmdstan_model("solanum-aa2.stan")

stan_sim = read_rds("stan_sim{n}.rds")
    
ad = c(stan_sim$n, stan_sim$n_comp)

init = list(
  c_a = array(rep(stan_sim$CO2_s, stan_sim$n_comp), dim = ad),
  log_gsw = array(rep(log(stan_sim$g_sw), stan_sim$n_comp), dim = ad),
  sd_A = rep(sd(stan_sim$A), stan_sim$n_comp),
  mean_A = rep(mean(stan_sim$A), stan_sim$n_comp)
)

fit_sim = solanum_aa$sample(
  data = stan_sim,
  chains = 2L,
  parallel_chains = 2L,
  init = list(init, init),
  seed = 20240203,
  iter_warmup = 1e3,
  iter_sampling = 1e3
)
    
fit_sim$save_object("fit_sim{n}.rds")
