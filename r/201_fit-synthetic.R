# Fit Stan model to synthetic data

source("r/header.R")

stan_sim = read_rds("objects/stan_sim.rds")

solanum_aa = cmdstan_model("stan/solanum-aa.stan", dir = "stan/bin")

fit_sim = solanum_aa$sample(
  data = stan_sim,
  chains = 1L,
  parallel_chains = 1L,
  init = list(
    list(g_sw = stan_sim$g_sw)
  )
)

fit_sim$save_object("objects/fit_sim.rds")

fit_sim$summary("S_c")
