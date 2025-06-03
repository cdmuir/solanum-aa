source("r/header.R")

stan_rh_curves = read_rds("data/stan_rh_curves.rds")
stan_rh_curves$acc1 = NULL

m5 = cmdstan_model("stan/solanum-aa5.stan")

fit5 = m5$sample(
  data = stan_rh_curves,
  seed = 508504744,
  chains = 4,
  parallel_chains = 4,
  refresh = 1e2,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  thin = 1
)

fit5$summary()
