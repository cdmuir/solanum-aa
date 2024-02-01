# Fit Stan model to synthetic data

source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")

stan_sim = list(
  n_pts = length(unique(df_sim$pts)),
  n_rep = length(unique(df_sim$rep)),
  n_leaf_type = length(unique(df_sim$leaf_type))
)

stan_sim$log_gsw = with(stan_sim, array(df_sim$log_gsw, dim = c(n_pts, n_rep, n_leaf_type)))

stan_sim$log_A = with(stan_sim, array(df_sim$log_A, dim = c(n_pts, n_rep, n_leaf_type)))

solanum_aa = cmdstan_model("stan/solanum-aa.stan", dir = "stan/bin")

fit_sim = solanum_aa$sample(
  data = stan_sim,
  chains = 4L,
  parallel_chains = 4L,
)

fit_sim$save_object("objects/fit_sim.rds")
