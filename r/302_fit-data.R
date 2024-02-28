# Fit model to actual data
source("r/header.R")

solanum_aa = cmdstan_model("stan/solanum-aa3.stan", dir = "stan/bin")

stan_rh_curves = read_rds("data/stan_rh_curves.rds")

ad = c(stan_rh_curves$n, stan_rh_curves$n_comp)

# started at 9:52 PM
# 1% at 9:55 PM, 11.1 hours
# 2% at 10:04 PM, 10 hours
# 3% at 10:16 PM, 13.3 hours
fit_dat = solanum_aa$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  init = list(list(
    c_a = array(rep(stan_rh_curves$CO2_s, stan_rh_curves$n_comp), dim = ad),
    log_gsw = array(rep(log(stan_rh_curves$g_sw), stan_rh_curves$n_comp), dim = ad),
    sd_A = rep(sd(stan_rh_curves$A), stan_rh_curves$n_comp),
    mean_A = rep(mean(stan_rh_curves$A), stan_rh_curves$n_comp)
  )),
  seed = 898932814,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  refresh = 2e1
  
)

fit_dat$save_object("objects/fit_dat.rds")
    # some checks I can remove later
    # fit_sim$profiles()
    # fit_sim$summary("w")$median -> tmp
    # hist(tmp[80 + 1:80])
