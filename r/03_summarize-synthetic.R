source("r/header.R")

fit_sim = read_rds("objects/fit_sim.rds")

fit_sim$summary(c("mu_intercept", "slope", "sigma_error_intercept", "sigma_error_resid"))
