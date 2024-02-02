# Compare actual to estimated hyperparameter values
source("r/header.R")

aa_hyperpars = read_rds("objects/aa_hyperpars.rds")
fit_sim = read_rds("objects/fit_sim.rds")

pdf("figures/fit_sim_areaplot.pdf",
    width = 4,
    height = 4)

c(
  "sigma_c",
  "sigma_w",
  "mu_intercept",
  "mu_slope",
  "sigma_error_intercept",
  "b_autocorr_c",
  "b_autocorr_w"
) |>
  map(\(.hpar) {
    fit_sim$draws(.hpar) |>
      mcmc_areas(prob = 0.95) +
      geom_vline(xintercept = aa_hyperpars[[.hpar]])
  })

dev.off()
