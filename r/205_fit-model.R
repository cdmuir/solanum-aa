# Fit phylogenetic model in Stan
source("r/header.R")

aa5 = cmdstan_model("stan/aa5.stan", dir = "stan/bin")

stan_data = read_rds("data/stan_data.rds")

focal_pars = c(
  "b0_aa",
  "b_aa_light_intensity_2000",
  "b_aa_light_treatment_high",
<<<<<<< HEAD
  "rhosq_aa_acc",
  "etasq_aa_acc",
  "log_sigma_aa_acc_id",
  "b0_log_sigma_aa",
  "b_log_sigma_aa_light_intensity_2000",
  "b_log_sigma_aa_light_treatment_high",
  "b_aa_2000_high",
  "rhosq_aa_light_intensity_2000_acc",
  "etasq_aa_light_intensity_2000_acc",
  "rhosq_aa_light_treatment_high_acc",
  "etasq_aa_light_treatment_high_acc",
  "rhosq_aa_2000_high_acc",
  "etasq_aa_2000_high_acc",
=======
  "b_aa_2000_high",
  "rhosq_aa_acc",
  "etasq_aa_acc",
  "sigma_aa_acc_id",
  "b0_log_sigma_aa" ,
  "b_log_sigma_aa_light_intensity_2000",
  "b_log_sigma_aa_light_treatment_high",
>>>>>>> e9e9fcdb876ee7ced60bb82dbeca09d02c4b2a4a
  "nu"
)

x = 5
ad = 0.8
n_divergent = Inf
converged = FALSE
while ((n_divergent > 10 | !converged) & (x < 11)) {
  fit_aa5 = aa5$sample(
    data = stan_data,
    seed = 987587829 + x,
    refresh = x * 2000 * 0.1,
    chains = 4L,
    parallel_chains = 4L,
    iter_warmup = x * 1000,
    iter_sampling = x * 1000,
    thin = x,
    adapt_delta = ad
  )
  
  n_divergent = nuts_params(fit_aa5) |>
    subset(Parameter == "divergent__") |>
    pull(Value) |>
    sum()
  s = fit_aa5$summary(focal_pars)
  converged = all(s$rhat < 1.01)
  x = x + 1
  ad = min(ad * 1.1, 0.99)
  
}

fit_aa5$save_object("objects/fit_aa5.rds")
