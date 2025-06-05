# Fit model on posterior estimate of AA ----
source("r/header.R")

d1 = read_rds("objects/stan_data_df.rds")
fit1 = brm(bf(
  aa ~ amphi_first + light_treatment * light_intensity + (1 + light_treatment | acc),
  sigma ~ light_intensity
  ), data = d1, backend = "cmdstanr", chains = 1)
summary(fit1)

aa1 = cmdstan_model("stan/aa1.stan", dir = "stan/bin")
aa5 = cmdstan_model("stan/aa5.stan", dir = "stan/bin")

stan_data = read_rds("data/stan_data.rds")

focal_pars = c(
  "b0_aa",
  "b_aa_amphi_first",
  "b_aa_light_intensity_2000",
  "b_aa_light_treatment_high",
  "b_aa_2000_high",
  "rhosq_aa_acc",
  "etasq_aa_acc",
  "sigma_aa_acc_id",
  "b0_log_sigma_aa" ,
  "b_log_sigma_aa_light_intensity_2000",
  "b_log_sigma_aa_light_treatment_high"
)

x = 5
ad = 0.8
n_divergent = Inf
converged = FALSE
while ((n_divergent > 10 | !converged) & (x < 11)) {
  fit_aa1 = aa1$sample(
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
  
  n_divergent = nuts_params(fit_aa1) |>
    subset(Parameter == "divergent__") |>
    pull(Value) |>
    sum()
  s = fit_aa1$summary(focal_pars)
  converged = all(s$rhat < 1.01)
  x = x + 1
  # x = 2 * x
  ad = min(ad * 1.1, 0.99)
  
  
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
  # x = 2 * x
  ad = min(ad * 1.1, 0.99)
}

# if (n_divergent <= 10 & converged) {
#   fit_aa1$save_object(paste0("objects/fit-aa1/", draw_id, ".rds"))
# }

log_lik1 <- fit_aa1$draws("log_lik", format = "draws_matrix")
log_lik5 <- fit_aa5$draws("log_lik", format = "draws_matrix")

loo1 <- loo(log_lik1)
loo5 <- loo(log_lik5)

loo_compare(loo1, loo5)

tibble()
fit1 = brm(aa ~ 1, data = stan_data)
