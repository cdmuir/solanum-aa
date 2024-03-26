# Fit model to actual data
source("r/header.R")
library(chkptstanr)

# m = cmdstan_model("stan/solanum-aa1.stan", dir = "stan/bin")

rh_curves = read_rds("data/prepared_rh_curves.rds")
stan_rh_curves = read_rds("data/stan_rh_curves.rds")

b = rh_curves |>
  mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
  split(~ curve) |>
  map_dfr(\(.x) {
    fit = lm(log(A) ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    b = unname(coef(fit))
    tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  })

Mu_curve = apply(b, 2, mean)

B_curve = b |>
  mutate(b0 = b0 - Mu_curve[1], b1 = b1 - Mu_curve[2], b2 = b2 - Mu_curve[3])

R_curve = cor(B_curve)

log_sigma_curve = log(diag(cov(B_curve)))

init = list(
  Mu_curve = Mu_curve,
  R_curve = R_curve,
  log_sigma_curve = log_sigma_curve,
  B_curve = B_curve
)

stan_code = read_lines("stan/solanum-aa1.stan") |> paste(collapse = "\n")
path = create_folder(folder_name = "chkpt_folder_aa1")
fit_aa1 = chkpt_stan(
  model_code = stan_code, 
  data = stan_rh_curves,
  iter_warmup = 10,
  iter_sampling = 10,
  thin = 2,
  iter_per_chkpt = 2,
  iter_typical = 5,
  chains = 4,
  parallel_chains = 4,
  chkpt_progress = TRUE,
  init = list(init, init, init, init),
  seed = 898932815,
  path = path
)


draws = combine_chkpt_draws(object = fit_aa1)


fit_aa1 = m$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  init = list(init, init),
  seed = 898932815,
  iter_warmup = 2e3,
  iter_sampling = 2e3,
  refresh = 4e1
)

fit_aa1$save_object("objects/fit_aa1.rds")
