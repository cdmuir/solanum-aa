# Fit model to actual data
source("r/header.R")

m = cmdstan_model("stan/solanum-aa1.stan", dir = "stan/bin")

rh_curves = read_rds("data/prepared_rh_curves.rds")
stan_rh_curves = read_rds("data/stan_rh_curves.rds")

init = rh_curves |>
  mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
  split(~ curve) |>
  map_dfr(\(.x) {
    fit = lm(log(A) ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    b = unname(coef(fit))
    tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  }) |>
  as.list()

# Started 8:59
# 10% at 35 minutes, ETA 5h 50m
# 35% at 123 minutes, ETA 5h 51m
# 50% at 175 minutes, ETA 5h 50m
# 72% at 284 minutes, ETA 6h 34m
fit_aa1 = m$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  init = list(init, init),
  seed = 898932814,
  iter_warmup = 2e3,
  iter_sampling = 2e3,
  refresh = 4e1
)

fit_aa1$save_object("objects/fit_aa1.rds")
