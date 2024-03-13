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

# Started 6:02
# 14% at 70 minutes, ETA 8h 20m
# 28% at 133 minutes, ETA 7h 55m
# 36% at 173 minutes, ETA 8h 0m
# 41% at 200 minutes, ETA 8h 8m
# 56% at 274 minutes, ETA 8h 9m
fit_aa1 = m$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  init = list(init, init),
  seed = 898932814,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  refresh = 2e1
)

fit_aa1$save_object("objects/fit_aa1.rds")
