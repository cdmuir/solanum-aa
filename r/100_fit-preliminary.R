# Fit preliminary model to each RH curve to identify AA outliers
source("r/header.R")

rh_curves = read_rds("data/prepared_rh_curves.rds") |>
  mutate(log_A = log(A))

# This model worked well for LA2172. Scaling up to all accessions
# started at 1:20 PM
# 100 / 4000 at 1:45 PM, so ETA is 40 * 25 = 1000 minutes
fit_aa0 = brm(
  bf(log_A ~ poly(log_gsw, 2) + (poly(log_gsw, 2) | curve),
     sigma ~ 1 + (1 | curve)),
  family = gaussian(),
  data = rh_curves,
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  silent = 0,
  iter = 4e3,
  thin = 2e0,
  max_treedepth = 12
)

write_rds(fit_aa0, "objects/fit_aa0.rds")

# Check convergence
summary(fit_aa0)

fit_aa0 |>
  spread_draws(r_curve[curve, term]) |>
  summarize(rhat = posterior::rhat(r_curve)) |>
  dplyr::filter(rhat > 1.05)
