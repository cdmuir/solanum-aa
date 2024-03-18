# Compare models using LOOIC
source("r/header.R")

options(mc.cores = 10)

list.files("objects/", pattern = "fit_aa[0-9]+.rds", full.names = TRUE) |>
  map(read_rds) |>
  map(\(.x) {
    ll = .x$draws("log_lik")
    loo(ll, r_eff = relative_eff(exp(ll)))
  }) |>
  write_rds("objects/loo_aa.rds")


# loo_compare(read_rds("objects/loo_aa.rds"))
