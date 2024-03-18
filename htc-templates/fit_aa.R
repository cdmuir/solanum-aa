library(cmdstanr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

set_cmdstan_path("cmdstan-2.34.1")

m = cmdstan_model("solanum-{m}.stan")

rh_curves = read_rds("prepared_rh_curves.rds")
stan_rh_curves = read_rds("stan_rh_curves.rds")

init = rh_curves |>
  mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
  split(~ curve) |>
  map_dfr(\(.x) {{
    fit = lm(log(A) ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    b = unname(coef(fit))
    tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  }}) |>
  as.list()

fit_m = m$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  init = list(init, init),
  seed = 898932814,
  iter_warmup = 2e3,
  iter_sampling = 2e3,
  thin = 2e0
)
    
fit_m$save_object("fit_{m}.rds")
