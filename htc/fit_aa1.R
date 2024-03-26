library(chkptstanr)
library(cmdstanr)
library(dplyr)
library(purrr)
library(readr)
library(stringr)
library(tibble)

args = commandArgs(trailingOnly = TRUE)

set_cmdstan_path("cmdstan-2.34.1")

# m = cmdstan_model("solanum-aa1.stan")
stan_code = read_lines("solanum-aa1.stan") |> paste(collapse = "\n")

chkpt_path = paste0("chkpt_folder_aa1_", args[1])
if (!dir.exists(chkpt_path)) {
  create_folder(folder_name = chkpt_path)
}

rh_curves = read_rds("prepared_rh_curves.rds")
stan_rh_curves = read_rds("stan_rh_curves.rds")

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

# Old version without checkpoint
# fit_m = m$sample(
#   data = stan_rh_curves,
#   chains = 4L,
#   parallel_chains = 4L,
#   init = list(init, init, init, init),
#   seed = 898932814,
#   iter_warmup = 4e3,
#   iter_sampling = 4e3,
#   thin = 4e0,
#   max_treedepth = 12L
# )
# 
# fit_m$save_object("fit_aa1.rds")

# New version with checkpoint
fit_m = chkpt_stan(
  model_code = stan_code, 
  data = stan_rh_curves,
  iter_warmup = 4e3,
  iter_sampling = 4e3,
  iter_per_chkpt = 1e2,
  iter_typical = 2e2,
  parallel_chains = 1,
  control = list(max_treedepth = 12L),
  init = init,
  path = chkpt_path
)

draws = combine_chkpt_draws(object = fit_m)

write_rds(draws, paste0("draws_aa1_", args[1], ".rds"))
tar(paste0(chkpt_path, ".tar"), chkpt_path)
