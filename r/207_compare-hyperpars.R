# Compare actual to estimated hyperparameter values
source("r/header.R")

aa_hyperpars = read_rds("objects/aa_hyperpars.rds")
n_sim = lapply(aa_hyperpars, length)

assert_true(length(unique(n_sim)) == 1L)
n_sim = n_sim[1]

fit = list.files("objects", "fit_sim[0-9]{4}.rds", full.names = TRUE)
n_fit = str_extract(fit, "[0-9]{4}")

assert_true(length(n_fit) == n_sim)

pars = c(
  "sigma_c",
  "sigma_w",
  "mu_intercept",
  "mu_intercept_low_light",
  "sigma_intercept_id",
  "sigma_intercept_low_light_id",
  "sigma_intercept_error",
  "mu_slope",
  "mu_slope_low_light",
  "sigma_slope_id",
  "sigma_slope_low_light_id",
  "b_autocorr_c",
  "b_autocorr_w"
)

full_join(
  # Estimated values
  n_fit |>
    map_dfr(\(.x) {
      fit_sim = read_rds(glue("objects/fit_sim{.x}.rds"))
      fit_sim$summary(pars) |>
        mutate(sim = glue("sim{.x}"))
    }),
  
  # Simulated values
  aa_hyperpars |>
    as_tibble() |>
    select(all_of(pars)) |>
    mutate(sim = glue("sim{n_fit}")) |>
    pivot_longer(-sim, names_to = "variable", values_to = "simulated"),
  
  by = join_by(sim, variable)
  
) |>
  write_rds("objects/fit_sim_summary_hyperpars.rds")

read_rds("objects/fit_sim_summary_hyperpars.rds") |>
  select(variable, simulated, q5, q95) |>
  filter(simulated < q5)

# Area plots (NOT SURE IF IʻLL USE FOR SIMULATED DATA, BUT MIGHT BE USEFUL FOR LATER)
# pdf("figures/fit_sim_areaplot.pdf",
#     width = 4,
#     height = 4)
# 
# pars |>
#   map(\(.hpar) {
#     fit_sim$draws(.hpar) |>
#       mcmc_areas(prob = 0.95) +
#       geom_vline(xintercept = aa_hyperpars[[.hpar]])
#   })
# 
# dev.off()
