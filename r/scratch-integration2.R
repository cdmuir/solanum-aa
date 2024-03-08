# TESTING METHOD OF ESTIMATING AA BY INTEGRATING OVER A-GSW CURVES
# Prepare actual data for Stan
# Right now, just trying a single accession (LA2172)
# Some changes to data should be made earlier on...
source("r/header.R")

rh_curves = read_rds("data/thinned_rh_curves.rds") |>
  mutate(
    RH = H2O_s / li6800_svp(Tair, Pa),
    leaf_type = case_when(
      curve_type == "1-sided RH" ~ "pseudohypo",
      curve_type == "2-sided RH" ~ "amphi"
    )
  ) |>
  separate_wider_delim(acc_id, "-", names = c("acc", "id")) |>
  # make unique ID for each light_intensity within id (change to acc_id later)
  unite("lightintensity_x_id", id, light_intensity, remove = FALSE) |>
  # make unique ID for each curve
  unite("curve", id, leaf_type, light_intensity, remove = FALSE) |>
  filter(acc == "LA2172") |>
  mutate(
    log_gsw = log(gsw),
    scaled_log_gsw = (log_gsw - mean(log_gsw)) / sd(log_gsw)
  ) |>
  dplyr::select(
    A,
    scaled_log_gsw,
    curve,
    id,
    leaf_type,
    light_intensity,
    lightintensity_x_id,
    light_treatment
  ) 

# check
mean(rh_curves$scaled_log_gsw)
sd(rh_curves$scaled_log_gsw)


# Compose data for Stan
stan_rh_curves = rh_curves |>
  compose_data()

# drop variables that will be indexed by higher-level groups
stan_rh_curves$light_intensity = NULL
stan_rh_curves$light_treatment = NULL

# Manual changes
stan_rh_curves = c(
  stan_rh_curves,
  
  ## Mapping of lightintensity_x_id to curve
  rh_curves |>
  mutate(
    # character to numeric
    across(c("lightintensity_x_id", "curve", "light_intensity"), \(.x) as.numeric(as.factor(.x))),
    # factor to numeric
    across(c("light_treatment"), as.numeric),
  ) |>
  summarize(
    curve = first(curve),
    light_intensity = first(light_intensity),
    light_treatment = first(light_treatment),
    .by = c(lightintensity_x_id, leaf_type)
  ) |>
  pivot_wider(names_from = "leaf_type", values_from = "curve") |>
  arrange(lightintensity_x_id) |>
  dplyr::select(amphi, pseudohypo, light_intensity, light_treatment) |>
  as.list(),

  ## min and max scaled_log_gsw by curve
  rh_curves |>
  mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
    summarize(
      min_scaled_log_gsw = min(scaled_log_gsw),
      max_scaled_log_gsw = max(scaled_log_gsw),
      .by = c(curve)
    ) |>
    arrange(curve) |>
    dplyr::select(min_scaled_log_gsw, max_scaled_log_gsw) |>
    as.list(),
  
  rh_curves |>
    mutate(across(c("curve", "leaf_type"), \(.x) as.numeric(as.factor(.x)))) |>
    summarize(
      leaf_type = first(leaf_type),
      .by = c(curve)
    ) |>
    arrange(curve) |>
    dplyr::select(tmp = leaf_type) |>
    as.list()
  
)

stan_rh_curves$n_pts = rh_curves |>
  summarise(n = n(), .by = "curve") |>
  pull(n)
stan_rh_curves$d1 = 0
stan_rh_curves$d2 = 0L

  
# Fit
m = cmdstan_model("stan/solanum-aa4.stan", dir = "stan/bin")

# init = list(
#   b0 = stan_rh_curves$tmp,
#   b1 = rep(0, stan_rh_curves$n_curve),
#   b2 = rep(0, stan_rh_curves$n_curve),
#   mu_aa = 0,
#   sigma_aa = 0.1,
#   aa = rep(0, stan_rh_curves$n_lightintensity_x_id)
# )

init = rh_curves |>
  mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
  split(~ curve) |>
  map_dfr(\(.x) {
    fit = lm(A ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    b = unname(coef(fit))
    tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  }) |>
  as.list()

fit = m$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  refresh = 2e1,
  init = list(init, init)
)

fit$save_object("objects/test-integration.rds")

fit = read_rds("objects/test-integration.rds")
fit$profiles()
