# Prepare actual data for Stan
# Right now, just trying a single accession (LA2172)
# Some changes to data should be made earlier on...
source("r/header.R")

rh_curves = read_rds("data/prepared_rh_curves.rds") |>
  dplyr::select(
    A,
    scaled_log_gsw,
    acc,
    curve,
    acc_id,
    leaf_type,
    light_intensity,
    lightintensity_x_acc_id,
    light_treatment
  ) 

stan_rh_curves = rh_curves |>
  compose_data()

  # commented variables are needed for full IRGA model if I go back to that
    # elapsed,
    # flow = Flow,
    # g_bw = gbw,
    # K,
    # P = Pa,
    # RH,
    # s = S,
    # T_air = Tair,
    # T_leaf = Tleaf,
    # CO2_r,
    # CO2_s,
    # H2O_r,
    # H2O_s
  
# drop variables that will be indexed by higher-level groups
stan_rh_curves$acc = NULL
stan_rh_curves$acc_id = NULL
stan_rh_curves$light_intensity = NULL
stan_rh_curves$light_treatment = NULL

# Manual changes
stan_rh_curves = c(
  stan_rh_curves,
  
  ## Mapping of lightintensity_x_acc_id to curve
  rh_curves |>
    mutate(
      # character to numeric
      across(c("lightintensity_x_acc_id", "acc", "acc_id", "curve", "light_intensity"),
             \(.x) as.numeric(as.factor(.x))),
      # factor to numeric
      across(c("light_treatment"), as.numeric),
    ) |>
    summarize(
      curve = first(curve),
      acc = first(acc),
      acc_id = first(acc_id),
      light_intensity = first(light_intensity),
      light_treatment = first(light_treatment),
      .by = c(lightintensity_x_acc_id, leaf_type)
    ) |>
    pivot_wider(names_from = "leaf_type", values_from = "curve") |>
    arrange(lightintensity_x_acc_id) |>
    dplyr::select(amphi, pseudohypo, acc, acc_id, light_intensity, light_treatment) |>
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
    as.list()
  #,
  
  # why do I have this?
  # rh_curves |>
  #   mutate(across(c("curve", "leaf_type"), \(.x) as.numeric(as.factor(.x)))) |>
  #   summarize(
  #     leaf_type = first(leaf_type),
  #     .by = c(curve)
  #   ) |>
  #   arrange(curve) |>
  #   dplyr::select(leaf_type) |>
  #   as.list()
  # 
)

stan_rh_curves$n_pts = rh_curves |>
  summarise(n = n(), .by = "curve") |>
  pull(n)

assert_false(any(duplicated(names(stan_rh_curves))))

write_rds(stan_rh_curves, "data/stan_rh_curves.rds")
