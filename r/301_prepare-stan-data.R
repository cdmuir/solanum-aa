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
  # make unique ID for each leaf_type within id (change to acc_id later)
  unite("leaftype_x_id", id, leaf_type, remove = FALSE) |>
  # make unique ID for each curve
  unite("curve", id, leaf_type, light_intensity, remove = FALSE) |>
  filter(acc == "LA2172")

stan_rh_curves = rh_curves |>
  # select variables needed for Stan
  dplyr::select(
    A,
    curve,
    id,
    leaf_type,
    leaftype_x_id,
    light_intensity,
    light_treatment,
    # elapsed,
    flow = Flow,
    g_bw = gbw,
    g_sw = gsw,
    K,
    P = Pa,
    RH,
    s = S,
    T_air = Tair,
    T_leaf = Tleaf,
    CO2_r,
    CO2_s,
    H2O_r,
    H2O_s
  ) |>
  compose_data()

# Manual changes
stan_rh_curves$n_comp = 2
stan_rh_curves$n_pts = rh_curves |>
  summarise(n = n(), .by = "curve") |>
  pull(n)

write_rds(stan_rh_curves, "data/stan_rh_curves.rds")
