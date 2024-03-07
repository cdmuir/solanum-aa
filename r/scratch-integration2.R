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
  )

# check
mean(rh_curves$scaled_log_gsw)
sd(rh_curves$scaled_log_gsw)

# WORKING HERE
# ALMOST THERE - JUST NEED TO COMPOSE CURVE BEFORE TRANSFERRING TO LIST
rh_curves |>
  # select variables needed for Stan
  dplyr::select(
    A,
    scaled_log_gsw,
    curve,
    id,
    leaf_type,
    light_intensity,
    lightintensity_x_id
  ) |>
  summarize(
    curve = first(curve),
    .by = c(light_intensity, id)
  )


stan_rh_curves = rh_curves |>
  # select variables needed for Stan
  dplyr::select(
    A,
    scaled_log_gsw,
    curve,
    id,
    leaf_type,
    light_intensity,
    lightintensity_x_id
  ) |>
  compose_data()

# Manual changes
stan_rh_curves$n_pts = rh_curves |>
  summarise(n = n(), .by = "curve") |>
  pull(n)
stan_rh_curves$d1 = 0
stan_rh_curves$d2 = 0L

# Fit
m = cmdstan_model("stan/solanum-aa4.stan", dir = "stan/bin")

fit = m$sample(
  data = stan_rh_curves,
  chains = 2L,
  parallel_chains = 2L,
  iter_warmup = 1e3,
  iter_sampling = 1e3,
  refresh = 2e1
)

fit$profiles()
