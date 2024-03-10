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
stan_rh_curves$id = NULL
stan_rh_curves$light_intensity = NULL
stan_rh_curves$light_treatment = NULL

# Manual changes
stan_rh_curves = c(
  stan_rh_curves,
  
  ## Mapping of lightintensity_x_id to curve
  rh_curves |>
  mutate(
    # character to numeric
    across(c("lightintensity_x_id", "id", "curve", "light_intensity"), \(.x) as.numeric(as.factor(.x))),
    # factor to numeric
    across(c("light_treatment"), as.numeric),
  ) |>
  summarize(
    curve = first(curve),
    id = first(id),
    light_intensity = first(light_intensity),
    light_treatment = first(light_treatment),
    .by = c(lightintensity_x_id, leaf_type)
  ) |>
  pivot_wider(names_from = "leaf_type", values_from = "curve") |>
  arrange(lightintensity_x_id) |>
  dplyr::select(amphi, pseudohypo, id, light_intensity, light_treatment) |>
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
# stan_rh_curves$d1 = 0
# stan_rh_curves$d2 = 0L

# Fit
# m = cmdstan_model("stan/solanum-aa4.stan", dir = "stan/bin")
m = cmdstan_model("stan/solanum-aa5.stan", dir = "stan/bin")

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
    # For solanum-aa4.stan
    # fit = lm(A ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    # For solanum-aa5.stan
    fit = lm(log(A) ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    b = unname(coef(fit))
    tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  }) |>
  as.list()

# started 8:23
# 25% at 8:25, ETA 8 min, so 10x faster
# 50% at 8:27, ETA 8 min, so 10x faster
# 100% at 8:31, ETA 8 min, so 10x faster
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
fit$summary(
  c(
    "rho_resid",
    "b0_aa",
    "b_aa_light_intensity_2000",
    "b_aa_light_treatment_high",
    "b0_log_sigma_aa",
    "log_sigma_resid",
    "log_sigma_aa_id",
    "b_log_sigma_aa_light_intensity_2000",
    "b_log_sigma_aa_light_treatment_high"
  )
) |>
  select(variable, median, q5, q95)

# move these to separate scripts
# calculate aa ----
b_string = "^(b[0-2])\\[([0-9]+)\\]$"
df_b = full_join(
  rh_curves |>
    mutate(curve = as.character(as.numeric(as.factor(curve)))) |>
    summarise(
      id = first(id),
      leaf_type = first(leaf_type),
      light_intensity = first(light_intensity),
      min_scaled_log_gsw = min(scaled_log_gsw),
      max_scaled_log_gsw = max(scaled_log_gsw),
      .by = "curve"
    ),
  
  fit$draws(c("b0", "b1", "b2")) |>
    as_draws_df() |>
    pivot_longer(matches(b_string)) |>
    mutate(
      term = str_replace(name, b_string, "\\1"),
      curve = str_replace(name, b_string, "\\2")
    ) |>
    select(-name) |>
    pivot_wider(names_from = "term"),
  
  by = join_by(curve)
)

df_aa = df_b |>
  mutate(ex_scaled_log_gsw = min_scaled_log_gsw * (leaf_type == "amphi") +
           max_scaled_log_gsw * (leaf_type == "pseudohypo")) |>
  select(.draw, id, light_intensity, leaf_type, b0:ex_scaled_log_gsw) |>
  pivot_wider(names_from = "leaf_type", values_from = b0:ex_scaled_log_gsw) |>
  mutate(
    aa1 = aa_int(ex_scaled_log_gsw_amphi, b0_amphi, b0_pseudohypo, b1_amphi, 
                 b1_pseudohypo, b2_amphi, b2_pseudohypo),
    aa2 = aa_int(ex_scaled_log_gsw_pseudohypo, b0_amphi, b0_pseudohypo, b1_amphi, 
                 b1_pseudohypo, b2_amphi, b2_pseudohypo),
    aa = aa2 - aa1
  ) |>
  select(id, light_intensity, 
         ex_scaled_log_gsw_pseudohypo, b0_amphi, b0_pseudohypo, b1_amphi, 
         b1_pseudohypo, b2_amphi, b2_pseudohypo, aa1, aa2, aa) |>
  group_by(id, light_intensity) |>
  point_interval(aa)

df_aa |>
  arrange(aa) |>
  print(n = 40) 

df_aa |>
  summarize(aa = mean(aa), .by = "light_intensity")

df_aa |>
  arrange(aa) |>
  mutate(r = row_number()) |>
  ggplot(aes(r, aa, ymin = .lower, ymax = .upper, color = light_intensity)) +
  geom_pointinterval()

# residuals ----
df_resid = full_join(
  crossing(
    rh_curves |>
      mutate(.row = row_number(),
             curve = as.character(as.numeric(as.factor(
               curve
             )))),
    .draw = seq_len(max(df_b$.draw))
  ),
  
  fit$draws(c("b0", "b1", "b2")) |>
    as_draws_df() |>
    pivot_longer(matches(b_string)) |>
    mutate(
      term = str_replace(name, b_string, "\\1"),
      curve = str_replace(name, b_string, "\\2")
    ) |>
    select(-name) |>
    pivot_wider(names_from = "term"),
  
  by = join_by(curve, .draw)
) |>
  mutate(log_A_hat = b0 + b1 * scaled_log_gsw + b2 * scaled_log_gsw ^ 2,
         resid = log_A_hat - log(A))

df_resid1 = full_join(
  df_resid |>
    summarize(resid = median(resid), .by = ".row") |>
    arrange(.row),
  
  rh_curves |>
    mutate(.row = row_number(),
           curve = as.character(as.numeric(as.factor(
             curve
           )))),
  by = join_by(.row)
)

ggplot(df_resid1, aes(A, resid)) +
  geom_point() +
  scale_x_log10()

ggplot(df_resid1, aes(scaled_log_gsw, resid)) +
  geom_point() 

ggplot(df_resid1, aes(sample = resid)) +
  facet_wrap(~ curve, scales = "free_y") +
  stat_qq() +
  stat_qq_line()

# plot curve fits ----
df_curve = df_b |>
  crossing(i = seq(0, 1, length.out = 1e1)) |>
  mutate(scaled_log_gsw = min_scaled_log_gsw + 
           i * (max_scaled_log_gsw - min_scaled_log_gsw),
        log_A = b0 + b1 * scaled_log_gsw + b2 * scaled_log_gsw ^ 2) |>
  group_by(id, curve, scaled_log_gsw, leaf_type, light_intensity) |>
  point_interval(log_A) |>
  mutate(A = exp(log_A), lower_ci = exp(.lower), upper_ci = exp(.upper))

ggplot(mapping = aes(scaled_log_gsw, A, color = light_intensity, linetype = leaf_type, fill = light_intensity)) +
  facet_wrap(~ id) +
  geom_lineribbon(
    data = df_curve,
    mapping = aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.5
  ) +
  geom_point(
    data = df_resid1
  ) +
  scale_y_log10()
