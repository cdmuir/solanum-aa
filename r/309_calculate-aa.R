# Calculated estimated AA for each replicate from posterior
source("r/header.R")

# MODEL SELECTION NEEDS TO BE DONE BEFORE THIS STEP
fit_aa = read_rds("objects/fit_aa1.rds")

rh_curves = read_rds("data/prepared_rh_curves.rds")
stan_rh_curves = read_rds("data/stan_rh_curves.rds")

# check that rh_curves and stan_rh_curves are consistent
assert_true(all(rh_curves$scaled_log_gsw == stan_rh_curves$scaled_log_gsw))

# calculate aa ----
b_string = "^(b[0-2])\\[([0-9]+)\\]$"
df_b = full_join(
  rh_curves |>
    mutate(curve = as.character(as.numeric(as.factor(
      curve
    )))) |>
    summarise(
      acc_id = first(acc_id),
      leaf_type = first(leaf_type),
      light_intensity = first(light_intensity),
      min_scaled_log_gsw = min(scaled_log_gsw),
      max_scaled_log_gsw = max(scaled_log_gsw),
      .by = "curve"
    ),
  
  fit_aa$draws(c("b0", "b1", "b2")) |>
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
  mutate(
    ex_scaled_log_gsw = min_scaled_log_gsw * (leaf_type == "amphi") +
      max_scaled_log_gsw * (leaf_type == "pseudohypo")
  ) |>
  dplyr::select(.draw,
                acc_id,
                light_intensity,
                leaf_type,
                b0:ex_scaled_log_gsw) |>
  pivot_wider(names_from = "leaf_type",
              values_from = b0:ex_scaled_log_gsw) |>
  mutate(
    aa1 = aa_int(
      ex_scaled_log_gsw_amphi,
      b0_amphi,
      b0_pseudohypo,
      b1_amphi,
      b1_pseudohypo,
      b2_amphi,
      b2_pseudohypo
    ),
    aa2 = aa_int(
      ex_scaled_log_gsw_pseudohypo,
      b0_amphi,
      b0_pseudohypo,
      b1_amphi,
      b1_pseudohypo,
      b2_amphi,
      b2_pseudohypo
    ),
    aa = aa2 - aa1
  ) |>
  dplyr::select(
    acc_id,
    light_intensity,
    ex_scaled_log_gsw_pseudohypo,
    b0_amphi,
    b0_pseudohypo,
    b1_amphi,
    b1_pseudohypo,
    b2_amphi,
    b2_pseudohypo,
    aa1,
    aa2,
    aa
  ) |>
  group_by(acc_id, light_intensity) |>
  point_interval(aa)

write_rds(df_b, "objects/df_b.rds")
write_rds(df_aa, "objects/df_aa.rds")
