# Calculated estimated AA for each replicate from posterior
# needs some work - I just copied it from a scratch script
source("r/header.R")

fit_dat = read_rds("objects/fit_dat.rds")

rh_curves = read_rds("data/prepared_rh_curves.rds")
stan_rh_curves = read_rds("data/stan_rh_curves.rds")

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
  
  fit_dat$draws(c("b0", "b1", "b2")) |>
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

