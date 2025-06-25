# Plot effect LMA on AA
source("r/header.R")

fit_aa = read_rds("objects/fit_aa2.rds")

# Estimated llma for each acc, light_intensity, and light_treatment
df_new1 = crossing(
    acc = unique(fit_aa$data$acc),
    light_intensity = unique(fit_aa$data$light_intensity),
    light_treatment = unique(fit_aa$data$light_treatment)
  ) |>
  mutate(row = row_number())

df_pred1 = posterior_epred(fit_aa, newdata = df_new1, resp = "llma") |>
  as_draws_df() |>
  pivot_longer(starts_with("..."), values_to = "llma") |>
  mutate(row = as.integer(str_remove(name, "...")), .keep = "unused") |>
  full_join(df_new1, by = join_by(row)) |>
  group_by(acc, light_intensity, light_treatment) |>
  point_interval(llma) |>
  mutate(row = row_number())

# Estimated AA for each acc, light_intensity, and light_treatment at estimated llma
df_pred2 = posterior_epred(fit_aa, newdata = df_pred1, resp = "aa") |>
  as_draws_df() |>
  pivot_longer(starts_with("..."), values_to = "aa") |>
  mutate(row = as.integer(str_remove(name, "...")), .keep = "unused") |>
  full_join(df_pred1, by = join_by(row)) |>
  group_by(acc, light_intensity, light_treatment) |>
  point_interval(aa)

# Join
df1 = full_join(
  dplyr::select(
    df_pred2,
    acc,
    light_intensity,
    light_treatment,
    aa,
    aa_lower = .lower,
    aa_upper = .upper
  ),
  dplyr::select(
    df_pred1,
    acc,
    light_intensity,
    light_treatment,
    llma,
    llma_lower = .lower,
    llma_upper = .upper
  ),
  by = join_by(acc, light_intensity, light_treatment)
) |>
  ungroup() |>
  mutate(across(contains("llma"), exp)) |>
  rename_with(~ str_replace(., "^llma", "lma"), .cols = starts_with("llma")) |>
  mutate(
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low"),
    Measurement = light_intensity |>
      factor(levels = c("150", "2000")) |>
      fct_recode(low = "150", high = "2000")
  )

# Regression lines
fit_aa$fit$draws

as_draws_df(fit_aa) %>%
  colnames()
  dplyr::select("b_aa_Intercept", "b_aa_light_intensity2000", "bsp_aa_millma",
                "bsp_aa_millma:light_intensity2000")
  
# Plot
ggplot(df1, aes(x = lma, y = aa, color = Growth, shape = Measurement)) +
  # geom_interval(aes(ymin = aa_lower, ymax = aa_upper), linewidth = 0.5) +
  # geom_interval(aes(xmin = lma_lower, xmax = lma_upper), linewidth = 0.5) +
  geom_point(fill = "white") +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  labs(
    x = expression(paste("leaf mass per area [g ", m^-2, "]")),,
    y = "amphi advantage"
  ) 
