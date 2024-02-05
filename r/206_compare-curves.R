# Compare simulated to estimated A-gsw curves
# NOT SURE HOW/IF I SHOULD SCALE THIS UP FOR MANY SIMULATIONS
source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim0001.rds")
fit_sim = read_rds("objects/fit_sim0001.rds")

# Overall slope and intercept
df_mu = c("mu_intercept",
          "mu_intercept_low_light",
          "mu_slope",
          "mu_slope_low_light") |>
  fit_sim$draws() |>
  as_draws_df()

# id level coefficients
df_id = fit_sim$draws(c("b_intercept_id", "b_intercept_low_light_id",
                        "b_slope_id", "b_slope_low_light_id")) |>
  as_draws_df() |>
  pivot_longer(starts_with("b_")) |>
  mutate(id = LETTERS[as.numeric(str_extract(name, "[0-9]+"))],
         name = str_remove(name, "\\[[0-9]+\\]")) |>
  pivot_wider()

# curve level coefficients
df_curve = fit_sim$draws("b_intercept_error") |>
  as_draws_df() |>
  pivot_longer(starts_with("b_intercept_error"), 
               values_to = "b_intercept_error") |>
  mutate(
    id = LETTERS[
      str_replace(name, "^b_intercept_error\\[([0-9]+),([0-9]+)\\]$", "\\1") |>
        as.numeric()
    ],
    lt = str_replace(name, "^b_intercept_error\\[([0-9]+),([0-9]+)\\]$", "\\2"),
    leaf_type = case_when(lt == 1 ~  "amphi",
                          lt == 2 ~ "pseudohypo")
  ) |>
  select(-name, -lt)

df_pred = df_mu |>
  full_join(df_id, by = join_by(.chain, .iteration, .draw)) |>
  full_join(df_curve, by = join_by(id, .chain, .iteration, .draw)) |>
  full_join(select(df_sim, id, pts, leaf_type, light_treatment, g_sw), 
            by = join_by(id, leaf_type),
            relationship = "many-to-many") |>
  mutate(
    A = mu_intercept + b_intercept_id + b_intercept_error + 
      (mu_intercept_low_light + b_intercept_low_light_id) * 
        (light_treatment == "low") +
      (mu_slope + b_slope_id + (mu_slope_low_light + b_slope_low_light_id) * 
         (light_treatment == "low")) * log(g_sw)
  ) |>
  summarize(
    g_sw = first(g_sw),
    .lower = matrix(quantile(A, probs = 0.025), ncol = 1)[,1],
    .upper = matrix(quantile(A, probs = 0.975), ncol = 1)[,1],
    A = median(A),
    .by = c(light_treatment, leaf_type, id, pts)
  )

ggplot(mapping = aes(x = g_sw, color = leaf_type, linetype = light_treatment)) +
  facet_wrap(~ id) +
  geom_ribbon(
    data = df_pred,
    mapping = aes(y = A, min = .lower, ymax = .upper, fill = leaf_type),
    alpha = 0.5
  ) +
  geom_point(
    data = df_sim, mapping = aes(y = A_hat, shape = light_treatment)
  ) +
  scale_x_log10()
