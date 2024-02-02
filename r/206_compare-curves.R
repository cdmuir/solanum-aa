# Compare simulated to estimated A-gsw curves
source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")
fit_sim = read_rds("objects/fit_sim.rds")

# Overall slope and intercept
df_mu = fit_sim$draws(c("mu_slope", "mu_intercept")) |>
  as_draws_df()

# Rep level coefficients
df_rep = fit_sim$draws("intercept") |>
  as_draws_df() |>
  pivot_longer(starts_with("intercept"), values_to = "intercept") |>
  mutate(
    rep = str_c(
      "r",
      str_replace(name, "^intercept\\[([0-9]+),([0-9]+)\\]$", "\\1") |>
        str_pad(2L, "left", "0")
    ),
    lt = str_replace(name, "^intercept\\[([0-9]+),([0-9]+)\\]$", "\\2"),
    leaf_type = case_when(lt == 1 ~  "amphi",
                          lt == 2 ~ "pseudohypo")
  ) |>
  select(-name, -lt)

df_pred = df_mu |>
  full_join(df_rep, by = join_by(.chain, .iteration, .draw)) |>
  full_join(select(df_sim, rep, pts, leaf_type, g_sw), by = join_by(rep, leaf_type),
            relationship = "many-to-many") |>
  mutate(A = exp(mu_intercept + intercept + mu_slope * log(g_sw))) |>
  summarize(
    g_sw = first(g_sw),
    .lower = matrix(quantile(A, probs = 0.025), ncol = 1)[,1],
    .upper = matrix(quantile(A, probs = 0.975), ncol = 1)[,1],
    A = median(A),
    .by = c(leaf_type, rep, pts)
  )

ggplot(mapping = aes(x = g_sw, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_ribbon(
    data = df_pred,
    mapping = aes(y = A, min = .lower, ymax = .upper, fill = leaf_type),
    alpha = 0.5
  ) +
  geom_point(
    data = df_sim, mapping = aes(y = A_hat)
  )
