# Compare simulated to estimated A values
source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")
fit_sim = read_rds("objects/fit_sim.rds")

par_string = "^g_sw\\[([0-9]+),([0-9]+),([0-9]+)\\]$"

df_gsw = full_join(
  # Simulated g_sw
  df_sim |>
    select(leaf_type, rep, pts, g_sw, gsw_sim = gsw_hat),
  
  # Estimated g_sw
  fit_sim$draws("g_sw") |>
    as_draws_df() |>
    pivot_longer(starts_with("g_sw"), values_to = "g_sw") |>
    mutate(
      pts = str_c(
        "p",
        str_replace(name, par_string, "\\1") |>
          str_pad(2L, "left", "0")
      ),
      rep = str_c(
        "r",
        str_replace(name, par_string, "\\2") |>
          str_pad(2L, "left", "0")
      ),
      lt = str_replace(name, par_string, "\\3"),
      leaf_type = case_when(lt == 1 ~  "amphi",
                            lt == 2 ~ "pseudohypo")
    ) |>
    select(-name, -lt) |>
    summarize(gsw_est = median(g_sw), .by = c(pts, rep, leaf_type)),
  by = join_by(leaf_type, rep, pts)
)

# True g_sw versus simulated g_sw
ggplot(df_gsw, aes(g_sw, gsw_sim, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# Simulated g_sw versus estimated g_sw
ggplot(df_gsw, aes(gsw_sim, gsw_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# True g_sw versus estimated g_sw
ggplot(df_gsw, aes(g_sw, gsw_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()
