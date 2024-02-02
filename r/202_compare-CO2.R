# Compare simulated to estimated [CO2] values
source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")
fit_sim = read_rds("objects/fit_sim.rds")

par_string = "^c_a\\[([0-9]+),([0-9]+),([0-9]+)\\]$"

df_ca = full_join(
  # Simulated c_a
  df_sim |>
    select(leaf_type, rep, pts, c_a, CO2s_sim = CO2_s),
  
  # Estimated c_a
  fit_sim$draws("c_a") |>
    as_draws_df() |>
    pivot_longer(starts_with("c_a"), values_to = "c_a") |>
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
    summarize(ca_est = median(c_a), .by = c(pts, rep, leaf_type)),
  by = join_by(leaf_type, rep, pts)
)

# True c_a versus simulated CO2_s
ggplot(df_ca, aes(c_a, CO2s_sim, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# Simulated CO2_s versus estimated c_a
ggplot(df_ca, aes(CO2s_sim, ca_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()

# True c_a versus estimated c_a
ggplot(df_ca, aes(c_a, ca_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()
