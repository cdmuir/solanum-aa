# Compare simulated to estimated A values
source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")
fit_sim = read_rds("objects/fit_sim.rds")

df_A = full_join(
  # Simulated A
  df_sim |>
    select(leaf_type, rep, pts, A, A_sim = A_hat),
  
  # Estimated A
  fit_sim$draws("A") |>
    as_draws_df() |>
    pivot_longer(starts_with("A"), values_to = "A") |>
    mutate(
      pts = str_c(
        "p",
        str_replace(name, "^A\\[([0-9]+),([0-9]+),([0-9]+)\\]$", "\\1") |>
          str_pad(2L, "left", "0")
      ),
      rep = str_c(
        "r",
        str_replace(name, "^A\\[([0-9]+),([0-9]+),([0-9]+)\\]$", "\\2") |>
          str_pad(2L, "left", "0")
      ),
      lt = str_replace(name, "^A\\[([0-9]+),([0-9]+),([0-9]+)\\]$", "\\3"),
      leaf_type = case_when(lt == 1 ~  "amphi",
                            lt == 2 ~ "pseudohypo")
    ) |>
    select(-name, -lt) |>
    summarize(A_est = median(A), .by = c(pts, rep, leaf_type)),
  by = join_by(leaf_type, rep, pts)
)

# True A versus simulated A
ggplot(df_A, aes(A, A_sim, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# Simulated A versus estimated A
ggplot(df_A, aes(A_sim, A_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# True A versus estimated A
ggplot(df_A, aes(A, A_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()
