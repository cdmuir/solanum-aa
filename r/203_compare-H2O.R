# Compare simulated to estimated [H2O] values
source("r/header.R")

df_sim = read_rds("synthetic-data/df_sim.rds")
fit_sim = read_rds("objects/fit_sim.rds")

# w_0 and H2O_r
par_string = "^w_0\\[([0-9]+),([0-9]+),([0-9]+)\\]$"

df_w0 = full_join(
  # Simulated w_0
  df_sim |>
    select(leaf_type, rep, pts, w_0, H2Or_sim = H2O_r),
  
  # Estimated w_0
  fit_sim$draws("w_0") |>
    as_draws_df() |>
    pivot_longer(starts_with("w_0"), values_to = "w_0") |>
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
    summarize(w0_est = median(w_0), .by = c(pts, rep, leaf_type)),
  by = join_by(leaf_type, rep, pts)
)

# True w_0 versus simulated H2O_r
ggplot(df_w0, aes(w_0, H2Or_sim, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# Simulated H2O_r versus estimated w_0
ggplot(df_w0, aes(H2Or_sim, w0_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# True w_0 versus estimated w_0
ggplot(df_w0, aes(w_0, w0_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# w_a and H2O_s
par_string = "^w_a\\[([0-9]+),([0-9]+),([0-9]+)\\]$"

df_wa = full_join(
  # Simulated w_a
  df_sim |>
    select(leaf_type, rep, pts, w_a, H2Os_sim = H2O_s),
  
  # Estimated w_a
  fit_sim$draws("w_a") |>
    as_draws_df() |>
    pivot_longer(starts_with("w_a"), values_to = "w_a") |>
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
    summarize(wa_est = median(w_a), .by = c(pts, rep, leaf_type)),
  by = join_by(leaf_type, rep, pts)
)

# True w_a versus simulated H2O_s
ggplot(df_wa, aes(w_a, H2Os_sim, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()

# Simulated H2O_s versus estimated w_a
ggplot(df_wa, aes(H2Os_sim, wa_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point()

# True w_a versus estimated w_a
ggplot(df_wa, aes(w_a, wa_est, color = leaf_type)) +
  facet_wrap(~ rep) +
  geom_abline(slope = 1, intercept = 0) +
  geom_point() +
  coord_equal()
