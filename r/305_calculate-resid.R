# Estimate and inspect residuals
# Still needs work, I just copied code from scratch script
source("r/header.R")

fit_dat = read_rds("objects/fit_dat.rds")

rh_curves = read_rds("data/prepared_rh_curves.rds")
stan_rh_curves = read_rds("data/stan_rh_curves.rds")

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
  
  fit_dat$draws(c("b0", "b1", "b2")) |>
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
