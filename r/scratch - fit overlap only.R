# I'm not sure if this is necessary...I think the idea was just to fit curves where they overlap, but I don't think this is necessary.
source("r/header.R")

aa_summary = read_rds("objects/aa_summary.rds")

# Summarize RH curves to get overlap of log_gsw values ----
trimmed_rh_curves = read_rds("data/trimmed_rh_curves.rds") 

df1 = trimmed_rh_curves |>
  summarize(
    min_log_gsw = min(log_gsw),
    max_log_gsw = max(log_gsw),
    .by = c(
      "acc",
      "acc_id",
      "light_treatment",
      "light_intensity",
      "leaf_type"
    )
  ) |>
  reframe(
    p = seq_len(20),
    log_gsw = seq(
      min_log_gsw,
      max_log_gsw,
      length.out = 20
    ),
    .by = c(
      "acc",
      "acc_id",
      "light_treatment",
      "light_intensity",
      "leaf_type"
    )
  )

# Draws from the curve fits ----
curve_fits_draws = read_rds("objects/curve-fits-draws.rds") |>
  rename(b0 = b_Intercept, b1 = b_polylog_gsw2rawEQTRUE1, b2 = b_polylog_gsw2rawEQTRUE2) |>
  dplyr::select(starts_with("."), matches("b[0-2]"), file) |>
  mutate(accid_leaftype_lightintensity = str_remove(file, ".rds"),
         .keep = "unused") |>
  separate_wider_delim(
    accid_leaftype_lightintensity,
    delim = "_",
    names = c("acc_id", "leaf_type", "light_intensity")
  )

df2 = full_join(df1, curve_fits_draws, by = join_by(acc_id, light_intensity, leaf_type)) |>
  mutate(log_A = b0 + b1 * log_gsw + b2 * log_gsw^2) |>
  select(-b0, -b1, -b2) |>
  group_by(acc, acc_id, light_treatment, light_intensity, leaf_type, p) |>
  point_interval(log_A) |>
  left_join(df1,
            by = join_by(acc, acc_id, light_treatment, light_intensity, leaf_type, p))


.x =  "LA0436-P"
aa_summary |>
  arrange(desc(sd))

ggplot(filter(df2, acc_id == .x), aes(x = log_gsw, y = log_A, color = leaf_type)) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.2) +
  geom_line() +
  geom_point(data = filter(trimmed_rh_curves, acc_id == .x)) +
  labs(
    title = paste("RH Curve Fits for", .x),
    subtitle = glue("AA = {.y}", .y = prettyNum(filter(aa_summary, acc_id == .x)$median, digits = 2)),
    x = "log_gsw",
    y = "log_A"
  ) +
  theme_minimal()
  
