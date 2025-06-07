# Plot RH curve fits
# needs work. I copied from scratch script. df_b is defined in 304_calculate-aa.R
# 
source("r/header.R")

df_curve = df_b |>
  crossing(i = seq(0, 1, length.out = 1e1)) |>
  mutate(scaled_log_gsw = min_scaled_log_gsw + 
           i * (max_scaled_log_gsw - min_scaled_log_gsw),
         log_A = b0 + b1 * scaled_log_gsw + b2 * scaled_log_gsw ^ 2) |>
  group_by(id, curve, scaled_log_gsw, leaf_type, light_intensity) |>
  point_interval(log_A) |>
  mutate(A = exp(log_A), lower_ci = exp(.lower), upper_ci = exp(.upper))

ggplot(mapping = aes(scaled_log_gsw, A, color = light_intensity, linetype = leaf_type, fill = light_intensity)) +
  facet_wrap(~ id) +
  geom_lineribbon(
    data = df_curve,
    mapping = aes(ymin = lower_ci, ymax = upper_ci),
    alpha = 0.5
  ) +
  geom_point(
    data = df_resid1
  ) +
  scale_y_log10()
