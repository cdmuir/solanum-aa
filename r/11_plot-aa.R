source("r/header.R")

fit_aa = read_rds("objects/fit_aa1.rds")
d1 = fit_aa$data

# code for brms model
df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000"),
  se_aa = 0
) |>
  mutate(row = row_number())

# Estimates of AA for each accession in each light treatment and light intensity
df_aa_pred1 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_new, by = "row") |>
  group_by(acc, light_treatment, light_intensity) |>
  point_interval(aa) |>
  group_by(light_treatment, light_intensity) |>
  arrange(aa) |>
  mutate(
    x = order(aa),
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low"),
    Measurement = light_intensity |>
      factor(levels = c("150", "2000")) |>
      fct_recode(low = "150", high = "2000")
  )

# Estimate difference in AA between high and low light treatment
# This would need to be changed for model where effect of light_treatment varies
# by accession
df_new = crossing(
  acc = first(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000"),
  se_aa = 0
) |>
  mutate(row = row_number())

df_aa_pred2 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  pivot_wider(names_from = light_treatment, values_from = aa) |>
  mutate(d_aa = high - low) |>
  group_by(light_intensity) |>
  point_interval(d_aa)

# Estimate difference in AA between high and low light intensity
df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000"),
  se_aa = 0
) |>
  mutate(row = row_number())

df_aa_pred3 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  pivot_wider(names_from = light_intensity, values_from = aa) |>
  mutate(d_aa = `2000` - `150`) |>
  group_by(acc, light_treatment) |>
  point_interval(d_aa)

df_aa_text = df_aa_pred1 |>
  ungroup() |>
  summarize(
    x = 0,
    label = sprintf("%.3f", round(mean(aa), 3)),
    aa = 0.30,
    .by = c("Growth", "Measurement")
  )

fig_aa = ggplot(
  df_aa_pred1,
  aes(
    x,
    aa,
    ymin = .lower,
    ymax = .upper,
    group = acc,
    color = Growth,
    shape = Measurement
  )
) +
  facet_grid(Measurement ~ Growth) +
  geom_pointinterval(fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = df_aa_text,
            aes(x, aa, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1) +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  ylab("amphi advantage") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

write_rds(fig_aa, "objects/fig_aa.rds")

ggsave(
  "figures/aa.pdf",
  fig_aa,
  width = 6,
  height = 4,
  device = cairo_pdf,
  bg = "transparent"
)

# Predictions for values in ms
write_rds(df_aa_pred1, "objects/df_aa_pred1.rds")
write_rds(df_aa_pred2, "objects/df_aa_pred2.rds")
write_rds(df_aa_pred3, "objects/df_aa_pred3.rds")

# Plot for amphistomy-proposal
# set.seed(07102025)
# amphi_advantage = ggplot(df_aa_pred1, aes(light_intensity, aa, color = light_treatment, shape = light_intensity)) +
#   facet_wrap(~Growth) +
#   geom_jitter(width = 0.1, height = 0, size = 2, fill = "white") +
#   stat_summary(
#     fun.data = mean_cl_boot,
#     fun.args = list(conf.int = 0.95),
#     geom = "pointinterval",
#     position = position_nudge(x = -0.3),
#     fill = "white", size = 12, linewidth = 1
#   ) +
#   geom_hline(yintercept = 0, linetype = "dashed") +
#   scale_color_manual(values = c("low" = "tomato4", "high" = "tomato")) +
#   scale_shape_manual(values = c("150" = 19, "2000" = 21)) +
#   xlab(expression(Light~intensity~(mu*mol~m^-2~s^-1))) +
#   ylab("Amphistomy advantage") +
#   theme(legend.position = "none", strip.background = element_blank())
# 
# write_rds(amphi_advantage, "objects/amphi_advantage.rds")
