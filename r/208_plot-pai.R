# Plot PAI against AA and LMA
source("r/header.R")

fit_aa1 = read_rds("objects/fit_aa1.rds")
fit_aa2 = read_rds("objects/fit_aa2.rds")
accession_gedi = read_rds("data/accession-gedi.rds") |>
  rename(acc = accession)
d1 = fit_aa1$data

# Some redundancy with 206_plot-rxn-norms. Need to streamline and condense.
df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

df_aa_pred1 = posterior_epred(fit_aa1, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_new, by = "row") |>
  group_by(acc, light_treatment, light_intensity) |>
  point_interval(aa) |>
  left_join(accession_gedi, by = join_by(acc)) |>
  mutate(
    x = order(aa),
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low"),
    Measurement = light_intensity |>
      factor(levels = c("150", "2000")) |>
      fct_recode(low = "150", high = "2000")
  )

fig_aa_pai = ggplot(df_aa_pred1, aes(pai, aa, color = Growth, shape = Measurement)) +
  facet_grid(Measurement ~ Growth) +
  geom_point() +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  scale_x_continuous(breaks = c(0.01, 0.1, 1), trans = reverselog10_trans()) +
  xlab(expression(paste("native plant area index [", m^2~m^-2, "]"))) +
  ylab("amphi advantage") +
  ylim(0, 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  theme(legend.position = "none")

write_rds(fig_aa_pai, "objects/fig_aa_pai.rds")

ggsave(
  "figures/pai-aa.pdf",
  fig_aa_pai,
  width = 6,
  height = 4,
  device = cairo_pdf,
  bg = "transparent"
)

df_new = crossing(acc = unique(d1$acc),
                  light_treatment = c("low", "high")) |>
  mutate(row = row_number())

df_aa_pred2 = posterior_epred(fit_aa2, newdata = df_new, resp = "llma") |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "llma") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  group_by(acc, light_treatment) |>
  point_interval(llma) |>
  left_join(accession_gedi, by = join_by(acc))  |>
  mutate(
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low")
  )

ggplot(df_aa_pred2, aes(pai, exp(llma), color = Growth)) +
  geom_point() +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_x_continuous(breaks = c(0.01, 0.1, 1), trans = reverselog10_trans()) +
  xlab(expression(paste("plant area index [", m^2~m^-2, "]"))) +
  ylab(expression(paste("leaf mass per area [g ", m^-2, "]")))
  
ggsave(
  "figures/pai-lma.pdf",
  width = 6,
  height = 4,
  device = cairo_pdf,
  bg = "transparent"
)
