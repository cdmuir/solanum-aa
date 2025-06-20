source("r/header.R")

fit_aa = read_rds("objects/fit_aa.rds")
accession_gedi = read_rds("data/accession-gedi.rds") |>
  rename(acc = accession)
d1 = fit_aa$data

# Some redundancy with 206_plot-rxn-norms. Need to streamline and condense.
df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

df_aa_pred1 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
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

ggplot(df_aa_pred1, aes(pai, aa, color = Growth, shape = Measurement)) +
  facet_grid(Measurement ~ Growth) +
  geom_point() +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  scale_x_log10() +
  scale_x_reverse() +
  xlab("plant area index") +
  ylab("amphi advantage") +
  theme(
    legend.position = "none"
  )

ggsave("figures/aa-plot2.pdf", width = 6, height = 4, device = cairo_pdf, bg = "transparent")

df_aa_pred2 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  pivot_wider(names_from = light_treatment, values_from = aa) |>
  mutate(d_aa = high - low) |>
  group_by(acc, light_intensity) |>
  point_interval(d_aa) |>
  left_join(accession_gedi, by = join_by(acc))

ggplot(df_aa_pred2, aes(pai, d_aa)) +
  facet_grid(light_intensity ~ .) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") 

df_pred2 |>
  ungroup() |>
  dplyr::select(acc, light_intensity, d_aa) |>
  pivot_wider(
    names_from = light_intensity,
    values_from = d_aa
  ) |>
  ggplot(aes(`150`, `2000`, label = acc)) +
  geom_point() +
  geom_label()
