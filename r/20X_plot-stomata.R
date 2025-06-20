source("r/header.R")

fit_aa = read_rds("objects/fit_aa.rds")
accession_gedi = read_rds("data/accession-gedi.rds") |>
  rename(acc = accession)
d1 = fit_aa$data

plant_info = read_rds("data/plant-info.rds") |>
  dplyr::select(accession, replicate, light_treatment)

df_stomata = read_rds("data/stomata.rds") |>
  left_join(plant_info, by = join_by(accession, replicate)) |>
  rename(acc = accession) |>
  filter(acc %in% unique(d1$acc))

df_gmax_ratio = df_stomata |>
  ungroup() |>
  summarize(gmax_ratio = mean(gmax_ratio, na.rm = TRUE),
            .by = c("acc", "light_treatment")) |>
  mutate(
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low")
  )

df_sd = df_stomata |>
  ungroup() |>
  summarize(total_stomatal_density_mm2 = mean(total_stomatal_density_mm2, na.rm = TRUE),
            .by = c("acc", "light_treatment"))

ggplot(df_gmax_ratio, aes(gmax_ratio, color = Growth, fill = Growth)) +
  facet_grid(Growth ~ . ) +
  geom_histogram(binwidth = 0.05, alpha = 0.5) +
  scale_fill_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  xlab("stomatal ratio") +
  xlim(0, 0.6) +
  # xlab(expression(italic(g)[paste(max, ",", ratio)])) +
  theme(
    legend.position = "none"
  )

ggsave("figures/aa-plot3.pdf", width = 6, height = 4, device = cairo_pdf, bg = "transparent")

ggplot(df_sd, aes(light_treatment, total_stomatal_density_mm2, group = acc)) +
  geom_line() + 
  geom_point() 

df_aa_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())


df_aa_pred1 = posterior_epred(fit_aa, newdata = df_aa_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_aa_new, by = "row") |>
  group_by(acc, light_treatment, light_intensity) |>
  point_interval(aa) |>
  left_join(df_gmax_ratio, by = join_by(acc, light_treatment)) |>
  left_join(df_sd, by = join_by(acc, light_treatment)) 

ggplot(df_aa_pred1, aes(gmax_ratio, aa, color = light_treatment)) +
  facet_grid(rows = vars(light_intensity)) +
  geom_point() 

ggplot(df_aa_pred1, aes(total_stomatal_density_mm2, aa, color = light_treatment)) +
  facet_grid(rows = vars(light_intensity)) +
  geom_point() 

df_gmax_ratio |>
  left_join(accession_gedi, by = "acc") |>
  ggplot(aes(pai, gmax_ratio, color = light_treatment)) +
  geom_point()

df_sd |>
  left_join(accession_gedi, by = "acc") |>
  ggplot(aes(pai, total_stomatal_density_mm2, color = light_treatment)) +
  geom_point()
