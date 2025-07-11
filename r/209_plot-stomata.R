source("r/header.R")

fit_aa = read_rds("objects/fit_aa1.rds")
fit_stomata = read_rds("objects/fit_stomata.rds")
accession_gedi = read_rds("data/accession-gedi.rds") |>
  rename(acc = accession)
accession_climate = read_rds("data/accession-climate.rds") |>
  rename(acc = accession)
d1 = fit_aa$data

plant_info = read_rds("data/plant-info.rds") |>
  dplyr::select(accession, replicate, light_treatment)

# Posterior of stomatal ratio for each accession in each light treatment
df_stomata_new = crossing(acc = unique(d1$acc),
                          light_treatment = c("low", "high")) |>
  mutate(row = row_number())

df_stomata_pred = c("lower_sd", "upper_sd", "lower_gcl", "upper_gcl") |>
  map(function(resp) {
    posterior_epred(fit_stomata, newdata = df_stomata_new, resp = str_remove(resp, "_")) |>
      t() |>
      as_draws_df() |>
      mutate(row = row_number()) |>
      pivot_longer(cols = -row,
                   names_to = "draw",
                   values_to = resp) |>
      full_join(df_stomata_new, by = "row")
  }) |>
  reduce(full_join, by = join_by(row, draw, acc, light_treatment)) |>
  mutate(
    gmax_ratio = exp(upper_gcl + upper_sd) / (exp(lower_gcl + lower_sd) + exp(upper_gcl + upper_sd)),
    sr = exp(upper_sd) / (exp(lower_sd) + exp(upper_sd))
  ) |>
  group_by(acc, light_treatment) |>
  point_interval(gmax_ratio) |>
  mutate(Growth = light_treatment |>
           factor(levels = c("low", "high")) |>
           fct_recode(sun = "high", shade = "low"))

ggplot(df_stomata_pred, aes(Growth, gmax_ratio, group = acc, color = Growth)) +
  geom_line(color = "black") +
  geom_point() +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  xlab("Growth light intensity") +
  ylab(expression(italic(g)[max]~ratio)) +
  ylim(0, 0.5) +
  theme(
    legend.position = "none", axis.title.x = element_blank() )
  
ggsave("~/Documents/grants/01_preparing/amiphistomy-proposal/figures/light_plasticity1.pdf", width = 3.25/1.5, height = 4/1.5)

ggsave(
  "figures/sr-rxnorm.pdf",
  width = 4,
  height = 4,
  device = cairo_pdf,
  bg = "transparent"
)


# STUFF OF PLOTTING STOMATAL TRAITS FOR AA AND CLIMATE. NOT SURE WHAT I AM GOING TO USE
# Posterior of AA for each accession in each light treatment and light intensity

full_join(
  df_stomata_pred |>
    dplyr::select(acc, light_treatment, sr),
  read_rds("objects/df_aa_pred1.rds") |>
    ungroup() |>
    dplyr::select(acc, aa, light_treatment, light_intensity),
  by = join_by(acc, light_treatment)
) |>
  ggplot(aes(sr, aa, color = light_treatment)) +
  facet_grid(light_intensity ~ .) +
  geom_point() +
  xlab("stomatal ratio") +
  ylab("amphi advantage") 



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
  left_join(accession_climate, by = "acc") |>
  ggplot(aes(ppfd_mol_m2, gmax_ratio, color = light_treatment)) +
  scale_x_log10() +
  geom_point()

df_sd |>
  left_join(accession_gedi, by = "acc") |>
  left_join(accession_climate, by = "acc") |>
  ggplot(aes(ppfd_mol_m2, total_stomatal_density_mm2, color = light_treatment)) +
  # scale_x_log10() +
  geom_point()
