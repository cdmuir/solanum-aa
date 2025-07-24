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

sr_rxnorm = ggplot(df_stomata_pred, aes(Growth, gmax_ratio, group = acc, color = Growth)) +
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
)

write_rds(sr_rxnorm, "figures/sr-rxnorm.rds")
