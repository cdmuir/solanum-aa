# revise to plot rxn norm in Amax
source("r/header.R")

# phylogenetic models
fit_aa = read_rds("objects/fit_aa.rds")
fit_aq = read_rds("objects/fit_aq.rds")
accession_gedi = read_rds("data/accession-gedi.rds") |>
  rename(acc = accession)
d1 = fit_aa$data

# old nonphylogenetic model
# fit_aq = read_rds("../solanum-aq/objects/fit1.rds")

# Analyze relationship between theta and AA ----

df_aa_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

df_aq_new = crossing(
  acc = fit_aq$data$acc,
  light_treatment = unique(fit_aq$data$light_treatment),
  .Qabs = 0
) |>
  mutate(row = row_number())

df_aa_pred1 = posterior_epred(fit_aa, newdata = df_aa_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  dplyr::select(-matches("^\\.[a-z]+$")) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_aa_new, by = join_by(row)) |>
  dplyr::select(-row)

df_aq_pred1 = posterior_epred(
  fit_aq,
  newdata = df_aq_new,
  re_formula = ~ (light_treatment |
                    ID2 | acc),
  nlpar = "logitThetaJ"
) |> t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  dplyr::select(-matches("^\\.[a-z]+$")) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "logitThetaJ") |>
  full_join(df_aq_new, by = join_by(row)) |>
  dplyr::select(-row, -.Qabs) |>
  mutate(thetaJ = plogis(logitThetaJ))

df_fit1 = full_join(df_aa_pred1, df_aq_pred1, by = join_by(draw, acc, light_treatment)) |>
  filter(if_all(c(aa, logitThetaJ), \(.x) ! is.na(.x))) |>
  split( ~ draw) |>
  map(\(.x) {
    lm(aa ~ light_treatment * light_intensity * thetaJ, data = .x)
    # lm(aa ~ light_treatment * light_intensity * logitThetaJ, data = .x)
  }) |>
  map(coef) |>
  map_dfr(\(.x) {
    tibble(
      b0_high_150 = .x["(Intercept)"],
      b0_high_2000 = .x["(Intercept)"] + .x["light_intensity2000"],
      b0_low_150 = .x["(Intercept)"] + .x["light_treatmentlow"],
      b0_low_2000 = .x["(Intercept)"] + .x["light_treatmentlow"] +
        .x["light_intensity2000"] + .x["light_treatmentlow:light_intensity2000"],
      b1_high_150 = .x["thetaJ"],
      b1_high_2000 = .x["thetaJ"] + .x["light_intensity2000:thetaJ"],
      b1_low_150 = .x["thetaJ"] + .x["light_treatmentlow:thetaJ"],
      b1_low_2000 = .x["thetaJ"] + .x["light_treatmentlow:thetaJ"] +
        .x["light_intensity2000:thetaJ"] +
        .x["light_treatmentlow:light_intensity2000:thetaJ"]
    )
    
  })

df_fit1 |>
  pivot_longer(
    cols = everything(),
    names_to = c("b", "light_treatment", "light_intensity"),
    names_pattern = "(b\\d)_(.+)_(.+)"
  ) |>
  group_by(b, light_treatment, light_intensity) |>
  point_interval(value) 

# Analyze relationship between change in and theta and change in AA ----

df_fit2 = full_join(
  df_aa_pred1 |>
    pivot_wider(names_from = "light_treatment", values_from = "aa") |>
    mutate(d_aa = high - low, .keep = "unused"),
  df_aq_pred1 |>
    dplyr::select(-logitThetaJ) |>
    pivot_wider(names_from = "light_treatment", values_from = "thetaJ") |>
    mutate(d_thetaJ = high - low, .keep = "unused"),
  by = join_by(draw, acc)
) |>
  filter(if_all(starts_with("d_"), \(.x) !is.na(.x))) |>
  split( ~ draw) |>
  map(\(.x) {
    lm(d_aa ~ light_intensity * d_thetaJ, data = .x)
    # lm(d_aa ~ light_intensity * d_logitThetaJ, data = .x)
  }) |>
  map(coef) |>
  map_dfr(\(.x) {
    tibble(
      b0_150 = .x["(Intercept)"],
      b0_2000 = .x["(Intercept)"] + .x["light_intensity2000"],
      b1_150 = .x["d_thetaJ"],
      b1_2000 = .x["d_thetaJ"] + .x["light_intensity2000:d_thetaJ"],
    )
  })


df_fit2 |>
  pivot_longer(
    cols = everything(),
    names_to = c("b", "light_intensity"),
    names_pattern = "(b\\d)_(.+)"
  ) |>
  group_by(b, light_intensity) |>
  point_interval(value) 

# Plot theta vs AA ----
# Some redundancy with 206_plot-rxn-norms. Need to streamline and condense.

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
  left_join(accession_gedi, by = join_by(acc))


df_aq_pred1 = posterior_epred(
  fit_aq,
  newdata = df_aq_new,
  re_formula = ~ (light_treatment |
                    ID2 | acc),
  nlpar = "logitThetaJ"
) |> t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "logitThetaJ") |>
  full_join(df_aq_new, by = "row") |>
  group_by(acc, light_treatment) |>
  point_interval(logitThetaJ)

full_join(
  dplyr::select(df_aq_pred1, acc, light_treatment, logitThetaJ),
  dplyr::select(df_aa_pred1, acc, light_treatment, light_intensity, aa),
  by = join_by(acc, light_treatment)
) |>
  mutate(
    Growth = light_treatment |>
      factor(levels = c("low", "high")) |>
      fct_recode(sun = "high", shade = "low"),
    Measurement = light_intensity |>
      factor(levels = c("150", "2000")) |>
      fct_recode(low = "150", high = "2000")
  ) |>
  ggplot(aes(plogis(logitThetaJ), aa, color = Growth, shape = Measurement)) +
  facet_grid(Measurement ~ .) +
  geom_point() +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  xlab(expression(theta)) +
  ylab("amphi advantage") +
  theme(
    legend.position = "none"
  )
  
ggsave("figures/aa-plot4.pdf", width = 4, height = 4, device = cairo_pdf, bg = "transparent")



df_aq_pred2 = posterior_epred(
  fit_aq,
  newdata = df_new,
  re_formula = ~ (light_treatment |
                    ID2 | acc),
  nlpar = "logitPhiJ"
) |> t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "logitPhiJ") |>
  full_join(df_aq_new, by = "row") |>
  group_by(acc, light_treatment) |>
  point_interval(logitPhiJ) |>
  group_by(light_treatment)

full_join(
  dplyr::select(df_aq_pred2, acc, light_treatment, logitPhiJ),
  dplyr::select(df_aa_pred1, acc, light_treatment, light_intensity, aa),
  by = join_by(acc, light_treatment)
) |>
  ggplot(aes(logitPhiJ, aa, color = light_treatment)) +
  facet_grid( ~ light_intensity) +
  geom_point()
