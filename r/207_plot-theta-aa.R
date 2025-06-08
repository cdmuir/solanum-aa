# not totally working yet, but will plot theta vs AA
source("r/header.R")

# phylogenetic model 
fit_aq = read_rds("objects/fit_aq.rds")

# current nonphylogenetic model
# fit_aq = read_rds("../solanum-aq/objects/fit1.rds")

df_new = crossing(
  acc = fit_aq$data$acc,
  light_treatment = unique(fit_aq$data$light_treatment),
  .Qabs = 0
) |>
  mutate(row = row_number())

formula(fit_aq)

df_aq_pred1 = posterior_epred(fit_aq, newdata = df_new, re_formula = ~ (light_treatment | ID2 | acc), nlpar = "logitThetaJ") |> t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "logitThetaJ") |>
  full_join(df_new, by = "row") |>
  group_by(acc, light_treatment) |>
  point_interval(logitThetaJ) |>
  group_by(light_treatment) 

full_join(
  dplyr::select(df_aq_pred1, acc, light_treatment, logitThetaJ),
  dplyr::select(df_aa_pred1, acc, light_treatment, light_intensity, aa),
  by = join_by(acc, light_treatment)
  ) |>
  ggplot(aes(plogis(logitThetaJ), aa, color = light_treatment)) +
  facet_grid(~ light_intensity) +
  geom_point()



df_aq_pred2 = posterior_epred(fit_aq, newdata = df_new, re_formula = ~ (light_treatment | ID2 | acc), nlpar = "logitPhiJ") |> t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "logitPhiJ") |>
  full_join(df_new, by = "row") |>
  group_by(acc, light_treatment) |>
  point_interval(logitPhiJ) |>
  group_by(light_treatment) 

full_join(
  dplyr::select(df_aq_pred2, acc, light_treatment, logitPhiJ),
  dplyr::select(df_aa_pred1, acc, light_treatment, light_intensity, aa),
  by = join_by(acc, light_treatment)
) |>
  ggplot(aes(logitPhiJ, aa, color = light_treatment)) +
  facet_grid(~ light_intensity) +
  geom_point()
