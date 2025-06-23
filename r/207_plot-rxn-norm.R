source("r/header.R")

fit_aa = read_rds("objects/fit_aa1.rds")
d1 = fit_aa$data

# code for brms model
df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

# Estimates of AA for each accession in each light treatment and light intensity
df_aa_pred1 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
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
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

df_aa_pred2 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  pivot_wider(    
    names_from = light_treatment,
    values_from = aa
  ) |>
  mutate(d_aa = high - low) |>
  group_by(light_intensity) |>
  point_interval(d_aa) 

# Estimate difference in AA between high and low light intensity
df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number())

df_aa_pred3 = posterior_epred(fit_aa, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  pivot_longer(cols = -row, names_to = "draw", values_to = "aa") |>
  full_join(df_new, by = "row") |>
  dplyr::select(-row) |>
  pivot_wider(    
    names_from = light_intensity,
    values_from = aa
  ) |>
  mutate(d_aa = `2000` - `150`) |>
  group_by(acc, light_treatment) |>
  point_interval(d_aa) 

df_aa_text = df_aa_pred1 |>
  ungroup() |>
  summarize(
    x = median(x),
    label = round(mean(aa), 3),
    aa = 0.20,
    .by = c("Growth", "Measurement")
  )

ggplot(df_aa_pred1, aes(x, aa, ymin = .lower, ymax = .upper, group = acc,
                        color = Growth, shape = Measurement)) +
  facet_grid(Measurement ~ Growth) +
  geom_pointinterval(fill = "white") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_text(data = df_aa_text, aes(x, aa, label = label), inherit.aes = FALSE) +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  ylab("amphi advantage") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none"
  )

ggsave("figures/aa-plot1.pdf", width = 6, height = 4, device = cairo_pdf, bg = "transparent")

# Predictions for values in ms
write_rds(df_aa_pred1, "objects/df_aa_pred1.rds")
write_rds(df_aa_pred2, "objects/df_aa_pred2.rds")
write_rds(df_aa_pred3, "objects/df_aa_pred3.rds")

# code I am not sure I need
df_pred2 = posterior_epred(fit_aa, newdata = df_new) |>
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
  group_by(light_intensity) |>
  arrange(d_aa) |>
  mutate(x = order(d_aa))

ggplot(df_pred2, aes(x, d_aa, ymin = .lower, ymax = .upper, group = acc)) +
  facet_grid(light_intensity ~ .) +
  geom_pointinterval() +
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

# code for Stan model
fit_aa = read_rds("objects/fit_aa5.rds")

post1 = fit_aa$draws() |>
  as_draws_df()

p1 = post1 |>
  select(
    starts_with("."),
    b0_aa,
    b_aa_light_intensity_2000,
    b_aa_light_treatment_high,
    b_aa_2000_high
  ) |>
  crossing(
    light_intensity = c("150", "2000"),
    light_treatment = c("low", "high"),
    acc_n = as.character(seq_len(29))
  )

p2 = post1 |>
  select(
    starts_with("."),
    starts_with("b_aa_light_intensity_2000_acc"),
    starts_with("b_aa_light_treatment_high_acc"),
    starts_with("b_aa_2000_high_acc"),
    matches("^b_aa_acc\\[\\d+\\]$")
  ) |>
  pivot_longer(cols = starts_with("b_aa_")) |>
  mutate(
    acc_n = str_extract(name, "(?<=\\[)\\d+(?=\\])"),
    term = str_remove(name, "\\[\\d+\\]"),
    .keep = "unused"
  ) |>
  pivot_wider(
    names_from = term,
    values_from = value
  ) 

p3 = full_join(p1, p2, by = join_by(.chain, .iteration, .draw, acc_n)) |>
  mutate(b_2000 = b_aa_light_intensity_2000 + b_aa_light_intensity_2000_acc,
         b_high = b_aa_light_treatment_high + b_aa_light_treatment_high_acc,
         b_2000_high = b_aa_2000_high + b_aa_2000_high_acc,
         .keep = "unused") |>
  mutate(aa = b0_aa + 
           b_2000 * (light_intensity == "2000") +
           b_high * (light_treatment == "high") + 
           b_2000_high * (light_intensity == "2000") * (light_treatment == "high") + b_aa_acc) |>
  group_by(acc_n, light_intensity, light_treatment) |>
  point_interval(aa)

ggplot(p3, aes(light_treatment, aa, group = acc_n)) +
  facet_wrap(~ light_intensity) +
  geom_line()

p4 = full_join(p1, p2, by = join_by(.chain, .iteration, .draw, acc_n)) |>
  mutate(b_2000 = b_aa_light_intensity_2000 + b_aa_light_intensity_2000_acc,
         b_high = b_aa_light_treatment_high + b_aa_light_treatment_high_acc,
         b_2000_high = b_aa_2000_high + b_aa_2000_high_acc,
         .keep = "unused") |>
  mutate(aa = b0_aa + 
           b_2000 * (light_intensity == "2000") +
           b_high * (light_treatment == "high") + 
           b_2000_high * (light_intensity == "2000") * (light_treatment == "high") + b_aa_acc)

p4 |>
  dplyr::select(.draw, acc_n, light_intensity, light_treatment, aa) |>
  pivot_wider(
    names_from = light_treatment,
    values_from = aa
  ) |>
  mutate(d_aa = high - low) |>
  group_by(acc_n, light_intensity) |>
  point_interval(d_aa) |>
  arrange(d_aa) |>
  print(n = Inf)
