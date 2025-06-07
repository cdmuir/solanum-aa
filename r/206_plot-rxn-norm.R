source("r/header.R")

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
