# Plot PAI against AA and LMA
source("r/header.R")

fit_aa1 = read_rds("objects/fit_aa1_sty.rds")
fit_aa2 = read_rds("objects/fit_aa2.rds")
d1 = fit_aa1$data
accession_gedi = read_rds("data/accession-gedi.rds") |>
  rename(acc = accession) |>
  filter(acc %in% unique(d1$acc))

df_new = crossing(
  acc = unique(d1$acc),
  light_treatment = c("low", "high"),
  light_intensity = c("150", "2000")
) |>
  mutate(row = row_number(), se_aa = 0)

df_aa_pred1 = posterior_epred(fit_aa1, newdata = df_new) |>
  t() |>
  as_draws_df() |>
  mutate(row = row_number()) |>
  dplyr::select(-.chain, -.iteration, -.draw) |>
  pivot_longer(cols = -row,
               names_to = "draw",
               values_to = "aa") |>
  full_join(df_new, by = "row")

df_coef = df_aa_pred1 |>
  left_join(accession_gedi, by = join_by(acc)) |>
  split( ~ draw + light_treatment + light_intensity) |>
  map(\(.x) {
    lm(aa ~ log(pai), data = .x)
  }) |>
  map(coef) |>
  map(as_tibble) |>
  map(mutate, term = c("intercept", "slope")) |>
  map(pivot_wider, names_from = term) |>
  imap_dfr(\(.x, .y) {
    mutate(
      .x,
      light_treatment = str_extract(.y, "low|high"),
      light_intensity = str_extract(.y, "150|2000")
    )
  })  |>
  mutate(draw = row_number())

df_coef_summary = df_coef |>
  pivot_longer(intercept:slope, names_to = "term", values_to = "value") |>
  group_by(term, light_treatment, light_intensity) |>
  point_interval(value)

df_aa_pred2 = df_aa_pred1 |>
  group_by(acc, light_treatment, light_intensity) |>
  point_interval(aa) |>
  left_join(accession_gedi, by = join_by(acc)) |>
  mutate(x = order(aa)) |>
  refactor_for_figure()

df_lineribbon = df_coef |>
  crossing(log_pai = log(seq(
    min(accession_gedi$pai),
    max(accession_gedi$pai),
    length.out = 10
  ))) |>
  mutate(aa = intercept + slope * log_pai) |>
  group_by(log_pai, light_treatment, light_intensity) |>
  point_interval(aa) |>
  mutate(pai = exp(log_pai)) |>
  refactor_for_figure()

df_text = df_coef_summary |>
  filter(term == "slope") |>
  refactor_for_figure() |>
  mutate(
    # across(value:.upper, \(.x) sprintf("%.3f", -.x)),
    # label = glue("-{term}: {value}\n95% CI: {.upper} – {.lower}"),
    label = case_when(
      sign(.lower) == sign(.upper) ~ "*",
      sign(.lower) != sign(.upper) ~ "n.s."
    ),
    pai = max(accession_gedi$pai),
    aa = 0.16
  )

fig_aa_pai = ggplot(df_aa_pred2, aes(pai, aa, color = Growth, shape = Measurement)) +
  facet_grid(Measurement ~ Growth) +
  geom_ribbon(
    data = df_lineribbon,
    aes(ymin = .lower, ymax = .upper),
    alpha = 0.2,
    linetype = "dashed"
  ) +
  geom_line(data = df_lineribbon) +
  geom_point() +
  geom_text(
    data = df_text,
    mapping = aes(label = label, size = label),
    hjust = 0,
    vjust = 1,
    color = "black"
  ) +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_shape_manual(values = c("low" = 19, "high" = 21)) +
  scale_size_manual(values = c("*" = 5, "n.s." = 3), guide = "none") +
  scale_x_continuous(breaks = c(0.01, 0.1, 1), trans = reverselog10_trans()) +
  xlab(expression(paste("native plant area index [", m^2 ~ m^-2, "]"))) +
  ylab("amphi advantage") +
  ylim(0, 0.16) +
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

write_rds(df_coef_summary, "objects/df_coef_summary.rds")

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
  mutate(Growth = light_treatment |>
           factor(levels = c("low", "high")) |>
           fct_recode(sun = "high", shade = "low"))

ggplot(df_aa_pred2, aes(pai, exp(llma), color = Growth)) +
  geom_point() +
  scale_color_manual(values = c("shade" = "tomato4", "sun" = "tomato")) +
  scale_x_continuous(breaks = c(0.01, 0.1, 1), trans = reverselog10_trans()) +
  xlab(expression(paste("plant area index [", m^2 ~ m^-2, "]"))) +
  ylab(expression(paste("leaf mass per area [g ", m^-2, "]")))

ggsave(
  "figures/pai-lma.pdf",
  width = 6,
  height = 4,
  device = cairo_pdf,
  bg = "transparent"
)
