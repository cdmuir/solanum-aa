# Fit preliminary model to each RH curve to identify and trim AA outliers
source("r/header.R")

rh_curves = read_rds("data/prepared_rh_curves.rds")

rh_curves1 = rh_curves |>
  summarize(
    fit_sty = list(lm(log_A ~ poly(log_gsw, 2, raw = TRUE))),
    fit_dyn = list(lm(log_Adyn ~ poly(log_gsw, 2, raw = TRUE))),
    min_log_gsw = min(log_gsw),
    max_log_gsw = max(log_gsw),
    licor_date = first(licor_date),
    S = unique(S),
    .by = c(
      "acc",
      "acc_id",
      "light_treatment",
      "light_intensity",
      "leaf_type"
    )
  ) |>
  rowwise() |>
  mutate(
    terms_sty = list(coef(fit_sty)),
    b0_sty = terms_sty[[1]],
    b1_sty = terms_sty[[2]],
    b2_sty = terms_sty[[3]],
    terms_dyn = list(coef(fit_dyn)),
    b0_dyn = terms_dyn[[1]],
    b1_dyn = terms_dyn[[2]],
    b2_dyn = terms_dyn[[3]],
    log_gsw = min_log_gsw * (leaf_type == "amphi") + max_log_gsw * (leaf_type == "pseudohypo")
  ) |>
  dplyr::select(-starts_with("fit_"), -starts_with("terms_"), -min_log_gsw, -max_log_gsw) |>
  pivot_wider(names_from = leaf_type,
              values_from = c(b0_sty, b1_sty, b2_sty, b0_dyn, b1_dyn, b2_dyn, log_gsw, S, licor_date)) |>
  mutate(
    upper_int_sty = aa_int(
      log_gsw_pseudohypo,
      b0_sty_amphi,
      b0_sty_pseudohypo,
      b1_sty_amphi,
      b1_sty_pseudohypo,
      b2_sty_amphi,
      b2_sty_pseudohypo
    ),
    lower_int_sty = aa_int(
      log_gsw_amphi,
      b0_sty_amphi,
      b0_sty_pseudohypo,
      b1_sty_amphi,
      b1_sty_pseudohypo,
      b2_sty_amphi,
      b2_sty_pseudohypo
    ),
    upper_int_dyn = aa_int(
      log_gsw_pseudohypo,
      b0_dyn_amphi,
      b0_dyn_pseudohypo,
      b1_dyn_amphi,
      b1_dyn_pseudohypo,
      b2_dyn_amphi,
      b2_dyn_pseudohypo
    ),
    lower_int_dyn = aa_int(
      log_gsw_amphi,
      b0_dyn_amphi,
      b0_dyn_pseudohypo,
      b1_dyn_amphi,
      b1_dyn_pseudohypo,
      b2_dyn_amphi,
      b2_dyn_pseudohypo
    ),
    aa_sty = (upper_int_sty - lower_int_sty) / (log_gsw_pseudohypo - log_gsw_amphi),
    aa_dyn = (upper_int_dyn - lower_int_dyn) / (log_gsw_pseudohypo - log_gsw_amphi),
    amphi_first = ymd(licor_date_amphi) < ymd(licor_date_pseudohypo)
  )

fit_sty = lm(aa_sty ~ amphi_first + acc * light_treatment * light_intensity,
             data = rh_curves1)
fit_dyn = lm(aa_dyn ~ amphi_first + acc * light_treatment * light_intensity,
             data = rh_curves1)
summary(aov(fit_sty))
summary(aov(fit_dyn))
rh_curves1$resid_sty = rstudent(fit_sty)
rh_curves1$resid_dyn = rstudent(fit_dyn)

# Using outliers from steady state data
acc_id_outlier = rh_curves1 |>
  filter(abs(resid_sty) > aa_args$aa_outlier_threshold) |>
  dplyr::select(acc_id, aa_sty, resid_sty) |>
  pull(acc_id)

# acc_id_outlier = rh_curves1 |>
#   filter(abs(resid_dyn) > aa_args$aa_outlier_threshold) |>
#   dplyr::select(acc_id, aa_dyn, resid_dyn) |>
#   pull(acc_id)

do.call(
  "add_to_stats",
  rh_curves |>
    filter(!(acc_id %in% acc_id_outlier)) |>
    dplyr::summarize(
      n_point = length(obs),
      .by = c("acc_id", "curve_type", "light_intensity")
    ) |>
    dplyr::summarize(
      n_rh_curve7 = n(),
      n_point_per_rh_curve7 = mean(n_point)
    ) |>
    as.list()
)

rh_curves |>
  filter(!(acc_id %in% acc_id_outlier)) |>
  write_rds("data/trimmed_rh_curves.rds")
