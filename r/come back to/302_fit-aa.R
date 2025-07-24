# Estimate AA for HI Solanum
source("r/header.R")

rh_curves = read_rds("data/prepared_rh_hi_curves.rds") |>
  mutate(log_A = log(A))

rh_curves1 = rh_curves |>
  summarize(
    fit = list(lm(log_A ~ poly(log_gsw, 2, raw = TRUE))),
    min_log_gsw = min(log_gsw),
    max_log_gsw = max(log_gsw),
    licor_date = first(licor_date),
    S = unique(S),
    .by = c("acc", "acc_id", "light_treatment", "light_intensity", "leaf_type")
  ) |>
  rowwise() |>
  mutate(
    terms = list(coef(fit)),
    b0 = terms[[1]],
    b1 = terms[[2]],
    b2 = terms[[3]],
    log_gsw = min_log_gsw * (leaf_type == "amphi") + max_log_gsw * (leaf_type == "pseudohypo")
  ) |>
  dplyr::select(-fit, -terms, -min_log_gsw, -max_log_gsw) |>
  pivot_wider(names_from = leaf_type, values_from = c(b0, b1, b2, log_gsw, S, licor_date)) |>
  mutate(
    upper_int = aa_int(log_gsw_pseudohypo, b0_amphi, b0_pseudohypo, b1_amphi, b1_pseudohypo, b2_amphi, b2_pseudohypo),
    lower_int = aa_int(log_gsw_amphi, b0_amphi, b0_pseudohypo, b1_amphi, b1_pseudohypo, b2_amphi, b2_pseudohypo),
    aa = (upper_int - lower_int) / (log_gsw_pseudohypo - log_gsw_amphi),
    amphi_first = ymd(licor_date_amphi) < ymd(licor_date_pseudohypo)
  ) 

fit1 = lm(aa ~ amphi_first + acc * light_treatment * light_intensity, data = rh_curves1)
summary(aov(fit1))

ggplot(filter(rh_curves1, light_intensity == "2000"), aes(x = light_treatment, y = aa)) +
  geom_point(position = position_jitter(width = 0.1)) +
  facet_wrap(~ acc) 
