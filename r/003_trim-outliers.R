# REVISIT THIS WITH PP CHECK
# Fit preliminary model to each RH curve to identify and trim AA outliers
source("r/header.R")

rh_curves = read_rds("data/prepared_rh_curves.rds") |>
  mutate(log_A = log(A))

rh_curves1 = rh_curves |>
  summarize(
    fit = list(lm(log_A ~ poly(log_gsw, 2, raw = TRUE))),
    min_log_gsw = min(log_gsw),
    max_log_gsw = max(log_gsw),
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
  select(-fit, -terms, -min_log_gsw, -max_log_gsw) |>
  pivot_wider(names_from = leaf_type, values_from = c(b0, b1, b2, log_gsw)) |>
  mutate(
    upper_int = aa_int(log_gsw_pseudohypo, b0_amphi, b0_pseudohypo, b1_amphi, b1_pseudohypo, b2_amphi, b2_pseudohypo),
    lower_int = aa_int(log_gsw_amphi, b0_amphi, b0_pseudohypo, b1_amphi, b1_pseudohypo, b2_amphi, b2_pseudohypo),
    aa = upper_int - lower_int
  )

fit1 = lm(aa ~ light_treatment * light_intensity * acc, data = rh_curves1)

rh_curves1$resid = rstudent(fit1)

acc_id_outlier = rh_curves1 |>
  filter(abs(resid) > aa_args$aa_outlier_threshold) |>
  select(acc_id, aa, resid) |>
  pull(acc_id)

do.call(
  "add_to_stats",
  rh_curves |>
    filter(!(acc_id %in% acc_id_outlier)) |>
    dplyr::summarize(
      n_point = length(obs),
      .by = c("acc_id", "curve_type", "light_intensity")
    ) |>
    dplyr::summarize(
      n_rh_curve6 = n(),
      n_point_per_rh_curve6 = mean(n_point)
    ) |>
    as.list()
)

rh_curves |> 
  filter(!(acc_id %in% acc_id_outlier)) |>
  write_rds("data/trimmed_rh_curves.rds")

# Writing rh_curves1 to use in preliminary analyses of AA and LMA, gmax_ratio, etc.
write_rds(rh_curves1, "objects/rh_curves1.rds")


# The rest of this is things I tried but probably won't use
# ggplot(rh_curves1, aes(light_treatment, aa, color = light_intensity)) +
#   geom_boxplot()
# 
# 
# fit1 = brm(aa ~ light_treatment * light_intensity + (1 | acc), data = rh_curves1,
#            backend = "cmdstanr", chains = 1, save_pars = save_pars(all = TRUE))
# # Estimate k-hat using the loo package
# loo1 = loo(fit1, moment_match = TRUE)
# loo1
# plot(loo1, "k")
# 
# pareto_k_values(loo1)
# pareto_k_influence_values(loo1)
# 
# rh_curves1$k_influence = pareto_k_influence_values(loo1)
# rh_curves1 |>
#   filter(k_influence > 0.7) |>
#   select(acc_id, aa)
# 
# rh_curves2 = filter(rh_curves1, acc_id != "LA1269-S")
# fit2 = brm(aa ~ light_treatment * light_intensity + (1 | acc), 
#            data = rh_curves2,
#            backend = "cmdstanr", chains = 1, save_pars = save_pars(all = TRUE))
# loo2 = loo(fit2, moment_match = TRUE)
# rh_curves2$k_hat = loo2$pointwise[, "influence_pareto_k"]
# rh_curves2 |>
#   filter(k_hat > 0.7) |>
#   select(acc_id, aa)
# plot(rh_curves2$k_hat)
# 
# # Linear model ----
# fit1 = lm(log_A ~ poly(log_gsw, 2, raw = TRUE), data = filter(rh_curves, curve == "LA0107-C_pseudohypo_150"))
# b1 = coef(fit1)
# 
# df1 = tibble(
#   log_gsw = seq(-3.25, -2.25, length.out = 100),
#   log_A = b1[1] + b1[2] * log_gsw + b1[3] * log_gsw^2
# )
# 
# int
# ggplot(
#   filter(rh_curves, curve == "LA0107-C_pseudohypo_150"),
#   aes(log_gsw, log_A)
# ) +
#   geom_point() +
#   geom_smooth(method = "lm", formula = y ~ poly(x, 2)) +
#   geom_line(data = df1, color = "red")
# 
# 
# # Bayesian approach ----
# # This model worked well for LA2172. Scaling up to all accessions
# # With all accessions, it took almost one week and still many parameters had not 
# # quite converged.
# fit_aa0 = brm(
#   bf(log_A ~ poly(log_gsw, 2) + (poly(log_gsw, 2) | curve),
#      sigma ~ 1 + (1 | curve)),
#   family = gaussian(),
#   data = rh_curves,
#   backend = "cmdstanr",
#   chains = 4,
#   cores = 4,
#   silent = 0,
#   iter = 4e3,
#   thin = 2e0,
#   max_treedepth = 12
# )
# 
# write_rds(fit_aa0, "objects/fit_aa0.rds")
# 
# # Check convergence
# summary(fit_aa0)
# 
# fit_aa0 |>
#   spread_draws(r_curve[curve, term]) |>
#   summarize(rhat = posterior::rhat(r_curve)) |>
#   dplyr::filter(rhat > 1.05)
