# Preliminary fit to inform simulations
source("r/header.R")
rh_curves = read_rds("data/rh_curves.rds") |>
  separate_wider_delim(acc_id, "-", names = c("acc", "id"), cols_remove = FALSE)
  
# Example curve to test with
tmp = rh_curves |>
  filter(acc_id == "LA1777-E", assumed_K == 0.5)

# M-M not working
pars = nls(
  A ~ V * gsw / (Km + gsw),
  data = filter(
    tmp,
    curve_type == "1-sided RH",
    light_intensity == "150",
    assumed_K == 0.5
  ),
  start = list(V = 50, Km = 0.3)
)

df_line = tibble(
  gsw = seq(0, 0.4, 0.01),
  A = coef(pars)["V"] * gsw / (coef(pars)["Km"] + gsw),
  curve_type = "1-sided RH",
  light_intensity = "150"
)

ggplot(tmp, aes(gsw, A, color = curve_type)) +
  geom_point() +
  # scale_x_log10() +
  geom_line(data = df_line)

# Function to calculate a test statistic for whether there is a nonlinear
# relationship between fitted values and residuals
test_stat = function(formula) {
  fit = lm(formula)
  dat = tibble(p = predict(fit), r = resid(fit))
  summary(lm(r ~ poly(p, 2), data  = dat))$r.squared
}

# Conclusion 1: A-log(gsw) is best, but still problematic, esp for some at 2000 PAR ----
df1 = rh_curves |>
  filter(assumed_K == 0.5, !is.na(gsw)) |>
  mutate(rsw = 1 / gsw) |>
  summarise(
    A_gsw = test_stat(A ~ gsw),
    logA_gsw = test_stat(log(A) ~ gsw),
    A_loggsw = test_stat(A ~ log(gsw)),
    logA_loggsw = test_stat(log(A) ~ log(gsw)),
    A_rsw = test_stat(A ~ rsw),
    logA_rsw = test_stat(log(A) ~ rsw),
    A_logrsw = test_stat(A ~ log(rsw)),
    logA_logrsw = test_stat(log(A) ~ log(rsw)),
    .by = c("curve_type", "light_intensity", "acc_id")
  ) 

df1 |>
  pivot_longer(A_gsw:logA_logrsw) |>
  ggplot(aes(value, color = name, fill = name)) +
  facet_wrap(light_intensity ~ name) +
  geom_density(alpha = 0.5)

# Conclusion 2: Is M-M better?


df1 = rh_curves |>
  filter(assumed_K == 0.5, !is.na(gsw)) |>
  summarise(
    # fit1 = list(try(nls(A ~ V * gsw / (Km + gsw)))),
    fit2 = list(lm(A ~ log(gsw))),
    keep = inherits(fit1[[1]], "nls"),
    .by = c("curve_type", "light_intensity", "acc_id")
  ) |>
  filter(keep)

df1$aic1 = sapply(seq_len(nrow(df1)), \(.i) AIC(df1$fit1[.i][[1]]))
df1$aic2 = sapply(seq_len(nrow(df1)), \(.i) AIC(df1$fit2[.i][[1]]))

# 2468 of 2566 fit MM
df1 |>
  # Negative: AIC1 > AIC2 (lm favored over nls)
  # Positive: AIC2 > AIC1 (nls favored over lm)
  mutate(dAIC = aic2 - aic1) |>
  # summarize(dAIC = median(dAIC), .by = c("curve_type", "light_intensity"))
  ggplot(aes(dAIC)) +
  facet_grid(rows = vars(curve_type), cols = vars(light_intensity)) +
  geom_histogram() +
  geom_vline(xintercept = 0)


# OLD
# fit1 = brm(A ~ curve_type + light_treatment + light_intensity + log(gsw) +
#              curve_type:light_treatment + curve_type:light_intensity +
#              light_treatment:log(gsw) + light_intensity:log(gsw) +
#              (curve_type + light_treatment + light_intensity + log(gsw) | acc) +
#              (log(gsw) | acc_id), data = rh_curves,
#           chains = 1, backend = "cmdstanr", silent = 0)
# 
# fit1$save_objects("objects/fit1-preliminary.rds")