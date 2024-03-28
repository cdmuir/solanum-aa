# Prepare initial values for solanum-aa models
source("r/header.R")

rh_curves = read_rds("data/trimmed_rh_curves.rds")

b = rh_curves |>
  mutate(across(c("curve"), \(.x) as.numeric(as.factor(.x)))) |>
  split(~ curve) |>
  map_dfr(\(.x) {
    fit = lm(log(A) ~ scaled_log_gsw + I(scaled_log_gsw ^ 2), data = .x)
    b = unname(coef(fit))
    tibble(b0 = b[1], b1 = b[2], b2 = b[3])
  })

Mu_curve = apply(b, 2, mean)

B_curve = b |>
  mutate(b0 = b0 - Mu_curve[1], b1 = b1 - Mu_curve[2], b2 = b2 - Mu_curve[3])

R_curve = cor(B_curve)

log_sigma_curve = log(diag(cov(B_curve)))

init = list(
  Mu_curve = Mu_curve,
  R_curve = R_curve,
  log_sigma_curve = log_sigma_curve,
  B_curve = B_curve
)

write_rds(init, "objects/init.rds")
