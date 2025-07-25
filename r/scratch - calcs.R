# simplified back-of the envelope calculation for how much water could be save by switching to amphistomy
source("r/header.R")

rh_curves = read_rds("data/trimmed_rh_curves.rds")
tr = read_rds("data/phylogeny.rds")
A = vcv(tr, corr = TRUE)

fit = brm(
  log_A ~ leaf_type + light_intensity * log_gsw + (light_intensity * log_gsw | acc_id),
  data = filter(rh_curves, acc == "LA0716"),
  # data2 = list(A = A),
  backend = "cmdstanr",
  chains = 1,
  family = student()
)

summary(fit)
conditional_effects(fit)
