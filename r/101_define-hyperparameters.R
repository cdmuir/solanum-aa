# SCRATCH
# 
# Get local second derivative on A-gsw curve to identify linear portion
# 
# 1. fit linear regression
# 2. simulate N synthetic data sets from fitted
# 3. Fit smooth function to each synthetic dataset
# 4. Calculate second derivative of each smooth function
# 5. Calculate second derivative cutoff based on null distribution
# Get error autocorrelation
# get_ac1 = function(x) {
#   cor(x[1:(length(x) - 1)], x[2:length(x)])
# }
# 
# rh_curves |>
#   filter(acc_id == "LA2172-AA") |>
#   reframe(resid = resid(lm(A ~ log(gsw))),
#           .by = c("acc", "acc_id", "light_treatment", "light_intensity", "curve_type")) |>
#   summarize(r = get_ac1(resid),
#             .by = c("acc", "acc_id", "light_treatment", "light_intensity", "curve_type"))
# 
# ggplot(filter(rh_curves, acc_id == "LA2172-AA"), aes(gsw, A, color = curve_type)) +
#   geom_point() +
#   scale_x_log10()
# 
# A = V * a / (K + a)
# A (K + a) = V * a 
# (A * K + A * a)/a = V
# A*K/a + A = V
# A*K/a = V - A


# Define hyperparameters for simulate synthetic data sets
source("r/header.R")

rh_curves = read_rds("data/thinned_rh_curves.rds") |>
  mutate(acc = str_replace(acc_id, id_string, "\\1")) |>
  # Focus on LA2172 as an example for now
  filter(acc == "LA2172", assumed_K == 0.5) |>
  # move this up in data prep
  mutate(leaf_type = case_when(
    curve_type == "1-sided RH" ~ "pseudohypo",
    curve_type == "2-sided RH" ~ "amphi"
  ))
fit_preliminary = read_rds("objects/fit_preliminary.rds")
post = as_draws(fit_preliminary)[[1]]

set.seed(20240202)

# Calculate realistic quantities from data
n_id = rh_curves |>
  summarise(n_id = length(unique(acc_id)), .by = "acc") |>
  pull(n_id) |>
  mean() |>
  round()

n_pts = 20 #rh_curves |>
  # summarise(n_pts = n(),
  #           .by = c("acc", "acc_id", "light_treatment", "light_intensity")) |>
  # pull(n_pts) |>
  # mean() |>
  # round()

log_gsw_range = rh_curves |>
  summarize(
    min_log_gsw = min(log(gsw)),
    max_log_gsw = max(log(gsw)),
    .by = c(
      "acc",
      "acc_id",
      "light_treatment",
      "light_intensity",
      "leaf_type"
    )
  ) |>
  summarize(
    across(min_log_gsw:max_log_gsw, mean),
    .by = c(
      "acc",
      "light_treatment",
      "light_intensity",
      "leaf_type"
    )
  )


# number of synthetic data sets
n_sim = 3

# hyper parameters
aa_hyperpars = list(
  
  # Experimental design hyperparameters
  n_acc = 1, # number of accessions
  n_id  = n_id, # number of replicates per accession per treatment
  n_pts = n_pts, # number of points per curve
  
  # Chamber environment hyperparameters
  c_a = 415, # CO2 [umol / mol]
  flow = 600, # chamber flow rate [umol / s]
  g_bw = 2.5, # boundary layer conductance to water vapor [mol / m^2 / s]
  interval = 10, # interval between measurements [s] 
  P = mean(rh_curves$Pa), # atmospheric pressure [kPa]
  RH = 0.5, # relative humidity
  s = 6, #  leaf area [cm^2]
  T_air = mean(rh_curves$Tair), # air temperature [degreeC]
  T_leaf = mean(rh_curves$Tleaf), # leaf temperature [degreeC]
  rho_error_c = runif(n_sim), # desired autocorrelation of measurement error between CO2
  rho_error_w = runif(n_sim), # desired autocorrelation of measurement error between H2O
  sigma_c = 0.1, # LI6800 IRGA SD of measurement error in CO2 [umol / mol]
  sigma_w = 0.1, # LI6800 IRGA SD of measurement error in H2O [mmol / mol]

  # Leaf hyperparameters
  K_amphi = 0.5, # stomatal conductance ratio (treating as constant, but could treat as variable),
  log_gsw_range = log_gsw_range,
  
  # intercept is amphi, light intensity = 150, light treatment = high
  mu_intercept = post$b_intercept_Intercept[seq_len(n_sim)], 
  b_intercept_pseudohypo = post$b_intercept_leaf_typepseudohypo[seq_len(n_sim)], 
  b_intercept_high_light = post$b_intercept_light_treatmenthigh[seq_len(n_sim)], 
  b_intercept_high_intensity = 
    post$b_intercept_light_intensity2000[seq_len(n_sim)], 
  `b_intercept_pseudohypo:high_intensity` =
    post$`b_intercept_leaf_typepseudohypo:light_intensity2000`[seq_len(n_sim)], 
  `b_intercept_pseudohypo:high_light` = 
    post$`b_intercept_leaf_typepseudohypo:light_treatmenthigh`[seq_len(n_sim)],
  `b_intercept_high_intensity:high_light` = 
    post$`b_intercept_light_intensity2000:light_treatmenthigh`[seq_len(n_sim)],
  `b_intercept_pseudohypo_high_intensity_high_light` = 
    post$`b_intercept_leaf_typepseudohypo:light_intensity2000:light_treatmenthigh`[seq_len(n_sim)],

  mu_slope = post$b_slope_Intercept[seq_len(n_sim)], 
  b_slope_pseudohypo = post$b_slope_leaf_typepseudohypo[seq_len(n_sim)], 
  b_slope_high_light = post$b_slope_light_treatmenthigh[seq_len(n_sim)], 
  b_slope_high_intensity = 
    post$b_slope_light_intensity2000[seq_len(n_sim)], 
  `b_slope_pseudohypo:high_intensity` =
    post$`b_slope_leaf_typepseudohypo:light_intensity2000`[seq_len(n_sim)], 
  `b_slope_pseudohypo:high_light` = 
    post$`b_slope_leaf_typepseudohypo:light_treatmenthigh`[seq_len(n_sim)],
  `b_slope_high_intensity:high_light` = 
    post$`b_slope_light_intensity2000:light_treatmenthigh`[seq_len(n_sim)],
  `b_slope_pseudohypo_high_intensity_high_light` = 
    post$`b_slope_leaf_typepseudohypo:light_intensity2000:light_treatmenthigh`[seq_len(n_sim)],
  
  mu_sigma_intercept_id = post$b_sigma_intercept_Intercept[seq_len(n_sim)],
  b_sigma_intercept_high_light_id =
    post$b_sigma_intercept_light_intensity2000[seq_len(n_sim)],

  mu_sigma_slope_id = post$b_sigma_slope_Intercept[seq_len(n_sim)],
  b_sigma_slope_high_light_id = 
    post$b_sigma_slope_light_intensity2000[seq_len(n_sim)],
  
  # No fit in model
  # day-to-day variation in intercept affecting both high and low light intensity
  sigma_intercept_error = 1
  
)

aa_hyperpars$b_autocorr_c = with(aa_hyperpars, calculate_corr_decay(rho_error_c, interval))
aa_hyperpars$b_autocorr_w = with(aa_hyperpars, calculate_corr_decay(rho_error_w, interval))

aa_hyperpars = c(aa_hyperpars,
                 aa_hyperpars$log_gsw_range |>
  # only simulating one accession at the moment
  dplyr::select(-acc) |>
  unite("x", light_treatment, light_intensity, leaf_type) |>
  pivot_longer(contains("log_gsw")) |>
  unite("name", name, x) |>
  pivot_wider() |>
  as.list()
) 

aa_hyperpars = aa_hyperpars[which(names(aa_hyperpars) != "log_gsw_range")]

lapply(aa_hyperpars, \(.x) {
  if (length(.x) == n_sim) {
    .x
  } else {
    rep(.x, n_sim)
  }
}) |>
  write_rds("objects/aa_hyperpars.rds")
