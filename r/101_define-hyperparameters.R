# SCRATCH
# 
# Get local secoond derivative on A-gsw curve to identify linear portion
# 
# 1. fit linear regression
# 2. simulate N synthetic data sets from fitted
# 3. Fit smooth function to each synthetic dataset
# 4. Calculate second derivative of each smooth function
# 5. Calculate second derivative cutoff based on null distribution
# Get error autocorrelation
get_ac1 = function(x) {
  cor(x[1:(length(x) - 1)], x[2:length(x)])
}

rh_curves |>
  filter(acc_id == "LA2172-AA") |>
  reframe(resid = resid(lm(A ~ log(gsw))),
          .by = c("acc", "acc_id", "light_treatment", "light_intensity", "curve_type")) |>
  summarize(r = get_ac1(resid),
            .by = c("acc", "acc_id", "light_treatment", "light_intensity", "curve_type"))

ggplot(filter(rh_curves, acc_id == "LA2172-AA"), aes(gsw, A, color = curve_type)) +
  geom_point() +
  scale_x_log10()

A = V * a / (K + a)
A (K + a) = V * a 
(A * K + A * a)/a = V
A*K/a + A = V
A*K/a = V - A


# Define hyperparameters for simulate synthetic data sets
source("r/header.R")

rh_curves = read_rds("data/rh_curves.rds") |>
  mutate(acc = str_replace(acc_id, id_string, "\\1")) |>
  # Focus on LA2172 as an example for now
  filter(acc == "LA2172", assumed_K == 0.5)
fit_preliminary = read_rds("objects/fit_preliminary.rds")

set.seed(20240202)

# Calculate realistic quantities from data
n_id = rh_curves |>
  summarise(n_id = length(unique(acc_id)),
            .by = c("acc", "light_treatment")) |>
  pull(n_id) |>
  mean() |>
  round()

n_pts = rh_curves |>
  summarise(n_pts = n(),
            .by = c("acc", "acc_id", "light_treatment", "light_intensity")) |>
  pull(n_pts) |>
  mean() |>
  round()



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
  P = 100, # atmospheric pressure [kPa]
  RH = 0.5, # relative humidity
  s = 6, #  leaf area [cm^2]
  T_air = 25, # air temperature [degreeC]
  T_leaf = 25, # leaf temperature [degreeC]
  rho_error_c = runif(n_sim), # desired autocorrelation of measurement error between CO2
  rho_error_w = runif(n_sim), # desired autocorrelation of measurement error between H2O
  sigma_c = 0.1, # LI6800 IRGA SD of measurement error in CO2 [umol / mol]
  sigma_w = 0.1, # LI6800 IRGA SD of measurement error in H2O [mmol / mol]

  # Leaf hyperparameters
  K_amphi = 0.5, # stomatal conductance ratio (treating as constant, but could treat as variable),
  min_gsw_amphi = 0.10,
  max_gsw_amphi = 0.50,
  min_gsw_pseudohypo = 0.05,
  max_gsw_pseudohypo = 0.25,
  
  # intercept is amphi, light intensity = 150, light treatment = low
  mu_intercept = 40, 
  mu_intercept_low_light = -20, 
  sigma_intercept_id = 2,
  sigma_intercept_error = 1,
  sigma_intercept_low_light_id = 3,
  
  mu_slope = 9, 
  mu_slope_low_light = -3, 
  sigma_slope_id = 1,
  sigma_slope_low_light_id = 1
)

# Calculate the decay, per s, of autocorrelation between data points given the
# desired autocorrelation and interval between data points
aa_hyperpars$b_autocorr_c = with(aa_hyperpars, calculate_corr_decay(rho_error_c, interval))
aa_hyperpars$b_autocorr_w = with(aa_hyperpars, calculate_corr_decay(rho_error_w, interval))

lapply(aa_hyperpars, \(.x) {
  if (length(.x) == n_sim) {
    .x
  } else {
    rep(.x, n_sim)
  }
}) |>
  write_rds("objects/aa_hyperpars.rds")
