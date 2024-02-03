# Define hyperparameters for simulate synthetic data sets
source("r/header.R")

set.seed(20240202)

# number of synthetic data sets
n_sim = 3

# hyper parameters
aa_hyperpars = list(
  
  # Experimental design hyperparameters
  n_acc = 1, # number of accessions
  n_rep = 1e1, # number of replicates per accession per treatment
  n_pts = 1e1, # number of points per curve
  
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
  mu_intercept = 0, # later on, change to mu_intercept, sigma_intercept
  mu_slope = 1,  # later on, change to mu_slope, sigma_slope
  sigma_error_intercept = 1
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
