# Simulate synthetic data sets to validate statistical models
source("r/header.R")

set.seed(20231227)

# hyper parameters
aa_hyperpars = list(
  
  # Experimental design hyperparameters
  n_acc = 1, # number of accessions
  n_rep = 1e1, # number of replicates per accession per treatment
  n_pts = 1e1, # number of points per curve
  
  # Chamber environment hyperparameters
  c_a = 415, # CO2 [umol / mol]
  flow = 600, # chamber flow rate [umol / s]
  leaf_area = 6, # leaf area [cm^2]
  P = 100, # atmospheric pressure [kPa]
  RH = 0.5, # relative humidity
  T_air = 25, # air temperature [degreeC]
  T_leaf = 25, # leaf temperature [degreeC]
  sigma_c = 0.1, # LI6800 IRGA SD of measurement error in CO2 [umol / mol]
  sigma_w = 0.1, # LI6800 IRGA SD of measurement error in H2O [mmol / mol]
  
  # Leaf hyperparameters
  K = 0.5, # stomatal conductance ratio (treating as constant, but could treat as variable),
  g_bw = 2.57, # boundary layer conductance. assume the same for both surfaces, all accessions
  min_gsw_amphi = 0.10,
  max_gsw_amphi = 0.50,
  min_gsw_pseudohypo = 0.05,
  max_gsw_pseudohypo = 0.25,
  intercept = 0, # later on, change to mu_intercept, sigma_intercept
  slope = 1,  # later on, change to mu_slope, sigma_slope
  rho_error_resid = 0.9, # correlation coefficient between data points (need to adjust to per second rate)
  sigma_error_intercept = 1,
  sigma_error_resid = 0.1
)


rep_vector = with(aa_hyperpars, str_c("r", str_pad(
  seq_len(n_rep), floor(log10(n_rep)) + 1, "left", "0"
)))
pts_vector = with(aa_hyperpars, str_c("p", str_pad(
  seq_len(n_pts), floor(log10(n_pts)) + 1, "left", "0"
)))
leaf_type_vector = c("amphi", "pseudohypo")


aa_pars = with(
  aa_hyperpars,
  crossing(nesting(
    crossing(rep = rep_vector,
             leaf_type = leaf_type_vector),
    error_intercept = rnorm(2 * n_rep, 0, sigma_error_intercept)
  ),
  pts = pts_vector) |>
    split( ~ rep + leaf_type) |>
    # this will need to be adjusted for time interval between points
    map_dfr(\(x) mutate(x, error_resid = arima.sim(
      list(ar = rho_error_resid), n = nrow(x), sd = sigma_error_resid
    ), c_a = c_a))
)

df_sim = with(
  aa_hyperpars,
  crossing(
    rep = rep_vector,
    leaf_type = leaf_type_vector,
    nesting(
      pts = pts_vector,
      i = seq(0, 1, length.out = n_pts)
    )
  ) |>
    full_join(aa_pars, by = join_by(rep, pts, leaf_type)) |>
    mutate(
      min_log_gsw = (leaf_type == "amphi") * log(min_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(min_gsw_pseudohypo),
      max_log_gsw = (leaf_type == "amphi") * log(max_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(max_gsw_pseudohypo),
      log_gsw_real = min_log_gsw + (max_log_gsw - min_log_gsw) * i,
      log_A_real = intercept + error_intercept + slope * log_gsw_real,
      gsw_real = exp(log_gsw_real),
      A_real = exp(log_A_real),
      K = K * 1 / (leaf_type == "amphi"),
))

df_sim
# scratch
leak = 1 # correction factor, AU
CO2_s = 415.775 # CW
CO2_r = 429.717 # CX
H2O_s = 18.2744 # CY
H2O_r = 15.899 # CZ
flow = 600.012 # DC19
leaf_area = 6 # $B$7
T_leaf = 27.90632 # V, using EB
T_air = 26.4793 # DF
P = 81.2318 # DD
w_i = 0.61365 * exp(17.502 * T_air / (240.97 + T_air)) # Do I want T_air or T_leaf here?

E_area = 1000 * flow * leak * (H2O_s - H2O_r) / (100 * leaf_area * (1000 - leak * H2O_s))
A = flow * leak * (CO2_r - CO2_s * (1000 - leak * H2O_r) / (1000 - leak * H2O_s)) / (100 * leaf_area)

# Solve for H2O_s if E_area is known
E_area = 1000 * flow * (H2O_s - H2O_r) / (100 * leaf_area * (1000 - H2O_s))
E_area * (100 * leaf_area * (1000 - H2O_s)) = 1000 * flow * (H2O_s - H2O_r)
E_area * 100 * leaf_area * 1000 - E_area * 100 * leaf_area * H2O_s = 1000 * flow * H2O_s - 1000 * flow * H2O_r
E_area * 100 * leaf_area * 1000 + 1000 * flow * H2O_r = 1000 * flow * H2O_s + E_area * 100 * leaf_area * H2O_s
E_area * 100 * leaf_area * 1000 + 1000 * flow * H2O_r = H2O_s * (1000 * flow + E_area * 100 * leaf_area)
(E_area * 100 * leaf_area * 1000 + 1000 * flow * H2O_r) / (1000 * flow + E_area * 100 * leaf_area) = H2O_s
(E_area * leaf_area * 1000 + 10 * flow * H2O_r) / (10 * flow + E_area * leaf_area) = H2O_s

RH = 0.428
C_a = 429.717 # C_a should be same as CO2_r
K = 0.5
g_bw = 2.57
A_real = 13.1
w_a = RH * w_i
g_tw = if_else(
  is.finite(K), 
  1 / ((1 + K) * (1 / gsw_real) + (1 / g_bw)) + K / ((1 + K) * (1 / gsw_real) + K / g_bw),
  1 / (1 / gsw_real + 1 / g_bw)
)
g_tc = if_else(
  is.finite(K), 
  1 / ((1 + K) * (1.6 / gsw_real) + (1.37 / g_bw)) + K / ((1 + K) * (1.6 / gsw_real) + K * 1.37 / g_bw),
  1 / (1.6 / gsw_real + 1.37 / g_bw)
)
E_area = g_tw * (w_i - w_a) # Transpiration per area [mmol / m^2 / s]
E = E_area * leaf_area / 1e4 # [mmol / s]
# True water vapour concentration in sample IRGA [mmol / mol]
w_0 = 1 - E * (1 - w_a) / (flow / 1e6)
  

u_a = u_0 + s * E_area
s * E_area = u_a * w_a - u_0 * w_0

s * E_area = u_0 + s * E_area * w_a - u_0 * w_0
s * E_area + u_0 * w_0 = u_0 + s * E_area * w_a 
E + u_0 * w_0 = u_0 + E * w_a 
u_0 * w_0 = u_0 + E * w_a - E
u_0 * w_0 = u_0 - E * (1 - w_a)
w_0 = 1 - E * (1 - w_a) / u_0

# input
# A : CO2 assimilation rate [umol / m^2 / s]
# c_a : [CO2] in [umol / mol] at sample IRGA
# flow : chamber flow rate [umol / s]
# g_bw : boundary layer conductance to water vapor [mol / m^2 / s]
# g_sw : stomatal conductance to water vapor [mol / m^2 / s]
# K : stomatal conductance ratio
# P : chamber pressure [kPa]
# RH : relative humidity
# s : leaf area [cm^2]
# T_air : air temperature [C] in chamber
# T_leaf : leaf temperature [C]
A = 12.93624 # L
c_a = 415.775 # CW
flow = 600.012 # DC
g_bw = 2.577362 # R
g_sw = 0.084753 # Q (using calculated value to check derivation, but usually this would be assumed)
K = 0.5 # $C$7
P = 81.2318 + 0.099916 # DD + DE
RH = 0.428183 # Y
s = 6 # $B$7
T_air = 26.4793 # DF
T_leaf = 27.90632 # V, using EB

# intermediate
# E : transpiration per area [mmol / m^2 / s]
# g_tC : total conductance to CO2 [mol / m^2 / s]
# g_tw : total conductance to water vapor [mol / m^2 / s]
# w_i : internal leaf [H2O] in [mmol / mol] (assuming saturation)
# w_sat : saturating chamber [H2O] in [mmol / mol]
w_i = 1000 * 0.61365 * exp(17.502 * T_leaf / (240.97 + T_leaf)) / P
w_sat = 1000 * 0.61365 * exp(17.502 * T_air / (240.97 + T_air)) / P
w_a = RH * w_sat # target = 18.2744
g_tc = if_else(
  is.finite(K), 
  1 / ((1 + K) * (1.6 / g_sw) + (1.37 / g_bw)) + K / ((1 + K) * (1.6 / g_sw) + K * 1.37 / g_bw),
  1 / (1.6 / g_sw + 1.37 / g_bw)
) # T
g_tw = if_else(
  is.finite(K), 
  1 / ((1 + K) * (1 / g_sw) + (1 / g_bw)) + K / ((1 + K) * (1 / g_sw) + K / g_bw),
  1 / (1 / g_sw + 1 / g_bw)
) # S

E = 1000 * g_tw * (w_i - w_a) / (1000 - (w_i + w_a) / 2) # K
E # target = 2.419665

# output
# w_a : [H2O] in [mmol / mol] at sample IRGA (intermediate calculation)
# w_0 : [H2O] in [mmol / mol] at reference IRGA
# c_a : [CO2] in [umol / mol] at sample IRGA (assumed)
# c_0 : [CO2] in [umol / mol] at reference IRGA

w_0 = w_a - (s * E * 100) * (1000 - w_a) / (1000 * flow) # target = 15.899, NEED TO MULTIPLY BY SCALARS TO GET CORRECT UNITS
w_0 # target = 15.899
c_0 = 1000 * (s * A * 100) / (1000 * flow) + c_a * ((1000 - w_0) / (1000 - w_a)) # CX
c_0 # target = 429.717

.vars = tibble(
  A = 12.93624, # L
  c_a = 415.775, # CW
  flow = 600.012, # DC
  g_bw = 2.577362, # R
  g_sw = 0.084753, # Q (using calculated value to check derivation, but usually this would be assumed)
  K = 0.5, # $C$7
  P = 81.2318 + 0.099916, # DD + DE
  RH = 0.428183, # Y
  s = 6, # $B$7
  T_air = 26.4793, # DF
  T_leaf = 27.90632 # V, using EB
  
)







ggplot(df_sim, aes(log_gsw, log_A, color = leaf_type)) +
  facet_wrap(~rep) +
  geom_point()

write_rds(df_sim, "synthetic-data/df_sim.rds")
