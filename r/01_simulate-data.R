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
  g_bw = 2.5, # boundary layer conductance to water vapor [mol / m^2 / s]
  P = 100, # atmospheric pressure [kPa]
  RH = 0.5, # relative humidity
  s = 6, #  leaf area [cm^2]
  T_air = 25, # air temperature [degreeC]
  T_leaf = 25, # leaf temperature [degreeC]
  sigma_c = 0.1, # LI6800 IRGA SD of measurement error in CO2 [umol / mol]
  sigma_w = 0.1, # LI6800 IRGA SD of measurement error in H2O [mmol / mol]

  # Leaf hyperparameters
  K_amphi = 0.5, # stomatal conductance ratio (treating as constant, but could treat as variable),
  g_bw = 2.5, # boundary layer conductance to water vapor [mol / m^2 / s]
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
    split(~ rep + leaf_type) |>
    # this will need to be adjusted for time interval between points
    map_dfr(\(x) {
      x |>
        mutate(
          error_resid = arima.sim(list(ar = rho_error_resid), n = nrow(x), sd = sigma_error_resid),
          c_a = c_a,
          flow = flow,
          g_bw = g_bw,
          P = P,
          RH = RH,
          s = s,
          T_air = T_air,
          T_leaf = T_leaf,
          sigma_c = sigma_c,
          sigma_w = sigma_w 
        )
    })
)


df_sim = with(
  aa_hyperpars,
  crossing(
    rep = rep_vector,
    leaf_type = leaf_type_vector,
    nesting(pts = pts_vector,
            i = seq(0, 1, length.out = n_pts))
  ) |>
    full_join(aa_pars, by = join_by(rep, pts, leaf_type)) |>
    mutate(
      min_log_gsw = (leaf_type == "amphi") * log(min_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(min_gsw_pseudohypo),
      max_log_gsw = (leaf_type == "amphi") * log(max_gsw_amphi) +
        (leaf_type == "pseudohypo") * log(max_gsw_pseudohypo),
      log_gsw_real = min_log_gsw + (max_log_gsw - min_log_gsw) * i,
      log_A_real = intercept + error_intercept + slope * log_gsw_real,
      g_sw = exp(log_gsw_real),
      A = exp(log_A_real),
      K = K_amphi * 1 / (leaf_type == "amphi")
    ) 
)

df_sim = df_sim |>
  li6800_simulate() |>
  # NEED TO ADD IN TIME-CORRELATED ERROR USING ARIMA
  li6800_add_error(select(df_sim, sigma_c, sigma_w)) |>
  li6800_estimate() 

ggplot(df_sim, aes(g_sw, A, color = leaf_type)) +
  facet_wrap(~rep) +
  geom_point()

ggplot(df_sim, aes(gsw_hat, A_hat, color = leaf_type)) +
  facet_wrap(~rep) +
  geom_point()

write_rds(df_sim, "synthetic-data/df_sim.rds")
