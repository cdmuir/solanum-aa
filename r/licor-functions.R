# Estimate parameters from simualted IRGA measurements with error following
# LI6800 manual
li6800_estimate = function(.sims, check = TRUE, ...) {
  # input
  # .sims : output from li6800_add_error()
  # flow : chamber flow rate [umol / s]
  # g_bw : boundary layer conductance to water vapor [mol / m^2 / s]
  # K : stomatal conductance ratio
  # P : chamber pressure [kPa]
  # RH : relative humidity
  # s : leaf area [cm^2]
  # T_air : air temperature [C] in chamber
  # T_leaf : leaf temperature [C]
  #
  # .sims : output from li6800_simulate()
  # CO2_r : apparent [CO2] in [umol / mol] at reference IRGA
  # CO2_s : apparent [CO2] in [umol / mol] at sample IRGA
  # H2O_r : apparent [H2O] in [mmol / mol] at reference IRGA
  # H2O_s : apparent [H2O] in [mmol / mol] at sample IRGA

  assert_data_frame(.sims)
  assert_flag(check)
  if (check) check_li6800_estimate_sims(.sims)

  # output
  # A_hat : estimated CO2 assimilation rate [umol / m^2 / s]

  .sims |>
    mutate(
      A_hat = li6800_A_hat(CO2_r, CO2_s, H2O_r, H2O_s, flow, s),
      E_hat = li6800_E_hat(H2O_r, H2O_s, flow, s),
      gtw_hat = li6800_gtw_hat(E, H2O_s, P, T_leaf),
      gsw_hat = li6800_gsw_hat(g_bw, g_tw, K)
    )
  
}

check_li6800_estimate_sims = function(.sims) {
  check_li6800_simulate_vars(.sims)
  assert_subset(names(.sims),
                choices = c("CO2_r",
                            "CO2_s",
                            "H2O_r",
                            "H2O_s"))
  assert_numeric(
    .sims$CO2_r,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .sims$CO2_s,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .sims$H2O_r,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .sims$H2O_s,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  
  invisible(.sims)
  
}

# Estimate A [umol / m^2 / s] from simulated LI6800 data following the LI6800
# manual
li6800_A_hat = function(CO2_r, CO2_s, H2O_r, H2O_s, flow, s) {
  flow * (CO2_r - CO2_s * ((1000 - H2O_r) / (1000 - H2O_s))) / (100 * s)
}

# Estimate E [mmol / m^2 / s] from simulated LI6800 data following the LI6800
# manual
li6800_E_hat = function(H2O_r, H2O_s, flow, s) {
  1000 * flow * (H2O_s - H2O_r) / (100 * s * (1000 - H2O_s))
}

# Estimate g_tw [mol / m^2 / s] from simulated LI6800 data following the LI6800
# manual
li6800_gtw_hat = function(E, H2O_s, P, T_leaf) {
  svp = li6800_svp(T_leaf, P)
  E / 1000 * (1000 - (svp + H2O_s) / 2) / (svp - H2O_s)
}

# Estimate g_sw [mol / m^2 / s] from simulated LI6800 data following the LI6800
# manual
li6800_gsw_hat = function(g_bw, g_tw, K) {
  r_bw = 1 / g_bw
  r_tw = 1 / g_tw
  if_else(is.finite(K),
          2 / ((r_tw - r_bw) + sign(g_tw) * sqrt((r_tw - r_bw) *
                                                   (r_tw - r_bw) + 4 * K / ((K + 1) * (K + 1)) * (2 * r_tw * 1 /
                                                                                                    g_bw - r_bw * r_bw)
          )),
          2 / ((r_tw - r_bw) + sign(g_tw) * sqrt((r_tw - r_bw) *
                                                   (r_tw - r_bw))))
}

# Add IRGA measurement error to simulated LI6800 data
li6800_add_error = function(.sims, .vars, check = TRUE, ...) {
  # input
  # .sims : output from li6800_simulate()
  # c_0 : [CO2] in [umol / mol] at reference IRGA
  # c_a : [CO2] in [umol / mol] at sample IRGA
  # w_0 : [H2O] in [mmol / mol] at reference IRGA
  # w_a : [H2O] in [mmol / mol] at sample IRGA
  #
  # .vars : variables defining error
  # sigma_c : standard deviation of error in [CO2] umol / mol
  # sigma_w : standard deviation of error in [H2O] mmol / mol
  
  assert_data_frame(.sims)
  assert_data_frame(.vars)
  assert_flag(check)
  if (check) {
    check_li6800_add_error_sims(.sims)
    check_li6800_add_error_vars(.vars)
  }
  n = nrow(.sims)
  .sims |>
    mutate(
      CO2_r = c_0 + rnorm(n, 0, .vars$sigma_c),
      CO2_s = c_a + rnorm(n, 0, .vars$sigma_c),
      H2O_r = w_0 + rnorm(n, 0, .vars$sigma_w),
      H2O_s = w_a + rnorm(n, 0, .vars$sigma_w)
    )
  
}

check_li6800_add_error_sims = function(.sims) {
  assert_subset(names(.sims),
                choices = c("c_0",
                            "c_a",
                            "w_0",
                            "w_a"))
  assert_numeric(
    .sims$c_0,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .sims$c_a,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .sims$w_0,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .sims$w_a,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  
  invisible(.sims)
  
}

check_li6800_add_error_vars = function(.vars) {
  assert_subset(names(.vars),
                choices = c("sigma_c",
                            "sigma_w"))
  assert_numeric(
    .vars$sigma_c,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$sigma_w,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  
  invisible(.vars)
  
}

# Simulate [H2O] and [CO2] at IRGAs in LI6800
li6800_simulate = function(.vars, check = TRUE, ...) {
  
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
  assert_data_frame(.vars)
  assert_flag(check)
  
  if (check) check_simulate_li6800_vars(.vars)
  
  # intermediate
  # E : transpiration per area [mmol / m^2 / s]
  # g_tC : total conductance to CO2 [mol / m^2 / s]
  # g_tw : total conductance to water vapor [mol / m^2 / s]
  # w_i : internal leaf [H2O] in [mmol / mol] (assuming saturation)
  # w_sat : saturating chamber [H2O] in [mmol / mol]

  # output
  # w_a : [H2O] in [mmol / mol] at sample IRGA (intermediate calculation)
  # w_0 : [H2O] in [mmol / mol] at reference IRGA
  # c_a : [CO2] in [umol / mol] at sample IRGA (assumed)
  # c_0 : [CO2] in [umol / mol] at reference IRGA
  
  .vars |>
    mutate(
      w_i = li6800_svp(T_leaf, P),
      w_a = RH * li6800_svp(T_air, P),
      g_tc = li6800_gtc(g_bw, g_sw, K),
      g_tw = li6800_gtw(g_bw, g_sw, K),
      E = li6800_E(g_tw, w_a, w_i),
      w_0 = li6800_w0(E, flow, s, w_a),
      c_0 = li6800_c0(A, c_a, flow, s, w_0)
    )
  
}

check_li6800_simulate_vars = function(.vars) {
  assert_subset(
    names(.vars),
    choices = c(
      "A",
      "c_a",
      "flow",
      "g_bw",
      "g_sw",
      "K",
      "P",
      "RH",
      "s",
      "T_air",
      "T_leaf"
    )
  )
  assert_numeric(.vars$A,
                 len = 1L,
                 finite = TRUE,
                 any.missing = FALSE)
  assert_numeric(
    .vars$c_a,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$flow,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$g_bw,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$g_sw,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$K,
    len = 1L,
    finite = FALSE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$P,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0
  )
  assert_numeric(
    .vars$RH,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0,
    upper = 1
  )
  assert_numeric(
    .vars$s,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0,
    upper = 6
  )
  assert_numeric(
    .vars$T_air,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0,
    upper = 50
  )
  assert_numeric(
    .vars$T_leaf,
    len = 1L,
    finite = TRUE,
    any.missing = FALSE,
    lower = 0,
    upper = 50
  )
  
  invisible(.vars)
  
}

# Calculate saturating vapor pressure following the LI6800 manual
li6800_svp = function(T_degreeC, P_kPa) {
  1000 * 0.61365 * exp(17.502 * T_degreeC / (240.97 + T_degreeC)) / P_kPa
}

# Calculate total conductance to CO2 following the LI6800 manual
li6800_gtc = function(g_bw, g_sw, K) {
  if_else(is.finite(K),
          1 / ((1 + K) * (1.6 / g_sw) + (1.37 / g_bw)) + K / ((1 + K) * (1.6 / g_sw) + K * 1.37 / g_bw),
          1 / (1.6 / g_sw + 1.37 / g_bw))
}

# Calculate total conductance to water vapor following the LI6800 manual
li6800_gtw = function(g_bw, g_sw, K) {
  r_bw = 1 / g_bw
  r_sw = 1 / g_sw
  if_else(is.finite(K),
          1 / ((1 + K) * r_sw + r_bw) + K / ((1 + K) * r_sw + K / g_bw),
          1 / (r_sw + r_bw))
}

# Calculate transpiration per area following the LI6800 manual
li6800_E = function(g_tw, w_a, w_i) {
  1000 * g_tw * (w_i - w_a) / (1000 - (w_i + w_a) / 2)
}

# Calculate 'true' reference [H2O] following the LI6800 manual
li6800_w0 = function(E, flow, s, w_a) {
  w_a - (s * E * 100) * (1000 - w_a) / (1000 * flow) 
}

# Calculate 'true' reference [CO2] following the LI6800 manual
li6800_c0 = function(A, c_a, flow, s, w_0) {
  1000 * (s * A * 100) / (1000 * flow) + c_a * ((1000 - w_0) / (1000 - w_a))
}
