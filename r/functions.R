# 0. Import data ----

# 1. Simulate synthetic data ----
# n.b. Many specialized functions related to simulating and fitting LI6800 data 
# are in r/licor-functions.R

# Function to calculate the decay, per s, in autocorrelation between data points
calculate_corr_decay = function(rho, t) {
  assert_numeric(rho, lower = -1, upper = 1, any.missing = FALSE, min.len = 1L)
  assert_numeric(t, lower = 0, any.missing = FALSE, min.len = 1L)
  assert_true(length(t) %in% c(1L, length(rho)))
  # rho = exp(-b*t)
  # log(rho) = - b * t
  # -log(rho) / t = b
  -log(rho) / t
}
