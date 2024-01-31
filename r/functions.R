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

# Function to generate correlation matrix for simulating correlated error
make_autocorr_matrix = function(elapsed, b_autocorr) {
  # `elapsed` should be vector of time points in a series
  assert_numeric(elapsed, lower = 0, any.missing = FALSE, finite = TRUE)
  # `b_autocorr` should be a number indicating the autocorrelation decay rate
  assert_number(b_autocorr, lower = 0, finite = TRUE)
  
  m1 = outer(elapsed, elapsed, "-")
  R = diag(1, length(elapsed))
  
  ltri = exp(-b_autocorr * m1[which(lower.tri(m1))])
  utri = exp(-b_autocorr * t(m1)[which(upper.tri(m1))])
  
  R[which(lower.tri(R))] = ltri
  R[which(upper.tri(R))] = utri
 
  assert_true(isSymmetric(R))
  
  R
  
}
