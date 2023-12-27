data {
  int<lower=0> n_pts;
  int<lower=0> n_rep;
  int<lower=0> n_leaf_type;
  array[n_pts,n_rep,n_leaf_type] real log_gsw;
  array[n_pts,n_rep,n_leaf_type] real log_A;
}
parameters {
  real mu_intercept;
  real slope;
  real<lower=0> sigma_error_intercept;
  real<lower=0> sigma_error_resid;
  array[n_rep,n_leaf_type] real intercept;
}
model {
  
  // placeholders
  real y;
  real yhat;
  
  // priors
  mu_intercept ~ normal(0, 10);
  slope ~ normal(0, 10);
  sigma_error_intercept ~ cauchy(1,1);
  sigma_error_resid ~ cauchy(1,1);
  for (k in 1:n_leaf_type) {
    intercept[,k] ~ normal(0, sigma_error_intercept);
  }
  
  // model
  for (i in 1:n_pts) {
    for (j in 1:n_rep) {
      for (k in 1:n_leaf_type) {
        y = log_A[i,j,k];
        yhat = mu_intercept + intercept[j,k] + slope * log_gsw[i,j,k];
        y ~ normal(yhat, sigma_error_resid);
      }
    }
  }
  
}

