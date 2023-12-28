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
  real<lower=-1,upper=1> rho_error_resid;
  array[n_rep,n_leaf_type] real error_resid1;
  real<lower=0> sigma_error_intercept;
  real<lower=0> sigma_error_resid;
  array[n_rep,n_leaf_type] real intercept;
}
model {
  
  // placeholders
  real y;
  array[n_pts,n_rep,n_leaf_type] real mu = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  // brms syntax (change to my own once I understand)
  // array storing lagged residuals
  array[n_pts,n_rep,n_leaf_type] real Err = rep_array(0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real err;  // actual residuals
  
  // priors
  mu_intercept ~ normal(0, 10);
  slope ~ normal(0, 10);
  rho_error_resid ~ normal(0, 1);
  sigma_error_intercept ~ cauchy(1,1);
  sigma_error_resid ~ cauchy(1,1);
  for (k in 1:n_leaf_type) {
    intercept[,k] ~ normal(0, sigma_error_intercept);
  }
  
  // model
  for (k in 1:n_leaf_type) {
    for (j in 1:n_rep) {
      
      mu[1,j,k] += mu_intercept + intercept[j,k] + slope * log_gsw[1,j,k];
      err[1,j,k] = log_A[1,j,k] - mu[1,j,k];
      y = log_A[1,j,k];
      y ~ normal(mu[1,j,k], sigma_error_resid);
        
      for (i in 2:n_pts) {
        mu[i,j,k] += mu_intercept + intercept[j,k] + slope * log_gsw[i,j,k];
        err[i,j,k] = log_A[i,j,k] - mu[i,j,k];
        Err[i,j,k] = err[i-1,j,k];
        mu[i,j,k] += Err[i,j,k] * rho_error_resid;
        y = log_A[i,j,k];
        y ~ normal(mu[i,j,k], sigma_error_resid);
      }
      
    }
  }
  
}

