functions {
  // indefinite integral of log(A_amphi) - log(A_hypo) based on parameters theta
  real aa_int(real x, array[] real theta) {
    
    return x ^ 3 * (theta[3] / 3 - theta[6] / 3) + 
      x ^ 2 * (theta[2] / 2 - theta[5] / 2) + 
      x * (theta[1] - theta[4]);
    
  }
  
}
data {
  
  // total number of rows
  int<lower=0> n;
  // total number of curves
  int<lower=0> n_curve; 
  
  int<lower=0> n_acc;
  int<lower=0> n_acc_id;
  int<lower=0> n_light_intensity;  
  int<lower=0> n_lightintensity_x_acc_id;
  int<lower=0> n_light_treatment;  
  
  // vector of integer lengths per curve
  array[n_curve] int<lower=0> n_pts; 

  // data indexed by row
  array[n] int<lower=1,upper=n_curve> curve;
  vector[n] A;
  vector[n] scaled_log_gsw;

  // index of amphi and pseudohypo curves for each light_intensity x id combination
  array[n_lightintensity_x_acc_id] int<lower=0,upper=n_curve> amphi;
  array[n_lightintensity_x_acc_id] int<lower=0,upper=n_curve> pseudohypo;

  // variables indexed by lightintensity_x_id groups
  array[n_lightintensity_x_acc_id] int<lower=1,upper=n_acc> acc;
  array[n_lightintensity_x_acc_id] int<lower=1,upper=n_acc_id> acc_id;
  array[n_lightintensity_x_acc_id] int<lower=0,upper=1> amphi_first;
  array[n_lightintensity_x_acc_id] int<lower=1,upper=n_light_intensity> light_intensity;
  array[n_lightintensity_x_acc_id] int<lower=1,upper=n_light_treatment> light_treatment;

  // min and max scaled_log_gsw by curve
  array[n_curve] real min_scaled_log_gsw;
  array[n_curve] real max_scaled_log_gsw;
  array[n_curve] real S;
  
}
transformed data {
  vector[n] log_A;
  log_A = log(A);
}
parameters {
  array[n_curve] vector[3] B_curve;
  vector[3] Mu_curve;
  vector[3] log_sigma_curve;
  corr_matrix[3] R_curve;
  real<lower=-1, upper=1> rho_resid;
  real b0_log_sigma_resid;
  real b_log_sigma_resid_S;
}
transformed parameters {
  vector[3] sigma_curve;
  sigma_curve = exp(log_sigma_curve);
  matrix[3,3] Sigma_curve;
  Sigma_curve = quad_form_diag(R_curve, sigma_curve);
}
model {

  // priors
  B_curve ~ multi_normal(Mu_curve, Sigma_curve);
  Mu_curve ~ normal(0, 10);
  log_sigma_curve ~ normal(0,1);
  R_curve ~ lkj_corr(2);
  b0_log_sigma_resid ~ normal(-3,5);
  b_log_sigma_resid_S ~ normal(0,1);

  // Estimate RH curve parameters
  vector[n] mu2;
  vector[n] resid;
  vector[n] sigma2;
  
  for (i in 1:n) {
      
    mu2[i] = Mu_curve[1] + B_curve[curve[i],1] + 
      (Mu_curve[2] + B_curve[curve[i],2]) * scaled_log_gsw[i] +
      (Mu_curve[3] + B_curve[curve[i],3]) * scaled_log_gsw[i] ^ 2;
    
    sigma2[i] = exp(b0_log_sigma_resid + (6 - S[curve[i]]) * b_log_sigma_resid_S);
  
  }
  
  resid[1] = log_A[1] - mu2[1];
  for (i in 2:n) {
    resid[i] = log_A[i] - mu2[i];
    mu2[i] += rho_resid * resid[i - 1] * (curve[i] == curve[i - 1]);
  }
  
  target += normal_lpdf(log_A | mu2, sigma2);
}
