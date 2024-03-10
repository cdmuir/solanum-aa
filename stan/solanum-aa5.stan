// This model fits a quadratic function of log_gsw to log_A for every curve, then
// integrates overlapping region to estimate AA. Then we estimate parameters 
// describing effects of variables on AA. I still need to check with log-log fit
// is reasonable.
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
  
  int<lower=0> n_id;
  int<lower=0> n_light_intensity;  
  int<lower=0> n_lightintensity_x_id;
  int<lower=0> n_light_treatment;  
  
  // vector of integer lengths per curve
  array[n_curve] int<lower=0> n_pts; 

  // data indexed by row
  array[n] int<lower=1,upper=n_curve> curve;
  vector[n] A;
  vector[n] scaled_log_gsw;

  // index of amphi and pseudohypo curves for each light_intensity x id combination
  array[n_lightintensity_x_id] int<lower=0,upper=n_curve> amphi;
  array[n_lightintensity_x_id] int<lower=0,upper=n_curve> pseudohypo;

  // variables indexed by lightintensity_x_id groups
  array[n_lightintensity_x_id] int<lower=1,upper=n_id> id;
  array[n_lightintensity_x_id] int<lower=1,upper=n_light_intensity> light_intensity;
  array[n_lightintensity_x_id] int<lower=1,upper=n_light_treatment> light_treatment;

  // min and max scaled_log_gsw by curve
  array[n_curve] real min_scaled_log_gsw;
  array[n_curve] real max_scaled_log_gsw;
  
}
transformed data {
  vector[n] log_A;
  log_A = log(A);
}
parameters {
  
  // quadratic regression coefficients
  vector[n_curve] b0;
  vector[n_curve] b1;
  vector[n_curve] b2;

  real log_sigma_resid;
  real<lower=-1,upper=1> rho_resid;
  
  // regression on aa
  real b0_aa;
  real b_aa_light_intensity_2000;
  real b_aa_light_treatment_high;
  vector[n_id] b_aa_id;
  real log_sigma_aa_id;
  
  // regression on sigma_aa
  real b0_log_sigma_aa;
  real b_log_sigma_aa_light_intensity_2000;
  real b_log_sigma_aa_light_treatment_high;
  
}
transformed parameters {
  real sigma_resid;
  real sigma_aa_id;
  
  sigma_resid = exp(log_sigma_resid);
  sigma_aa_id = exp(log_sigma_aa_id);
}
model {
  
  profile("priors") {
  // priors ----
    
    // quadratic regression coefficients
    b0 ~ normal(0, 10);
    b1 ~ normal(0, 10);
    b2 ~ normal(0, 10);
  
    log_sigma_resid ~ normal(0, 1); 
    rho_resid ~ normal(0, 1);
    
    // regression on aa
    b0_aa ~ normal(0, 1);
    b_aa_light_intensity_2000 ~ normal(0, 1);
    b_aa_light_treatment_high ~ normal(0, 1);
    b_aa_id ~ normal(0, sigma_aa_id);
    log_sigma_aa_id ~ normal(0, 1);
    
    // regression on sigma_aa
    b0_log_sigma_aa ~ normal(0, 1);
    b_log_sigma_aa_light_intensity_2000 ~ normal(0, 1);
    
  }
  
      profile("Estimate AA") {
  // Estimate AA
  
  for (i in 1:n_lightintensity_x_id) {
    
    real mu1;
    real sigma;

    // regression on aa
    mu1 = b0_aa + 
      b_aa_light_intensity_2000 * (light_intensity[i] == 2) +
      b_aa_light_treatment_high * (light_treatment[i] == 2) +
      b_aa_id[id[i]];
    
    // regression on sigma_aa
    sigma = exp(b0_log_sigma_aa +
      b_log_sigma_aa_light_intensity_2000 * (light_intensity[i] == 2) +
      b_log_sigma_aa_light_treatment_high * (light_treatment[i] == 2));

    real aa_i;
    int amphi_curve; 
    int pseudohypo_curve;
    
    real a; // min of amphi curve
    real b; // max of hypo
    array[6] real theta; // parameters
    
    amphi_curve = amphi[i];
    pseudohypo_curve = pseudohypo[i];
    
    a = min_scaled_log_gsw[amphi_curve];
    b = max_scaled_log_gsw[pseudohypo_curve];
    
    theta[1] = b0[amphi_curve];      // b0_amphi;
    theta[2] = b1[amphi_curve];      // b1_amphi;
    theta[3] = b2[amphi_curve];      // b2_amphi;
    theta[4] = b0[pseudohypo_curve]; // b0_hypo;
    theta[5] = b1[pseudohypo_curve]; // b1_hypo;
    theta[6] = b2[pseudohypo_curve]; // b2_hypo;
    
    aa_i = aa_int(b, theta) - aa_int(a, theta);

    target += normal_lpdf(aa_i | mu1, sigma);

  }
  }

  profile("likelihood") {
  // likelihood ----
  vector[n] resid;
  vector[n] mu2;
  
    for (i in 1:n) {
      
      mu2[i] = b0[curve[i]] + 
        b1[curve[i]] * scaled_log_gsw[i] +
        b2[curve[i]] * scaled_log_gsw[i] ^ 2;
  
    }
    
    resid[1] = log_A[1] - mu2[1];
    for (i in 2:n) {
      resid[i] = log_A[i] - mu2[i];
      mu2[i] += rho_resid * resid[i - 1] * (curve[i] == curve[i - 1]);
    }
    
  target += normal_lpdf(log_A | mu2, sigma_resid);
  
  }
}
