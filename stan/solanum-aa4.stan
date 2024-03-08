functions {
  // Estimate AA by integrating over amphi- and pseudohypo-curves
  real aa_int(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    return log((theta[1] + theta[2] * x + theta[3] * x ^ 2) / (theta[4] + theta[5] * x + theta[6] * x ^ 2));
  }
}
data {
  
  // total number of rows
  int<lower=0> n;
  // total number of curves
  int<lower=0> n_curve; 
  
  int<lower=0> n_id;
  int<lower=0> n_leaf_type;  
  int<lower=0> n_light_intensity;  
  int<lower=0> n_lightintensity_x_id;
  int<lower=0> n_light_treatment;  
  
  // vector of integer lengths per curve
  array[n_curve] int<lower=0> n_pts; 

  array[n] int<lower=1,upper=n_curve> curve;
  
  array[n] int<lower=1,upper=n_id> id;
  array[n] int<lower=1,upper=n_leaf_type> leaf_type;
  array[n] int<lower=1,upper=n_lightintensity_x_id> lightintensity_x_id;
  
  vector[n] A;
  vector[n] scaled_log_gsw;

  // index of amphi and pseudohypo curves for each light_intensity x id combination
  array[n_lightintensity_x_id] int<lower=0,upper=n_curve> amphi;
  array[n_lightintensity_x_id] int<lower=0,upper=n_curve> pseudohypo;

  // variables indexed by lightintensity_x_id groups
  array[n_lightintensity_x_id] int<lower=1,upper=n_light_intensity> light_intensity;
  array[n_lightintensity_x_id] int<lower=1,upper=n_light_treatment> light_treatment;

  // min and max scaled_log_gsw by curve
  array[n_curve] real min_scaled_log_gsw;
  array[n_curve] real max_scaled_log_gsw;
  
  array[1] real d1;
  array[1] int d2;

}
parameters {
  
  // quadratic regression coefficients
  vector[n_curve] b0;
  vector[n_curve] b1;
  vector[n_curve] b2;

  real<lower=0> sigma_resid;

  // regression on aa
  real b0_aa;
  real b_aa_light_treatment_high;
  
  // regression on sigma_aa
  real b0_log_sigma_aa;
  real b_log_sigma_aa_light_intensity_2000;
  
}
model {
  
  profile("priors") {
  // priors ----
    
    // quadratic regression coefficients
    b0 ~ normal(0, 10);
    b1 ~ normal(0, 10);
    b2 ~ normal(0, 10);
  
    sigma_resid ~ normal(0, 1); 
    
    // regression on aa
    b0_aa ~ normal(0, 1);
    b_aa_light_treatment_high ~ normal(0, 1);

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
      b_aa_light_treatment_high * (light_treatment[i] == 2);
    
    // regression on sigma_aa
    sigma = b0_log_sigma_aa +
      b_log_sigma_aa_light_intensity_2000 * (light_intensity[i] == 2);

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
    
    aa_i = integrate_1d(aa_int, a, b, theta, d1, d2);

    target += normal_lpdf(aa_i | mu1, sigma);

  }
  }

  profile("likelihood") {
  // likelihood ----
    for (i in 1:n) {
      
      real mu2;
      mu2 = b0[curve[i]] + 
        b1[curve[i]] * scaled_log_gsw[i] +
        b2[curve[i]] * scaled_log_gsw[i] ^ 2;
      target += normal_lpdf(A[i] | mu2, sigma_resid);
  
    }
  }
  
}
