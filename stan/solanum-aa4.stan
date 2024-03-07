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
  
  // vector of integer lengths per curve
  array[n_curve] int<lower=0> n_pts; 

  array[n] int<lower=1,upper=n_curve> curve;
  
  array[n] int<lower=1,upper=n_id> id;
  array[n] int<lower=1,upper=n_leaf_type> leaf_type;
  array[n] int<lower=1,upper=n_light_intensity> light_intensity;
  array[n] int<lower=1,upper=n_lightintensity_x_id> lightintensity_x_id;
  
  vector[n] A;
  vector[n] scaled_log_gsw;

  array[1] real d1;
  array[1] int d2;

}
parameters {
  
  // quadratic regression coefficients
  vector<lower=0>[n_curve] b0;
  vector<lower=0>[n_curve] b1;
  vector<lower=0>[n_curve] b2;

  real<lower=0> sigma;

}
model {
  
  profile("priors") {
  // priors ----
    b0 ~ normal(0, 10);
    b1 ~ normal(0, 10);
    b2 ~ normal(0, 10);
  
    sigma ~ normal(0, 1); 
  }
  
  profile("likelihood") {
  // likelihood ----
    for (i in 1:n) {
      
      real mu;
      mu = b0[curve[i]] + 
        b1[curve[i]] * scaled_log_gsw[i] +
        b2[curve[i]] * scaled_log_gsw[i] ^ 2;
      target += normal_lpdf(A[i] | mu, sigma);
  
    }
  }
  
}
generated quantities {
  
  # Calculate AA
  vector[n_lightintensity_x_id] aa;
  
  for (i in 1:n_lightintensity_x_id) {
    
    real a; // min of amphi curve
    real b; // max of hypo
    array[6] real theta; // parameters
    
    a = min(scaled_log_gsw[???]);
    b = max(scaled_log_gsw[???]);
    
    theta[1] = b0[???] // b0_amphi;
    theta[2] = b1[???] // b1_amphi;
    theta[3] = b2[???] // b2_amphi;
    theta[4] = b0[???] // b0_hypo;
    theta[5] = b1[???] // b1_hypo;
    theta[6] = b2[???] // b2_hypo;
    
    aa[i] = integrate_1d(aa_int, a, b, theta, d1, d2);

  }
}
