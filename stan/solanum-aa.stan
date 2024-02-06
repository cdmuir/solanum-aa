functions {
  // Calculate saturating vapor pressure following the LI6800 manual
  real li6800_svp(real T_degreeC, real P_kPa) {
    return 1000 * 0.61365 * exp(17.502 * T_degreeC / (240.97 + T_degreeC)) / P_kPa;
  }
  
  // Calculate total conductance to CO2 following the LI6800 manual
  real li6800_gtc_amphi(real g_bw, real g_sw, real K) {
    return 1 / ((1 + K) * (1.6 / g_sw) + (1.37 / g_bw)) + K / ((1 + K) * (1.6 / g_sw) + K * 1.37 / g_bw);
  }
  real li6800_gtc_hypo(real g_bw, real g_sw) {
    return 1 / (1.6 / g_sw + 1.37 / g_bw);
  }

  // Calculate total conductance to water vapor following the LI6800 manual
  real li6800_gtw_amphi(real g_bw, real g_sw, real K) {
    real r_bw = 1 / g_bw;
    real r_sw = 1 / g_sw;
    return 1 / ((1 + K) * r_sw + r_bw) + K / ((1 + K) * r_sw + K / g_bw);
  }
  real li6800_gtw_hypo(real g_bw, real g_sw) {
    real r_bw = 1 / g_bw;
    real r_sw = 1 / g_sw;
    return 1 / (r_sw + r_bw);
  }

  // Calculate transpiration per area following the LI6800 manual
  real li6800_E(real g_tw, real w_a, real w_i) {
    return 1000 * g_tw * (w_i - w_a) / (1000 - (w_i + w_a) / 2);
  }

  // Calculate 'true' reference [H2O] following the LI6800 manual
  real li6800_w0(real E, real flow, real s, real w_a) {
    return w_a - (s * E * 100) * (1000 - w_a) / (1000 * flow);
  }

  // Calculate 'true' reference [CO2] following the LI6800 manual
  real li6800_c0(real A, real c_a, real flow, real s, real w_a, real w_0) {
    return 1000 * (s * A * 100) / (1000 * flow) + c_a * ((1000 - w_0) / (1000 - w_a));
  }
}
data {
  int<lower=0> n_pts;
  int<lower=0> n_id;
  int<lower=0> n_leaf_type;  
  int<lower=0> n_light_treatment;  
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real elapsed;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real flow;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real g_bw;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real K;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real P;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real RH;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real s;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real T_air;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real T_leaf;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real CO2_r;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real CO2_s;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real H2O_r;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real H2O_s;
}
parameters {
  // hyperparameters
  real<lower=0> b_autocorr_c; 
  real<lower=0> b_autocorr_w; 
  real<lower=0> sigma_c; 
  real<lower=0> sigma_w; 

  real mu_intercept;
  real mu_intercept_low_light;
  real<lower=0> sigma_intercept_id;
  real<lower=0> sigma_intercept_low_light_id;
  real<lower=0> sigma_intercept_error;

  real mu_slope;
  real mu_slope_low_light;
  real<lower=0> sigma_slope_id;
  real<lower=0> sigma_slope_low_light_id;
  
  // parameters
  array[n_id] real b_intercept_id;
  array[n_id] real b_intercept_low_light_id;
  array[n_id,n_leaf_type,n_light_treatment] real b_intercept_error;

  array[n_id] real b_slope_id;
  array[n_id] real b_slope_low_light_id;

  array[n_pts,n_id,n_leaf_type,n_light_treatment] real log_gsw;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real c_a;
  
}
transformed parameters {
  
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real A;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real g_sw  = exp(log_gsw);
  
  // calculated quantities
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real w_i;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real g_tc;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real g_tw;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real E;

  // real [CO2] and [H2O]
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real c_0;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real w_0;
  array[n_pts,n_id,n_leaf_type,n_light_treatment] real w_a;
  
  // error covariance matrices
  array[n_pts,n_pts,n_id,n_leaf_type,n_light_treatment] real R_c;
  array[n_pts,n_pts,n_id,n_leaf_type,n_light_treatment] real R_w;
  
  // calculations
  for (l in 1:n_light_treatment) {
    for (k in 1:n_leaf_type) {
      for (j in 1:n_id) {
        for (i2 in 1:n_pts) {
        
          A[i2,j,k,l] = mu_intercept + 
            (mu_intercept_low_light + b_intercept_low_light_id[j]) * (l - 1) +
            b_intercept_id[j] + b_intercept_error[j,k,l] + 
            (mu_slope + (mu_slope_low_light + b_slope_low_light_id[j]) * (l - 1) + 
              b_slope_id[j]) * log_gsw[i2,j,k,l];

          w_i[i2,j,k,l]  = li6800_svp(T_leaf[i2,j,k,l], P[i2,j,k,l]);
          w_a[i2,j,k,l]  = RH[i2,j,k,l] * li6800_svp(T_air[i2,j,k,l], P[i2,j,k,l]);
        
          // amphi leaves
          if (k == 1) {
            g_tc[i2,j,k,l] = li6800_gtc_amphi(g_bw[i2,j,k,l], g_sw[i2,j,k,l], K[i2,j,k,l]);
            g_tw[i2,j,k,l] = li6800_gtw_amphi(g_bw[i2,j,k,l], g_sw[i2,j,k,l], K[i2,j,k,l]);
          } 
          // pseudohypo leaves
          if (k == 2) {
            g_tc[i2,j,k,l] = li6800_gtc_hypo(g_bw[i2,j,k,l], g_sw[i2,j,k,l]);
            g_tw[i2,j,k,l] = li6800_gtw_hypo(g_bw[i2,j,k,l], g_sw[i2,j,k,l]);
          } 

          E[i2,j,k,l]    = li6800_E(g_tw[i2,j,k,l], w_a[i2,j,k,l], w_i[i2,j,k,l]);
          w_0[i2,j,k,l]  = li6800_w0(E[i2,j,k,l], flow[i2,j,k,l], s[i2,j,k,l], w_a[i2,j,k,l]);
          c_0[i2,j,k,l]  = li6800_c0(A[i2,j,k,l], c_a[i2,j,k,l], flow[i2,j,k,l], s[i2,j,k,l], w_a[i2,j,k,l], w_0[i2,j,k,l]);

          for (i1 in 1:n_pts) {
            R_c[i1, i2, j, k, l] = exp(-b_autocorr_c * abs(elapsed[i2, j, k, l] - elapsed[i1, j, k, l]));
            R_w[i1, i2, j, k, l] = exp(-b_autocorr_w * abs(elapsed[i2, j, k, l] - elapsed[i1, j, k, l]));
          }
        }
      }
    }
  }
}
model {
  
  // placeholders
  // nothing right now
  
  // priors on hyperparameters
  b_autocorr_c ~ normal(0, 1); 
  b_autocorr_w ~ normal(0, 1); 
  // sigma_c ~ gamma(2, 0); // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  sigma_c ~ normal(0, 1); 
  sigma_w ~ normal(0, 1); 
  
  mu_intercept ~ normal(30, 10);
  mu_intercept_low_light ~ normal(0, 20);
  sigma_intercept_id ~ normal(0, 10);
  sigma_intercept_low_light_id ~ normal(0, 10);
  sigma_intercept_error ~ normal(0, 10);

  mu_slope ~ normal(10, 5);
  mu_slope_low_light ~ normal(0, 5);
  sigma_slope_id ~ normal(0, 10);
  sigma_slope_low_light_id ~ normal(0, 10);

  // priors on parameters
  b_intercept_id ~ normal(0, sigma_intercept_id);
  b_intercept_low_light_id ~ normal(0, sigma_intercept_low_light_id);
  b_slope_id ~ normal(0, sigma_slope_id);
  b_slope_low_light_id ~ normal(0, sigma_slope_low_light_id);

  for (l in 1:n_light_treatment) {
    for (k in 1:n_leaf_type) {
      b_intercept_error[,k,l] ~ normal(0, sigma_intercept_error);
      for (j in 1:n_id) {
        for (i in 1:n_pts) {
          c_a[i,j,k,l] ~ normal(415, 1);
          log_gsw[i,j,k,l] ~ normal(-1, 1); // is this the right choice?
        }
      }
    }
  }
  
  // likelihood
  for (l in 1:n_light_treatment) {
    for (k in 1:n_leaf_type) {
      for (j in 1:n_id) {
      
        // placeholders
        row_vector[n_pts] x;
        row_vector[n_pts] y;
        matrix[n_pts, n_pts] R;
        vector[n_pts] s_vec;
      
        // CO2_r
        for (i2 in 1:n_pts) {
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; c_0[i,j,k] = ", c_0[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; A[i,j,k] = ", A[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; c_a[i,j,k] = ", c_a[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; flow[i,j,k] = ", flow[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; s[i,j,k] = ", s[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; w_a[i,j,k] = ", w_a[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; w_0[i,j,k] = ", w_0[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; E[i,j,k] = ", E[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; w_i[i,j,k] = ", w_i[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; g_tw[i,j,k] = ", g_tw[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; g_sw[i,j,k] = ", g_sw[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; g_bw[i,j,k] = ", g_bw[i2,j,k]);
          // print("i = ", i2, "; j = ", j, "; k = ", k, "; K[i,j,k] = ", K[i2,j,k]);
          x[i2] = c_0[i2,j,k,l];
          y[i2] = CO2_r[i2,j,k,l];
          s_vec[i2] = sigma_c;
          for (i1 in 1:n_pts) {
            R[i1, i2] = R_c[i1, i2, j, k, l];
          }
        }
        y ~ multi_normal(x, quad_form_diag(R, s_vec));

        // CO2_s
        for (i2 in 1:n_pts) {
          x[i2] = c_a[i2,j,k,l];
          y[i2] = CO2_s[i2,j,k,l];
          s_vec[i2] = sigma_c;
          for (i1 in 1:n_pts) {
            R[i1, i2] = R_c[i1, i2, j, k, l];
          }
        }
        y ~ multi_normal(x, quad_form_diag(R, s_vec));

        // H2O_r
        for (i2 in 1:n_pts) {
          x[i2] = w_0[i2,j,k,l];
          y[i2] = H2O_r[i2,j,k,l];
          s_vec[i2] = sigma_w;
          for (i1 in 1:n_pts) {
            R[i1, i2] = R_w[i1, i2, j, k, l];
          }
        }
        y ~ multi_normal(x, quad_form_diag(R, s_vec));

        // H2O_s
        for (i2 in 1:n_pts) {
          x[i2] = w_a[i2,j,k,l];
          y[i2] = H2O_s[i2,j,k,l];
          s_vec[i2] = sigma_w;
          for (i1 in 1:n_pts) {
            R[i1, i2] = R_w[i1, i2, j, k, l];
          }
        }
        y ~ multi_normal(x, quad_form_diag(R, s_vec));
      
      }
    }
  }
}

