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

  // Function to calculate the mean of an 4-dimensional array
  real array_mean_4d(array[,,,] real x) {
    int n_dim1 = dims(x)[1];
    int n_dim2 = dims(x)[2];
    int n_dim3 = dims(x)[3];
    int n_dim4 = dims(x)[4];

    real mu = 0.0;

    // Calculate the mean
    for (i in 1:n_dim1) {
      for (j in 1:n_dim2) {
        for (k in 1:n_dim3) {
          for (l in 1:n_dim4) {
            mu += x[i, j, k, l];
          }
        }
      }
    }
    mu /= (n_dim1 * n_dim2 * n_dim3 * n_dim4);

    return mu;
    
  }

  // Function to calculate the sd of an 4-dimensional array
  real array_sd_4d(array[,,,] real x) {
    int n_dim1 = dims(x)[1];
    int n_dim2 = dims(x)[2];
    int n_dim3 = dims(x)[3];
    int n_dim4 = dims(x)[4];

    real sum_sq_diff = 0.0;
    real mu = array_mean_4d(x);

    // Calculate the sum of squared differences
    for (i in 1:n_dim1) {
      for (j in 1:n_dim2) {
        for (k in 1:n_dim3) {
          for (l in 1:n_dim4) {
            sum_sq_diff += (x[i, j, k, l] - mu) ^ 2;
          }
        }
      }
    }

    // Calculate the standard deviation
    real std_dev = sqrt(sum_sq_diff / (n_dim1 * n_dim2 * n_dim3 * n_dim4));

    return std_dev;
    
  }

}
data {
  int<lower=0> n_pts;
  int<lower=0> n_id;
  int<lower=0> n_leaf_type;  
  int<lower=0> n_light_treatment;  
  int<lower=0> n_comp;  
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
  // variable scaling
  vector<lower=0>[n_comp] sd_A;
  vector<lower=0>[n_comp] mean_A;
  
  // hyperparameters
  vector<lower=0>[n_comp] b_autocorr_c; 
  vector<lower=0>[n_comp] b_autocorr_w; 
  vector<lower=0>[n_comp] sigma_c; 
  vector<lower=0>[n_comp] sigma_w; 

  vector[n_comp] mu_intercept;
  vector[n_comp] mu_intercept_low_light;
  vector<lower=0>[n_comp] sigma_intercept_id;
  vector<lower=0>[n_comp] sigma_intercept_low_light_id;
  vector<lower=0>[n_comp] sigma_intercept_error;

  vector[n_comp] mu_slope;
  vector[n_comp] mu_slope_low_light;
  vector<lower=0>[n_comp] sigma_slope_id;
  vector<lower=0>[n_comp] sigma_slope_low_light_id;
  
  // parameters
  array[n_id,n_comp] real b_intercept_id;
  array[n_id,n_comp] real b_intercept_low_light_id;
  array[n_id,n_leaf_type,n_light_treatment,n_comp] real b_intercept_error;

  array[n_id,n_comp] real b_slope_id;
  array[n_id,n_comp] real b_slope_low_light_id;

  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real log_gsw;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real c_a;
 
  // Component 2 parameters
  real mu_b2;
  real mu_b2_low_light;
  vector[n_id] b_b2_low_light_id;
  vector[n_id] b_b2_id;
  real<lower=0> sigma_b2_id;
  real<lower=0> sigma_b2_low_light_id;

  // component weights
  array[n_id,n_leaf_type,n_light_treatment,n_comp] simplex[n_comp] w;
  
}
transformed parameters {
  
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real A;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real g_sw = exp(log_gsw);
  real sd_log_gsw[n_comp];
  real mean_log_gsw[n_comp];
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real scaled_log_gsw;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real scaled_A;
  
  // calculated quantities
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real w_i;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real g_tc;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real g_tw;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real E;

  // real [CO2] and [H2O]
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real c_0;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real w_0;
  array[n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real w_a;
  
  // error covariance matrices
  array[n_pts,n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real R_c;
  array[n_pts,n_pts,n_id,n_leaf_type,n_light_treatment,n_comp] real R_w;
  
  // calculations
  for (z in 1:n_comp) {
    sd_log_gsw[z] = array_sd_4d(log_gsw[,,,,z]);
    mean_log_gsw[z] = array_mean_4d(log_gsw[,,,,z]);
    for (l in 1:n_light_treatment) {
      for (k in 1:n_leaf_type) {
        for (j in 1:n_id) {
          for (i2 in 1:n_pts) {
        
            scaled_log_gsw[i2,j,k,l,z] = (log_gsw[i2,j,k,l,z] - mean_log_gsw[z]) / sd_log_gsw[z];
        
            // Component 1: linear A-log(gsw)
            if (z == 1) {
              scaled_A[i2,j,k,l,z] = mu_intercept[z] + 
                (mu_intercept_low_light[z] + b_intercept_low_light_id[j,z]) * (l - 1) +
                b_intercept_id[j,z] + b_intercept_error[j,k,l,z] + 
                (mu_slope[z] + (mu_slope_low_light[z] + b_slope_low_light_id[j,z]) * (l - 1) + 
                  b_slope_id[j,z]) * scaled_log_gsw[i2,j,k,l,z];
            }

            // Component 2: quadratic A-log(gsw)
            if (z == 2) {
              scaled_A[i2,j,k,l,z] = mu_intercept[z] + 
                (mu_intercept_low_light[z] + b_intercept_low_light_id[j,z]) * (l - 1) +
                b_intercept_id[j,z] + b_intercept_error[j,k,l,z] + 
                (mu_slope[z] + (mu_slope_low_light[z] + b_slope_low_light_id[j,z]) * (l - 1) + 
                  b_slope_id[j,z]) * scaled_log_gsw[i2,j,k,l,z] +
                (mu_b2 + (mu_b2_low_light + b_b2_low_light_id[j]) * (l - 1) + 
                  b_b2_id[j]) * scaled_log_gsw[i2,j,k,l,z] ^ 2;
            }

            A[i2,j,k,l,z] = scaled_A[i2,j,k,l,z] * sd_A[z] + mean_A[z];
            
            w_i[i2,j,k,l,z] = li6800_svp(T_leaf[i2,j,k,l], P[i2,j,k,l]);
            w_a[i2,j,k,l,z] = RH[i2,j,k,l] * li6800_svp(T_air[i2,j,k,l], P[i2,j,k,l]);
        
            // amphi leaves
            if (k == 1) {
              g_tc[i2,j,k,l,z] = li6800_gtc_amphi(g_bw[i2,j,k,l], g_sw[i2,j,k,l,z], K[i2,j,k,l]);
              g_tw[i2,j,k,l,z] = li6800_gtw_amphi(g_bw[i2,j,k,l], g_sw[i2,j,k,l,z], K[i2,j,k,l]);
            } 
            // pseudohypo leaves
            if (k == 2) {
              g_tc[i2,j,k,l,z] = li6800_gtc_hypo(g_bw[i2,j,k,l], g_sw[i2,j,k,l,z]);
              g_tw[i2,j,k,l,z] = li6800_gtw_hypo(g_bw[i2,j,k,l], g_sw[i2,j,k,l,z]);
            } 

            E[i2,j,k,l,z]    = li6800_E(g_tw[i2,j,k,l,z], w_a[i2,j,k,l,z], w_i[i2,j,k,l,z]);
            w_0[i2,j,k,l,z]  = li6800_w0(E[i2,j,k,l,z], flow[i2,j,k,l], s[i2,j,k,l], w_a[i2,j,k,l,z]);
            c_0[i2,j,k,l,z]  = li6800_c0(A[i2,j,k,l,z], c_a[i2,j,k,l,z], flow[i2,j,k,l], s[i2,j,k,l], w_a[i2,j,k,l,z], w_0[i2,j,k,l,z]);

            for (i1 in 1:n_pts) {
              R_c[i1, i2, j, k, l, z] = exp(-b_autocorr_c[z] * abs(elapsed[i2, j, k, l] - elapsed[i1, j, k, l]));
              R_w[i1, i2, j, k, l, z] = exp(-b_autocorr_w[z] * abs(elapsed[i2, j, k, l] - elapsed[i1, j, k, l]));
            }
          }
        }
      }
    }
  }
}
model {
  
  // placeholders
  // nothing right now
  
  // priors on variable scaling
  sd_A ~ normal(0, 10);
  mean_A ~ normal(15, 10);
  
  // priors on hyperparameters
  b_autocorr_c ~ normal(0, 1); 
  b_autocorr_w ~ normal(0, 1); 
  // sigma_c ~ gamma(2, 0); // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  sigma_c ~ normal(0, 1); 
  sigma_w ~ normal(0, 1); 
  
  // priors for unstandardized variables
  // mu_intercept ~ normal(30, 10);
  // mu_intercept_low_light ~ normal(0, 20);
  mu_intercept ~ normal(0, 10);
  mu_intercept_low_light ~ normal(0, 10);
  sigma_intercept_id ~ normal(0, 10);
  sigma_intercept_low_light_id ~ normal(0, 10);
  sigma_intercept_error ~ normal(0, 10);

  // mu_slope ~ normal(10, 5); // priors for unstandardized variables
  mu_slope ~ normal(0, 10); // priors for unstandardized variables
  mu_slope_low_light ~ normal(0, 5);
  sigma_slope_id ~ normal(0, 10);
  sigma_slope_low_light_id ~ normal(0, 10);

  // priors on parameters
  for (z in 1:n_comp) {
    b_intercept_id[z] ~ normal(0, sigma_intercept_id[z]);
    b_intercept_low_light_id[z] ~ normal(0, sigma_intercept_low_light_id[z]);
    b_slope_id[z] ~ normal(0, sigma_slope_id[z]);
    b_slope_low_light_id[z] ~ normal(0, sigma_slope_low_light_id[z]);
  }
  
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

