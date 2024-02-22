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
  
  // total number of rows
  int<lower=0> n;
  // total number of curves
  // this should be n_id * n_leaf_type * n_light_treatment if perfectly balanced
  int<lower=0> n_curve; 
  
  int<lower=0> n_id;
  int<lower=0> n_leaf_type;  
  int<lower=0> n_light_treatment;  
  
  int<lower=1> n_comp;  
  
  // vector of integer lengths per curve
  int<lower=0> n_pts[n_curve]; 
  
  int<lower=1,upper=n_id> id[n];
  int<lower=1,upper=n_leaf_type> leaf_type[n];
  int<lower=1,upper=n_light_treatment> light_treatment[n];
  
  vector[n] elapsed;
  vector[n] flow;
  vector[n] g_bw;
  vector[n] K;
  vector[n] P;
  vector[n] RH;
  vector[n] s;
  vector[n] T_air;
  vector[n] T_leaf;
  vector[n] CO2_r;
  vector[n] CO2_s;
  vector[n] H2O_r;
  vector[n] H2O_s;
  
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
  array[n_curve,n_comp] real b_intercept_error;

  array[n_id,n_comp] real b_slope_id;
  array[n_id,n_comp] real b_slope_low_light_id;

  array[n,n_comp] real log_gsw;
  array[n,n_comp] real c_a;
 
  // Component 2 parameters
  real mu_b2;
  real mu_b2_low_light;
  vector[n_id] b_b2_low_light_id;
  vector[n_id] b_b2_id;
  real<lower=0> sigma_b2_id;
  real<lower=0> sigma_b2_low_light_id;

  // component weights
  array[n_curve] simplex[n_comp] w;
  
}
model {
  
  // placeholders ----
  int ir;
  
  array[n,n_comp] real A;
  array[n,n_comp] real g_sw = exp(log_gsw);
  vector[n_comp] sd_log_gsw;
  vector[n_comp] mean_log_gsw;
  array[n,n_comp] real scaled_log_gsw;
  array[n,n_comp] real scaled_A;
  
  // calculated quantities
  array[n,n_comp] real w_i;
  array[n,n_comp] real g_tc;
  array[n,n_comp] real g_tw;
  array[n,n_comp] real E;

  // real [CO2] and [H2O]
  array[n,n_comp] real c_0;
  array[n,n_comp] real w_0;
  array[n,n_comp] real w_a;
  
  // calculations ----
  for (z in 1:n_comp) {
    sd_log_gsw[z] = sd(log_gsw[1:n,z]);
    mean_log_gsw[z] = mean(log_gsw[1:n,z]);
    for (k in 1:n_curve) {
      for (j in 1:n_pts[k]) {
            
        // row position
        ir = sum(n_pts[0:(k - 1)]) + j;
            
        scaled_log_gsw[ir,z] = (log_gsw[ir,z] - mean_log_gsw[z]) / sd_log_gsw[z];
        
        // Component 1: linear A-log(gsw)
        scaled_A[ir,z] = mu_intercept[z] + 
          (mu_intercept_low_light[z] + b_intercept_low_light_id[id[ir],z]) * (light_treatment[ir] == 1) +
          b_intercept_id[id[ir],z] + b_intercept_error[k,z] + 
          (mu_slope[z] + (mu_slope_low_light[z] + b_slope_low_light_id[id[ir],z]) * (light_treatment[ir] == 1) + 
            b_slope_id[id[ir],z]) * scaled_log_gsw[ir,z];

        // Component 2: quadratic A-log(gsw)
        if (z == 2) {
          scaled_A[ir,z] += 
            (mu_b2 + (mu_b2_low_light + b_b2_low_light_id[id[ir]]) * (light_treatment[ir] == 1) + 
              b_b2_id[id[ir]]) * scaled_log_gsw[ir,z] ^ 2;
        }

        A[ir,z] = scaled_A[ir,z] * sd_A[z] + mean_A[z];
            
        w_i[ir,z] = li6800_svp(T_leaf[ir], P[ir]);
        w_a[ir,z] = RH[ir] * li6800_svp(T_air[ir], P[ir]);
        
        // amphi leaves
        if (leaf_type[ir] == 1) {
          g_tc[ir,z] = li6800_gtc_amphi(g_bw[ir], g_sw[ir,z], K[ir]);
          g_tw[ir,z] = li6800_gtw_amphi(g_bw[ir], g_sw[ir,z], K[ir]);
        } 
        // pseudohypo leaves
        if (leaf_type[ir] == 2) {
          g_tc[ir,z] = li6800_gtc_hypo(g_bw[ir], g_sw[ir,z]);
          g_tw[ir,z] = li6800_gtw_hypo(g_bw[ir], g_sw[ir,z]);
        } 

        E[ir,z]    = li6800_E(g_tw[ir,z], w_a[ir,z], w_i[ir,z]);
        w_0[ir,z]  = li6800_w0(E[ir,z], flow[ir], s[ir], w_a[ir,z]);
        c_0[ir,z]  = li6800_c0(A[ir,z], c_a[ir,z], flow[ir], s[ir], w_a[ir,z], w_0[ir,z]);

      }
    }
  }
  
  // priors ----
  // priors on variable scaling
  sd_A ~ normal(0, 10);
  mean_A ~ normal(15, 10);
  
  // priors on hyperparameters
  b_autocorr_c ~ normal(0, 1); 
  b_autocorr_w ~ normal(0, 1); 
  // sigma_c ~ gamma(2, 0); // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  sigma_c ~ normal(0, 1); 
  sigma_w ~ normal(0, 1); 
  
  mu_intercept ~ normal(0, 10);
  mu_intercept_low_light ~ normal(0, 10);
  sigma_intercept_id ~ normal(0, 10);
  sigma_intercept_low_light_id ~ normal(0, 10);
  sigma_intercept_error ~ normal(0, 10);

  mu_slope ~ normal(0, 10);
  mu_slope_low_light ~ normal(0, 5);
  sigma_slope_id ~ normal(0, 10);
  sigma_slope_low_light_id ~ normal(0, 10);

  // priors on parameters and weights
  mu_b2 ~ normal(0, 10);
  mu_b2_low_light ~ normal(0, 10);
  sigma_b2_id ~ normal(0, 10);
  sigma_b2_low_light_id ~ normal(0, 10);

  b_b2_id ~ normal(0, sigma_b2_id);
  b_b2_low_light_id ~ normal(0, sigma_b2_low_light_id);

  vector[n_comp] alpha;
  for (z in 1:n_comp) {
    alpha[z] = 1;
  }
  
  for (z in 1:n_comp) {

    b_intercept_id[z] ~ normal(0, sigma_intercept_id[z]);
    b_intercept_low_light_id[z] ~ normal(0, sigma_intercept_low_light_id[z]);
    b_slope_id[z] ~ normal(0, sigma_slope_id[z]);
    b_slope_low_light_id[z] ~ normal(0, sigma_slope_low_light_id[z]);

    for (k in 1:n_curve) {
      b_intercept_error[k,z] ~ normal(0, sigma_intercept_error[z]);
      w[k] ~ dirichlet(alpha);
      for (j in 1:n_pts[k]) {
        c_a[sum(n_pts[0:(k - 1)]) + j,z] ~ normal(415, 1);
        log_gsw[sum(n_pts[0:(k - 1)]) + j,z] ~ normal(-1, 1); // is this the right choice?
      }
    }
  }

  // likelihood
  for (k in 1:n_curve) {

    // placeholders
    array[n_comp] real lambda;
    array[n_comp] row_vector[n_pts[k]] x;
    array[n_comp] row_vector[n_pts[k]] y;
    array[n_comp] matrix[n_pts[k], n_pts[k]] R;
    array[n_comp] vector[n_pts[k]] s_vec;
      
    // CO2_r
    for (z in 1:n_comp) {
      lambda[z] = w[k,z];
      for (j in 1:n_pts[k]) {
        x[z,j] = c_0[j,z];
        y[z,j] = CO2_r[j];
        s_vec[z,j] = sigma_c[z];
        for (i in 1:n_pts[k]) {
          R[z,i,j] = exp(-b_autocorr_c[z] * abs(elapsed[sum(n_pts[0:(k - 1)]) + j] - elapsed[sum(n_pts[0:(k - 1)]) + i]));
        }
      }
    }
        
    target += log_mix(
      lambda[1],
      multi_normal_lpdf(y[1] | x[1], quad_form_diag(R[1], s_vec[1])),
      multi_normal_lpdf(y[2] | x[2], quad_form_diag(R[2], s_vec[2]))
    );

    // CO2_s
    for (z in 1:n_comp) {
      lambda[z] = w[k,z];
      for (j in 1:n_pts[k]) {
        x[z,j] = c_a[j,z];
        y[z,j] = CO2_s[j];
        s_vec[z,j] = sigma_c[z];
        for (i in 1:n_pts[k]) {
          R[z,i,j] = exp(-b_autocorr_c[z] * abs(elapsed[sum(n_pts[0:(k - 1)]) + j] - elapsed[sum(n_pts[0:(k - 1)]) + i]));
        }
      }
    }

    target += log_mix(
      lambda[1],
      multi_normal_lpdf(y[1] | x[1], quad_form_diag(R[1], s_vec[1])),
      multi_normal_lpdf(y[2] | x[2], quad_form_diag(R[2], s_vec[2]))
    );

    // H2O_r
    for (z in 1:n_comp) {
      lambda[z] = w[k,z];
      for (j in 1:n_pts[k]) {
        x[z,j] = w_0[j,z];
        y[z,j] = H2O_r[j];
        s_vec[z,j] = sigma_w[z];
        for (i in 1:n_pts[k]) {
          R[z,i,j] = exp(-b_autocorr_w[z] * abs(elapsed[sum(n_pts[0:(k - 1)]) + j] - elapsed[sum(n_pts[0:(k - 1)]) + i]));
        }
      }
    }

    target += log_mix(
      lambda[1],
      multi_normal_lpdf(y[1] | x[1], quad_form_diag(R[1], s_vec[1])),
      multi_normal_lpdf(y[2] | x[2], quad_form_diag(R[2], s_vec[2]))
    );

    // H2O_s
    for (z in 1:n_comp) {
      lambda[z] = w[k,z];
      for (j in 1:n_pts[k]) {
        x[z,j] = w_a[j,z];
        y[z,j] = H2O_s[j];
        s_vec[z,j] = sigma_w[z];
        for (i in 1:n_pts[k]) {
          R[z,i,j] = exp(-b_autocorr_w[z] * abs(elapsed[sum(n_pts[0:(k - 1)]) + j] - elapsed[sum(n_pts[0:(k - 1)]) + i]));
        }
      }
    }

    target += log_mix(
      lambda[1],
      multi_normal_lpdf(y[1] | x[1], quad_form_diag(R[1], s_vec[1])),
      multi_normal_lpdf(y[2] | x[2], quad_form_diag(R[2], s_vec[2]))
    );
        
  }
}

