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
  int<lower=0> n_rep;
  int<lower=0> n_leaf_type;  
  array[n_pts,n_rep,n_leaf_type] real elapsed;
  array[n_pts,n_rep,n_leaf_type] real flow;
  array[n_pts,n_rep,n_leaf_type] real g_bw;
  array[n_pts,n_rep,n_leaf_type] real K;
  array[n_pts,n_rep,n_leaf_type] real P;
  array[n_pts,n_rep,n_leaf_type] real RH;
  array[n_pts,n_rep,n_leaf_type] real s;
  array[n_pts,n_rep,n_leaf_type] real T_air;
  array[n_pts,n_rep,n_leaf_type] real T_leaf;
  array[n_pts,n_rep,n_leaf_type] real CO2_r;
  array[n_pts,n_rep,n_leaf_type] real CO2_s;
  array[n_pts,n_rep,n_leaf_type] real H2O_r;
  array[n_pts,n_rep,n_leaf_type] real H2O_s;
}
parameters {
  // hyperparameters
  real<lower=0> b_autocorr_c; 
  real<lower=0> b_autocorr_w; 
  real<lower=0> sigma_c; 
  real<lower=0> sigma_w; 
  real mu_intercept;
  real mu_slope;
  real<lower=0> sigma_error_intercept;
  
  // parameters
  array[n_rep,n_leaf_type] real intercept;
  array[n_pts,n_rep,n_leaf_type] real log_gsw;
  array[n_pts,n_rep,n_leaf_type] real log_A;
}
transformed parameters {
  
  array[n_pts,n_rep,n_leaf_type] real g_sw = exp(log_gsw);
  array[n_pts,n_rep,n_leaf_type] real A = exp(log_A);

  // calculated quantities
  array[n_pts,n_rep,n_leaf_type] real w_i = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real g_tc = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real g_tw = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real E = rep_array(0.0, n_pts,n_rep,n_leaf_type);

  // real [CO2] and [H2O]
  array[n_pts,n_rep,n_leaf_type] real c_0 = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real c_a = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real w_0 = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real w_a = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  
  // error covariance matrices
  array[n_pts,n_pts,n_rep,n_leaf_type] real S_c;
  array[n_pts,n_pts,n_rep,n_leaf_type] real S_w;
  
  // calculations
  for (k in 1:n_leaf_type) {
    for (j in 1:n_rep) {
      for (i2 in 1:n_pts) {
        
        w_i[i2,j,k]  += li6800_svp(T_leaf[i2,j,k], P[i2,j,k]);
        w_a[i2,j,k]  += RH[i2,j,k] * li6800_svp(T_air[i2,j,k], P[i2,j,k]);
        
        // amphi leaves
        if (k == 1) {
          g_tc[i2,j,k] += li6800_gtc_amphi(g_bw[i2,j,k], g_sw[i2,j,k], K[i2,j,k]);
          g_tw[i2,j,k] += li6800_gtw_amphi(g_bw[i2,j,k], g_sw[i2,j,k], K[i2,j,k]);
        } 
        // pseudohypo leaves
        if (k == 2) {
          g_tc[i2,j,k] += li6800_gtc_hypo(g_bw[i2,j,k], g_sw[i2,j,k]);
          g_tw[i2,j,k] += li6800_gtw_hypo(g_bw[i2,j,k], g_sw[i2,j,k]);
        } 

        E[i2,j,k]    += li6800_E(g_tw[i2,j,k], w_a[i2,j,k], w_i[i2,j,k]);
        w_0[i2,j,k]  += li6800_w0(E[i2,j,k], flow[i2,j,k], s[i2,j,k], w_a[i2,j,k]);
        c_0[i2,j,k]  += li6800_c0(A[i2,j,k], c_a[i2,j,k], flow[i2,j,k], s[i2,j,k], w_a[i2,j,k], w_0[i2,j,k]);

        for (i1 in 1:n_pts) {
          S_c[i1, i2, j, k] = sigma_c * exp(-b_autocorr_c * abs(elapsed[i2, j, k] - elapsed[i1, j, k]));
          S_w[i1, i2, j, k] = sigma_w * exp(-b_autocorr_w * abs(elapsed[i2, j, k] - elapsed[i1, j, k]));

        }
      }
    }
  }
}
model {
  
  // placeholders
  array[n_pts,n_rep,n_leaf_type] real mu = rep_array(0.0, n_pts,n_rep,n_leaf_type);
  // brms syntax (change to my own once I understand)
  // array storing lagged residuals
  array[n_pts,n_rep,n_leaf_type] real Err = rep_array(0, n_pts,n_rep,n_leaf_type);
  array[n_pts,n_rep,n_leaf_type] real err;  // actual residuals
  
  // priors on hyperparameters
  b_autocorr_c ~ normal(0, 1); 
  b_autocorr_w ~ normal(0, 1); 
  sigma_c ~ gamma(2, 0); // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  sigma_w ~ gamma(2, 0); 
  mu_intercept ~ normal(0, 10);
  mu_slope ~ normal(0, 10);
  sigma_error_intercept ~ gamma(2, 0);

  // priors on parameters (WHAT GOES HERE?)
  // array[n_rep,n_leaf_type] real intercept;
  // array[n_pts,n_rep,n_leaf_type] real log_gsw;
  // array[n_pts,n_rep,n_leaf_type] real log_A;
  // array[n_pts,n_rep,n_leaf_type] real c_0;
  // array[n_pts,n_rep,n_leaf_type] real c_a;
  // array[n_pts,n_rep,n_leaf_type] real w_0;
  // array[n_pts,n_rep,n_leaf_type] real w_a;
  
  // likelihood
  for (k in 1:n_leaf_type) {
    intercept[,k] ~ normal(0, sigma_error_intercept);
  }
  
  // model
  for (k in 1:n_leaf_type) {
    for (j in 1:n_rep) {
      
      // placeholders
      row_vector[n_pts] x;
      row_vector[n_pts] y;
      matrix[n_pts, n_pts] S;
      
      // CO2_r
      for (i2 in 1:n_pts) {
        x[i2] = c_0[i2,j,k];
        y[i2] = CO2_r[i2,j,k];
        for (i1 in 1:n_pts) {
          S[i1, i2] = S_c[i1, i2, j, k];
        }
      }
      y ~ multi_normal(x, S);

      // mu[1,j,k] += mu_intercept + intercept[j,k] + slope * log_gsw[1,j,k];
      // err[1,j,k] = log_A[1,j,k] - mu[1,j,k];
      // y = log_A[1,j,k];
      // y ~ normal(mu[1,j,k], sigma_error_resid);
      //   
      // for (i in 2:n_pts) {
      //   mu[i,j,k] += mu_intercept + intercept[j,k] + slope * log_gsw[i,j,k];
      //   err[i,j,k] = log_A[i,j,k] - mu[i,j,k];
      //   Err[i,j,k] = err[i-1,j,k];
      //   mu[i,j,k] += Err[i,j,k] * rho_error_resid;
      //   y = log_A[i,j,k];
      //   y ~ normal(mu[i,j,k], sigma_error_resid);
      // }
      
    }
  }
  
}

