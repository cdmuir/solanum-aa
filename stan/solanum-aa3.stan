// Current working version using LICOR calculations, ragged arrays, but no autocorrelation.
// Currently paused because I am trying out fitting on calculated A and g_sw rather
// than recalculating during estimation
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
  int<lower=0> n_curve; 
  
  int<lower=0> n_id;
  int<lower=0> n_leaf_type;  
  int<lower=0> n_leaftype_x_id;
  int<lower=0> n_light_treatment;  
  int<lower=0> n_light_intensity;  
  
  int<lower=1> n_comp;  
  
  // vector of integer lengths per curve
  array[n_curve] int<lower=0> n_pts; 

  array[n] int<lower=1,upper=n_curve> curve;
  
  array[n] int<lower=1,upper=n_id> id;
  array[n] int<lower=1,upper=n_leaf_type> leaf_type;
  array[n] int<lower=1,upper=n_leaftype_x_id> leaftype_x_id;
  array[n] int<lower=1,upper=n_light_treatment> light_treatment;
  array[n] int<lower=1,upper=n_light_intensity> light_intensity;
  
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
  vector<lower=0>[n_comp] sigma_c; 
  vector<lower=0>[n_comp] sigma_w; 

  vector[n_comp] mu_intercept;
  vector<lower=0>[n_comp] sigma_intercept_id;
  vector<lower=0>[n_comp] sigma_intercept_curve;
  vector<lower=0>[n_comp] sigma_intercept_leaftype_x_id;

  vector[n_comp] mu_slope;
  vector<lower=0>[n_comp] sigma_slope_id;
  vector<lower=0>[n_comp] sigma_slope_curve;

  // 'fixed' effects on intercept
  vector[n_comp] b_intercept_pseudohypo;
  vector[n_comp] b_intercept_high_light;
  vector[n_comp] b_intercept_high_intensity;
  vector[n_comp] b_intercept_pseudohypo__high_intensity;
  vector[n_comp] b_intercept_pseudohypo__high_light;
  vector[n_comp] b_intercept_high_intensity__high_light;
  vector[n_comp] b_intercept_pseudohypo__high_intensity__high_light;

  // 'fixed' effects on intercept
  vector[n_comp] b_slope_pseudohypo;
  vector[n_comp] b_slope_high_light;
  vector[n_comp] b_slope_high_intensity;
  vector[n_comp] b_slope_pseudohypo__high_intensity;
  vector[n_comp] b_slope_pseudohypo__high_light;
  vector[n_comp] b_slope_high_intensity__high_light;
  vector[n_comp] b_slope_pseudohypo__high_intensity__high_light;

  // parameters
  array[n_id,n_comp] real b_intercept_id;
  array[n_curve,n_comp] real b_intercept_curve;
  // This is the 'error' associated with day-to-day changes in leaf physiology,
  // leaf area in the chamber, etc. that affect both high and low intensity curves,
  // but are not actually caused by amphi advantage/disadvantage.
  array[n_leaftype_x_id,n_comp] real b_intercept_leaftype_x_id;

  array[n_id,n_comp] real b_slope_id;
  array[n_curve,n_comp] real b_slope_curve;

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
transformed parameters{

  // placeholders ----
  real intercept;
  real slope;
  
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
    for (ir in 1:n) {

        scaled_log_gsw[ir,z] = (log_gsw[ir,z] - mean_log_gsw[z]) / sd_log_gsw[z];
        
  profile("regression") {
        // Component 1: linear A-log(gsw)
        intercept = mu_intercept[z] + 
          b_intercept_pseudohypo[z] * (leaf_type[ir] == 2) +
          b_intercept_high_light[z] * (light_treatment[ir] == 2) +
          b_intercept_high_intensity[z] * (light_intensity[ir] == 2) +
          b_intercept_pseudohypo__high_intensity[z] *
            (leaf_type[ir] == 2 && light_intensity[ir] == 2) +
          b_intercept_pseudohypo__high_light[z] *
            (leaf_type[ir] == 2 && light_treatment[ir] == 2) +
          b_intercept_high_intensity__high_light[z] *
            (light_intensity[ir] == 2 && light_treatment[ir] == 2) +
          b_intercept_pseudohypo__high_intensity__high_light[z] *
            (leaf_type[ir] == 2 && light_intensity[ir] == 2 && light_treatment[ir] == 2) +
          b_intercept_id[id[ir],z] + b_intercept_curve[curve[ir],z] +
          b_intercept_leaftype_x_id[leaftype_x_id[ir],z];
        
        slope = mu_slope[z] + 
          b_slope_pseudohypo[z] * (leaf_type[ir] == 2) +
          b_slope_high_light[z] * (light_treatment[ir] == 2) +
          b_slope_high_intensity[z] * (light_intensity[ir] == 2) +
          b_slope_pseudohypo__high_intensity[z] * 
            (leaf_type[ir] == 2 && light_intensity[ir] == 2) +
          b_slope_pseudohypo__high_light[z] *
            (leaf_type[ir] == 2 && light_treatment[ir] == 2) +
          b_slope_high_intensity__high_light[z] *
            (light_intensity[ir] == 2 && light_treatment[ir] == 2) +
          b_slope_pseudohypo__high_intensity__high_light[z] * 
            (leaf_type[ir] == 2 && light_intensity[ir] == 2 && light_treatment[ir] == 2) +
          b_slope_id[id[ir],z] + b_slope_curve[curve[ir],z];
  }
        scaled_A[ir,z] = intercept + slope * scaled_log_gsw[ir,z];

        // Component 2: quadratic A-log(gsw)
        if (z == 2) {
          
          scaled_A[ir,z] += 
            (mu_b2 + (mu_b2_low_light + b_b2_low_light_id[id[ir]]) * (light_treatment[ir] == 1) + 
              b_b2_id[id[ir]]) * scaled_log_gsw[ir,z] ^ 2;
              
        }

        A[ir,z] = scaled_A[ir,z] * sd_A[z] + mean_A[z];

  profile("LICOR calcs") {
            
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

  }
model {
  
  profile("priors") {
  // priors ----
  // priors on variable scaling
  sd_A ~ normal(0, 10);
  mean_A ~ normal(15, 10);
  
  // priors on hyperparameters
  // sigma_c ~ gamma(2, 0); // https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations
  sigma_c ~ normal(0, 1); 
  sigma_w ~ normal(0, 1); 
  
  mu_intercept ~ normal(0, 10);
  sigma_intercept_id ~ normal(0, 10);
  sigma_intercept_curve ~ normal(0, 10);
  sigma_intercept_leaftype_x_id ~ normal(0, 10);

  mu_slope ~ normal(0, 10);
  sigma_slope_id ~ normal(0, 10);
  sigma_slope_curve ~ normal(0, 10);

  // priors on 'fixed' effects on intercept
  b_intercept_pseudohypo ~ normal(0, 10);
  b_intercept_high_light ~ normal(0, 10);
  b_intercept_high_intensity ~ normal(0, 10);
  b_intercept_pseudohypo__high_intensity ~ normal(0, 10);
  b_intercept_pseudohypo__high_light ~ normal(0, 10);
  b_intercept_high_intensity__high_light ~ normal(0, 10);
  b_intercept_pseudohypo__high_intensity__high_light ~ normal(0, 10);

  // priors on 'fixed' effects on slope
  b_slope_pseudohypo ~ normal(0, 10);
  b_slope_high_light ~ normal(0, 10);
  b_slope_high_intensity ~ normal(0, 10);
  b_slope_pseudohypo__high_intensity ~ normal(0, 10);
  b_slope_pseudohypo__high_light ~ normal(0, 10);
  b_slope_high_intensity__high_light ~ normal(0, 10);
  b_slope_pseudohypo__high_intensity__high_light ~ normal(0, 10);
  
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
    b_intercept_curve[z] ~ normal(0, sigma_intercept_curve[z]);
    b_intercept_leaftype_x_id[z] ~ normal(0, sigma_intercept_leaftype_x_id[z]);
    b_slope_curve[z] ~ normal(0, sigma_slope_curve[z]);

    for (k in 1:n_curve) {
      // b_intercept_error[k,z] ~ normal(0, sigma_intercept_error[z]);
      w[k] ~ dirichlet(alpha);
      for (j in 1:n_pts[k]) {
        c_a[cumulative_sum(n_pts[1:k])[k] - n_pts[k] + j,z] ~ normal(415, 1);
        log_gsw[cumulative_sum(n_pts[1:k])[k] - n_pts[k] + j,z] ~ normal(-1, 1); // is this the right choice?
      }
    }
  }
}
  // likelihood ----
  for (i in 1:n) {

    // placeholders
    real lambda;

    profile("CO2_r likelihood") {  
    // CO2_r
    target += log_mix(
      w[curve[i], 1],
      normal_lpdf(CO2_r[i] | c_0[i,1], sigma_c[1]),
      normal_lpdf(CO2_r[i] | c_0[i,2], sigma_c[2])
    );
    }
    
    profile("CO2_s likelihood") {  
    // CO2_s
    target += log_mix(
      w[curve[i], 1],
      normal_lpdf(CO2_s[i] | c_a[i,1], sigma_c[1]),
      normal_lpdf(CO2_s[i] | c_a[i,2], sigma_c[2])
    );
    }

    profile("H2O_r likelihood") {  
    // H2O_r
    target += log_mix(
      w[curve[i], 1],
      normal_lpdf(H2O_r[i] | w_0[i,1], sigma_w[1]),
      normal_lpdf(H2O_r[i] | w_0[i,2], sigma_w[2])
    );
    }
    
    profile("H2O_s likelihood") {  
    // H2O_s
    target += log_mix(
      w[curve[i], 1],
      normal_lpdf(H2O_s[i] | w_a[i,1], sigma_w[1]),
      normal_lpdf(H2O_s[i] | w_a[i,2], sigma_w[2])
    );
    }
    
  }
}

