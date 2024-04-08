functions {
  // indefinite integral of log(A_amphi) - log(A_hypo) based on parameters theta
  real aa_int(real x, array[] real theta) {
    
    return x ^ 3 * (theta[3] / 3 - theta[6] / 3) + 
      x ^ 2 * (theta[2] / 2 - theta[5] / 2) + 
      x * (theta[1] - theta[4]);
    
  }
  
  // OU process from McElreath Rethinking v2
  matrix cov_GPL1(matrix x, real sq_alpha, real sq_rho, real delta) {
    int N = dims(x)[1];
    matrix[N, N] K;
    for (i in 1:(N-1)) {
      K[i, i] = sq_alpha + delta;
      for (j in (i + 1):N) {
        K[i, j] = sq_alpha * exp(-sq_rho * x[i,j] );
        K[j, i] = K[i, j];
      }
    }
      K[N, N] = sq_alpha + delta;
      return K;
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
  array[n_lightintensity_x_acc_id] int<lower=1,upper=n_light_intensity> light_intensity;
  array[n_lightintensity_x_acc_id] int<lower=1,upper=n_light_treatment> light_treatment;

  // min and max scaled_log_gsw by curve
  array[n_curve] real min_scaled_log_gsw;
  array[n_curve] real max_scaled_log_gsw;
  array[n_curve] real S;
  
  // SPLASH data
  vector[n_acc] scaled_ppfd_mol_m2;

  // distance matrix for Gaussian Process
  matrix[n_acc,n_acc] Dmat;
  
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
  real b0_log_sigma_resid;
  real b_log_sigma_resid_S;
  real<lower=-1, upper=1> rho_resid;
  real b0_aa;
  real b_aa_light_intensity_2000;
  real b_aa_light_treatment_high;
  vector[n_acc] b_aa_acc;
  real<lower=0> rhosq_aa_acc;
  real<lower=0> etasq_aa_acc;
  vector[n_acc_id] b_aa_acc_id;
  real log_sigma_aa_acc_id;
  real b0_log_sigma_aa;
  real b_log_sigma_aa_light_intensity_2000;
  real b_log_sigma_aa_light_treatment_high;
  real b0_ppfd_aa;
  real b1_ppfd_aa;
  real<lower=0> rhosq_ppfd_aa;
  real<lower=0> etasq_ppfd_aa;
  vector[n_acc] b_aa_light_treatment_high_acc;
  real<lower=0> rhosq_aa_light_treatment_high_acc;
  real<lower=0> etasq_aa_light_treatment_high_acc;
}
transformed parameters {
  real sigma_aa_acc_id;
  vector[3] sigma_curve;
  sigma_aa_acc_id = exp(log_sigma_aa_acc_id);
  sigma_curve = exp(log_sigma_curve);
  matrix[3,3] Sigma_curve;
  Sigma_curve = quad_form_diag(R_curve, sigma_curve);
}
model {
  // priors on phylogenetic structure
  matrix[n_acc,n_acc] Sigma_aa_acc;
  matrix[n_acc,n_acc] Sigma_aa_light_treatment_high_acc;
  Sigma_aa_acc = cov_GPL1(Dmat, etasq_aa_acc, rhosq_aa_acc, 0);
  Sigma_aa_light_treatment_high_acc = cov_GPL1(Dmat, etasq_aa_light_treatment_high_acc, rhosq_aa_light_treatment_high_acc, 0);

  // priors
  B_curve ~ multi_normal(Mu_curve, Sigma_curve);
  Mu_curve ~ normal(0, 10);
  log_sigma_curve ~ normal(0,1);
  R_curve ~ lkj_corr(2);
  b0_log_sigma_resid ~ normal(-3,5);
  b_log_sigma_resid_S ~ normal(0,1);
  rho_resid ~ normal(0,1);
  b0_aa ~ normal(0,1);
  b_aa_light_intensity_2000 ~ normal(0,1);
  b_aa_light_treatment_high ~ normal(0,1);
  b_aa_acc ~ multi_normal(rep_vector(0.0, n_acc), Sigma_aa_acc);
  rhosq_aa_acc ~ normal(0,10);
  etasq_aa_acc ~ normal(0,10);
  b_aa_acc_id ~ normal(0,sigma_aa_acc_id);
  log_sigma_aa_acc_id ~ normal(-3,5);
  b0_log_sigma_aa ~ normal(-3,5);
  b_log_sigma_aa_light_intensity_2000 ~ normal(0,1);
  b_log_sigma_aa_light_treatment_high ~ normal(0,1);
  b0_ppfd_aa ~ normal(0,1);
  b1_ppfd_aa ~ normal(0,1);
  rhosq_ppfd_aa ~ normal(0,10);
  etasq_ppfd_aa ~ normal(0,10);
  b_aa_light_treatment_high_acc ~ multi_normal(rep_vector(0.0, n_acc), Sigma_aa_light_treatment_high_acc);
  rhosq_aa_light_treatment_high_acc ~ normal(0, 10);
  etasq_aa_light_treatment_high_acc ~ normal(0, 10);

  // Estimate AA
  for (i in 1:n_lightintensity_x_acc_id) {
    
    real b_2000;
    real b_high;
    real mu1;
    real sigma;
    
    // regression on aa
    b_2000 = b_aa_light_intensity_2000;

    b_high = b_aa_light_treatment_high + b_aa_light_treatment_high_acc[acc[i]];
      
    mu1 = b0_aa + 
      b_2000 * (light_intensity[i] == 2) +
      b_high * (light_treatment[i] == 2) +
      
      b_aa_acc[acc[i]] +
      b_aa_acc_id[acc_id[i]];
    
    // regression on sigma_aa
    sigma = exp(b0_log_sigma_aa +
      b_log_sigma_aa_light_intensity_2000 * (light_intensity[i] == 2) +
      b_log_sigma_aa_light_treatment_high * (light_treatment[i] == 2));

    real aa_i;
    int amphi_curve; 
    int pseudohypo_curve;
    
    real a; // min of amphi curve
    real b; // max of hypo curve
    array[6] real theta; // parameters
    
    amphi_curve = amphi[i];
    pseudohypo_curve = pseudohypo[i];
    
    a = min_scaled_log_gsw[amphi_curve];
    b = max_scaled_log_gsw[pseudohypo_curve];
    
    theta[1] = Mu_curve[1] + B_curve[curve[amphi_curve],1]; // b0_amphi;
    theta[2] = Mu_curve[2] + B_curve[curve[amphi_curve],2]; // b1_amphi;
    theta[3] = Mu_curve[3] + B_curve[curve[amphi_curve],3]; // b2_amphi;
    theta[4] = Mu_curve[1] + B_curve[curve[pseudohypo_curve],1]; // b0_hypo;
    theta[5] = Mu_curve[2] + B_curve[curve[pseudohypo_curve],2]; // b1_hypo;
    theta[6] = Mu_curve[3] + B_curve[curve[pseudohypo_curve],3]; // b2_hypo;
    
    aa_i = aa_int(b, theta) - aa_int(a, theta);

    target += normal_lpdf(aa_i | mu1, sigma);

  }
  // Estimate RH curve parameters
  vector[n] resid;
  vector[n] mu2;
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
  // regression of scaled_ppfd_mol_m2 on aa_acc
  vector[n_acc] aa_acc;
  matrix[n_acc,n_acc] Sigma_ppfd;
  aa_acc = b0_aa + b_aa_acc + b_aa_light_treatment_high + b_aa_light_treatment_high_acc;;
  Sigma_ppfd = cov_GPL1(Dmat, etasq_ppfd_aa, rhosq_ppfd_aa, 0);
  aa_acc ~ multi_normal(b0_ppfd_aa + b1_ppfd_aa * scaled_ppfd_mol_m2, Sigma_ppfd);
}
generated quantities {
  // calculated log-likelihood to estimate LOOIC for model comparison
  vector[n_lightintensity_x_acc_id] log_lik;
  // Estimate AA
  for (i in 1:n_lightintensity_x_acc_id) {
    
    real b_2000;
    real b_high;
    real mu1;
    real sigma;
    
    // regression on aa
    b_2000 = b_aa_light_intensity_2000;

    b_high = b_aa_light_treatment_high + b_aa_light_treatment_high_acc[acc[i]];
      
    mu1 = b0_aa + 
      b_2000 * (light_intensity[i] == 2) +
      b_high * (light_treatment[i] == 2) +
      
      b_aa_acc[acc[i]] +
      b_aa_acc_id[acc_id[i]];
    
    // regression on sigma_aa
    sigma = exp(b0_log_sigma_aa +
      b_log_sigma_aa_light_intensity_2000 * (light_intensity[i] == 2) +
      b_log_sigma_aa_light_treatment_high * (light_treatment[i] == 2));

    real aa_i;
    int amphi_curve; 
    int pseudohypo_curve;
    
    real a; // min of amphi curve
    real b; // max of hypo curve
    array[6] real theta; // parameters
    
    amphi_curve = amphi[i];
    pseudohypo_curve = pseudohypo[i];
    
    a = min_scaled_log_gsw[amphi_curve];
    b = max_scaled_log_gsw[pseudohypo_curve];
    
    theta[1] = Mu_curve[1] + B_curve[curve[amphi_curve],1]; // b0_amphi;
    theta[2] = Mu_curve[2] + B_curve[curve[amphi_curve],2]; // b1_amphi;
    theta[3] = Mu_curve[3] + B_curve[curve[amphi_curve],3]; // b2_amphi;
    theta[4] = Mu_curve[1] + B_curve[curve[pseudohypo_curve],1]; // b0_hypo;
    theta[5] = Mu_curve[2] + B_curve[curve[pseudohypo_curve],2]; // b1_hypo;
    theta[6] = Mu_curve[3] + B_curve[curve[pseudohypo_curve],3]; // b2_hypo;
    
    aa_i = aa_int(b, theta) - aa_int(a, theta);

    log_lik[i] = normal_lpdf(aa_i | mu1, sigma);

  }
}
