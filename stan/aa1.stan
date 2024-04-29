functions {

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
  
  int<lower=0> n;
  int<lower=0> n_acc;
  int<lower=0> n_acc_id;
  int<lower=0> n_light_intensity;  
  int<lower=0> n_light_treatment;  
  
  // variables indexed by lightintensity_x_id groups
  array[n] real scaled_aa;
  array[n] int<lower=1,upper=n_acc> acc;
  array[n] int<lower=1,upper=n_acc_id> acc_id;
  array[n] int<lower=1,upper=n_light_intensity> light_intensity;
  array[n] int<lower=1,upper=n_light_treatment> light_treatment;

  // SPLASH data
  vector[n_acc] scaled_ppfd_mol_m2;

  // distance matrix for Gaussian Process
  matrix[n_acc,n_acc] Dmat;
  
}
parameters {
  real b0_aa;
  real b_aa_light_intensity_2000;
  real b_aa_light_treatment_high;
  vector[n_acc] b_aa_acc;
  real<lower=0> rhosq_aa_acc;
  real<lower=0> etasq_aa_acc;
  vector[n_acc_id] b_aa_acc_id;
  real<lower=0> sigma_aa_acc_id;
  real b0_log_sigma_aa;
  real b_log_sigma_aa_light_intensity_2000;
  real b_log_sigma_aa_light_treatment_high;
}
model {
  // priors on phylogenetic structure
  matrix[n_acc,n_acc] Sigma_aa_acc;
  Sigma_aa_acc = cov_GPL1(Dmat, etasq_aa_acc, rhosq_aa_acc, 0);

  // priors
  b0_aa ~ normal(0,1);
  b_aa_light_intensity_2000 ~ normal(0,1);
  b_aa_light_treatment_high ~ normal(0,1);
  b_aa_acc ~ multi_normal(rep_vector(0.0, n_acc), Sigma_aa_acc);
  rhosq_aa_acc ~ normal(0,10);
  etasq_aa_acc ~ normal(0,10);
  b_aa_acc_id ~ normal(0,sigma_aa_acc_id);
  // log_sigma_aa_acc_id ~ normal(-3,5);
  sigma_aa_acc_id ~ student_t(3, 0, 2.5);
  b0_log_sigma_aa ~ normal(-3,5);
  b_log_sigma_aa_light_intensity_2000 ~ normal(0,1);
  b_log_sigma_aa_light_treatment_high ~ normal(0,1);

  for (i in 1:n) {
    
    real b_2000;
    real b_high;
    real mu1;
    real sigma;
    
    // regression on aa
    b_2000 = b_aa_light_intensity_2000;

    b_high = b_aa_light_treatment_high;
      
    mu1 = b0_aa + 
      b_2000 * (light_intensity[i] == 2) +
      b_high * (light_treatment[i] == 2) +
      b_aa_acc[acc[i]] +
      b_aa_acc_id[acc_id[i]];
    
    // regression on sigma_aa
    sigma = exp(b0_log_sigma_aa +
      b_log_sigma_aa_light_intensity_2000 * (light_intensity[i] == 2) +
      b_log_sigma_aa_light_treatment_high * (light_treatment[i] == 2));

    target += normal_lpdf(scaled_aa[i] | mu1, sigma);

  }

}
generated quantities {
  
  // calculated log-likelihood to estimate LOOIC for model comparison
  vector[n] log_lik;

  for (i in 1:n) {
    
    real b_2000;
    real b_high;
    real mu1;
    real sigma;
    
    // regression on aa
    b_2000 = b_aa_light_intensity_2000;

    b_high = b_aa_light_treatment_high;
      
    mu1 = b0_aa + 
      b_2000 * (light_intensity[i] == 2) +
      b_high * (light_treatment[i] == 2) +
      
      b_aa_acc[acc[i]] +
      b_aa_acc_id[acc_id[i]];
    
    // regression on sigma_aa
    sigma = exp(b0_log_sigma_aa +
      b_log_sigma_aa_light_intensity_2000 * (light_intensity[i] == 2) +
      b_log_sigma_aa_light_treatment_high * (light_treatment[i] == 2));

    log_lik[i] = normal_lpdf(scaled_aa[i] | mu1, sigma);

  }
  
}
