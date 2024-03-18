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
