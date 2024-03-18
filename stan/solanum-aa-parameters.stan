parameters {
  vector[n_curve] b0;
  vector[n_curve] b1;
  vector[n_curve] b2;
  real b0_log_sigma_resid;
  real b_log_sigma_resid_S;
  real rho_resid;
  real b0_aa;
  real b_aa_light_intensity_2000;
  real b_aa_light_treatment_high;
  vector[n_acc] b_aa_acc;
  real rho_aa_acc;
  real etasq_aa_acc;
  vector[n_acc_id] b_aa_acc_id;
  real log_sigma_aa_acc_id;
  real b0_log_sigma_aa;
  real b_log_sigma_aa_light_intensity_2000;
  real b_log_sigma_aa_light_treatment_high;
  real b0_ppfd_aa;
  real b1_ppfd_aa;
  real rhosq_ppfd_aa;
  real etasq_ppfd_aa;
  real b_aa_2000_high;
  vector[n_acc] b_aa_light_intensity_2000_acc;
  real rhosq_aa_light_intensity_2000_acc;
  real etasq_aa_light_intensity_2000_acc;
  vector[n_acc] b_aa_light_treatment_high_acc;
  real rhosq_aa_light_treatment_high_acc;
  real etasq_aa_light_treatment_high_acc;
}
transformed parameters {
  real sigma_aa_acc_id;
  sigma_aa_acc_id = exp(log_sigma_aa_acc_id);
}
