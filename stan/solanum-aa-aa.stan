  // Estimate AA
  for (i in 1:n_lightintensity_x_acc_id) {
    
    real b_2000;
    real b_high;
    real mu1;
    real sigma;
    
    // regression on aa
    b_2000 = {b_2000};

    b_high = {b_high};
      
    mu1 = b0_aa + 
      b_aa_amphi_first * (amphi_first[i]) +
      b_2000 * (light_intensity[i] == 2) +
      b_high * (light_treatment[i] == 2) +
      {b_aa_2000_high}
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
