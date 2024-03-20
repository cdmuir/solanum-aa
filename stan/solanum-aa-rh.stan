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
