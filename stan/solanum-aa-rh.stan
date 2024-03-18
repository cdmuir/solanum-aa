  // Estimate RH curve parameters
  vector[n] resid;
  vector[n] mu2;
  vector[n] sigma2;
  
  for (i in 1:n) {
      
    mu2[i] = b0[curve[i]] + 
      b1[curve[i]] * scaled_log_gsw[i] +
      b2[curve[i]] * scaled_log_gsw[i] ^ 2;
    
    sigma2[i] = exp(b0_log_sigma_resid + (6 - S[curve[i]]) * b_log_sigma_resid_S);
  
  }
    
  resid[1] = log_A[1] - mu2[1];
  for (i in 2:n) {
    resid[i] = log_A[i] - mu2[i];
    mu2[i] += rho_resid * resid[i - 1] * (curve[i] == curve[i - 1]);
  }
    
  target += normal_lpdf(log_A | mu2, sigma2);
