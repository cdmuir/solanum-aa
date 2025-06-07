source("r/header.R")


df_parameters = read_csv("data/parameters1.csv", col_types = "cccc", n_max = 8) 

min(which(df_parameters$par_type == "ags"))

read_csv("data/parameters1.csv", col_types = "cccc", n_max = 6) 

  
  # 1. Estimate \agcurve{} curve parameters for each leaf
    # B_curve, n_curve x 3 array of random \agcurve{} curve-level coefficients (b_{0,j}, b_{1,j}, b_{2,j}). B_curve ~ MVN(0, S_curve)
    # Mu_curve, vector of mean quadratic coefficients (\Beta_0, \Beta_1, \Beta_2) 
    # 3 S_curve, covariance matrix of B_curve                  
    # 5 b0_log_sigma_resid       
    # 6 b_log_sigma_resid_S      
    # 7 rho_resid                
    2. Estimate \aax{} for each light intensity with leaf using \agcurve{} curve parameters
  3. Estimate the effects of light intensity, light treatment, and accession on \aax{} (assimilatory and plasticity hypotheses)
  4. Estimate the effects of native light habitat on accession-level \aax{} (constitutive hypothesis)
  
