  // regression of scaled_pai on aa_acc
  vector[n_acc] aa_acc;
  matrix[n_acc,n_acc] Sigma_pai;
  aa_acc = {aa_acc};
  Sigma_pai = cov_GPL1(Dmat, etasq_pai_aa, rhosq_pai_aa, 0);
  aa_acc ~ multi_normal(b0_pai_aa + b1_pai_aa * scaled_pai, Sigma_pai);
