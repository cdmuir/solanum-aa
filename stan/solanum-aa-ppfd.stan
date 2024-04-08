  // regression of scaled_ppfd_mol_m2 on aa_acc
  vector[n_acc] aa_acc;
  matrix[n_acc,n_acc] Sigma_ppfd;
  aa_acc = {aa_acc};
  Sigma_ppfd = cov_GPL1(Dmat, etasq_ppfd_aa, rhosq_ppfd_aa, 0);
  aa_acc ~ multi_normal(b0_ppfd_aa + b1_ppfd_aa * scaled_ppfd_mol_m2, Sigma_ppfd);
