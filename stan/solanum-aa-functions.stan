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