functions{
  real AA(real x, real xc, array[] real theta, array[] real x_r, array[] int x_i) {
    return log((theta[1] + theta[2] * x) / (theta[3] + theta[4] * x));
  }
}
data {
  int<lower=0> n1;
  int<lower=0> n2;
  vector[n1] x1;
  vector[n2] x2;
  vector[n1] y1;
  vector[n2] y2;
  array[1] real d1;
  array[1] int d2;
}
parameters {
  real b0_a;
  real b1_a;
  real b0_h;
  real b1_h;
  real<lower=0> sigma;
}
model {
  y1 ~ normal(b0_a + b1_h * x1, sigma);
  y2 ~ normal(b0_h + b1_a * x2, sigma);
}
generated quantities{
  real A;
  array[4] real theta;
  theta[1] = b0_a;
  theta[2] = b1_a;
  theta[3] = b0_h;
  theta[4] = b1_h;
  A = integrate_1d(AA, 0.5, 1.0, theta, d1, d2);

}

