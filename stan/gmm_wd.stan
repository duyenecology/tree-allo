data {
  int<lower=0> N;                    // num trees
  int<lower=1> J;                    // num sp
  int<lower=1> K;                    // num of tree-level predictors
  int<lower=1> L;                    // num sp-level predictor
  array[N] int<lower=1, upper=J> jj; // sp indicator for each tree
  vector[N] x;                       // DBH (non-log-scale)
  matrix[J, L] u;                    // sp-level predictors
  vector[N] log_y;                   // outcomes
}

parameters {
  vector<lower=0>[K] gamma_int_hat;
  vector[K] gamma_slope_hat;
  matrix[J, K] z;
  vector<lower=0>[K] tau;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[J, K] beta;
  row_vector[K] gamma_int;
  row_vector[K] gamma_slope;
  gamma_int[1] = gamma_int_hat[1] * 10;
  gamma_int[2] = gamma_int_hat[2];
  gamma_int[3] = gamma_int_hat[3] * 1000;
  gamma_slope[1] = gamma_slope_hat[1] * 10;
  gamma_slope[2] = gamma_slope_hat[2];
  gamma_slope[3] = gamma_slope_hat[3] * 1000;
  matrix[2, K] gamma = append_row(gamma_int, gamma_slope);

  // Corrected matrix multiplication and addition
  for (k in 1:K) {
    beta[, k] = u * gamma[, k] + tau[k] * z[, k];
  }
}

model {
  vector[N] log_mu;
  to_vector(z) ~ std_normal();
  tau ~ cauchy(0, 2.5);
  gamma_int_hat ~ normal(0, 2.5);
  gamma_slope_hat ~ normal(0, 1.25);
  for (n in 1:N) {
    log_mu[n] = beta[jj[n], 1] + beta[jj[n], 2] * log(x[n]) -
      log(beta[jj[n], 3] + pow(x[n], beta[jj[n], 2]));
  }
  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] |
      beta[jj[n], 1] +
      beta[jj[n], 2] * log(x[n]) -
      log(beta[jj[n], 3] + pow(x[n], beta[jj[n], 2])),
      sigma);
  }
}
