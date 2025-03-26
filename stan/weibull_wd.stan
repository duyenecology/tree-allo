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
  row_vector<lower=0>[K] gamma_int;
  row_vector[K] gamma_slope;
  matrix[J, K] z;
  vector<lower=0>[K] tau;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[J, K] beta;
  matrix[2, K] gamma = append_row(gamma_int, gamma_slope);

  // Corrected matrix multiplication and addition
  for (k in 1:K) {
    beta[, k] = u * gamma[, k] + tau[k] * z[, k];
  }
}

model {
  vector[N] log_mu;
  to_vector(z) ~ std_normal();
  to_vector(gamma_int) ~ normal(0, 2.5);
  to_vector(gamma_slope) ~ normal(0, 1.25);
  to_vector(tau) ~ cauchy(0, 2.5);
  for (n in 1:N) {
    log_mu[n] = beta[jj[n], 1] + log1m_exp(-beta[jj[n], 2] * pow(x[n], beta[jj[n], 3]));
  }
  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] |
      beta[jj[n], 1] + log1m_exp(-beta[jj[n], 2] * pow(x[n], beta[jj[n], 3])),
      sigma);
  }
}
