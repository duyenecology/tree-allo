data {
  int<lower=0> N;                    // num trees
  int<lower=1> J;                    // num sp
  int<lower=1> K;                    // num of tree-level predictors
  array[N] int<lower=1, upper=J> jj; // sp indicator for each tree
  vector[N] xs;                       // DBH / max(DBH)
  real dbh_max;
  vector[N] log_y;                   // outcomes
}

parameters {
  vector<lower=0>[K] gamma;
  matrix[J, K] z;
  vector<lower=0>[K] tau;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[J, K] beta;
  // Corrected matrix multiplication and addition
  for (k in 1:K) {
    beta[, k] =  gamma[k] + tau[k] * z[, k];
  }
}

model {
  vector[N] log_mu;
  sigma ~ normal(0, 1);
  to_vector(z) ~ std_normal();
  to_vector(tau) ~ normal(0, 1);
  to_vector(gamma) ~ normal(0, 2.5);
  log_mu = beta[jj, 1] + log1m_exp(-beta[jj, 2] .* pow(xs, beta[jj, 3]));
  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  real gamma2 = gamma[2] / pow(dbh_max, gamma[3]);
  vector[J] beta2 = beta[, 2] ./ pow(dbh_max, beta[, 3]);
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] |
      beta[jj[n], 1] + log1m_exp(-beta[jj[n], 2] * pow(xs[n], beta[jj[n], 3])),
      sigma);
  }
}
