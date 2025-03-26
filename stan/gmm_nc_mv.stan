data {
  int<lower=0> N;                    // num trees
  int<lower=1> J;                    // num sp
  int<lower=1> K;                    // num of tree-level predictors
  array[N] int<lower=1, upper=J> jj; // sp indicator for each tree
  vector[N] x;                       // DBH (non-log-scale)
  vector[N] log_y;                   // outcomes
}

transformed data {
  matrix[1, J] u = rep_matrix(1, 1, J);
}

parameters {
  matrix[K, 1] gamma;
  matrix[K, J] z;
  cholesky_factor_corr[K] L_Omega;
  vector<lower=0>[K] tau;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, J] beta;
  beta = gamma * u + diag_pre_multiply(tau, L_Omega) * z;
}

model {
  vector[N] log_mu;
  to_vector(z) ~ std_normal();
  to_vector(tau) ~ cauchy(0, 1);
  to_vector(gamma) ~ normal(0, 5);
  // log_mu = beta[1, jj] + beta[2, jj] .* log(x) -
  //     log(beta[3, jj] + pow(x, beta[2, jj]));
  log_mu = to_vector(beta[1, jj]) +
    to_vector(beta[2, jj]) .* log(x) -
    log(to_vector(beta[3, jj]) + pow(x, to_vector(beta[2, jj])));

  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] |
      beta[1, jj[n]] +
      beta[2, jj[n]] * log(x[n]) -
      log(beta[3, jj[n]] + pow(x[n], beta[2, jj[n]])),
      sigma);
  }
}
