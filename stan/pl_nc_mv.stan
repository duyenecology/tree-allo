data {
  int<lower=0> N;                    // num trees
  int<lower=1> K;                    // num of tree-level predictors
  int<lower=1> J;                    // num sp
  array[N] int<lower=1, upper=J> jj; // sp indicator for each trees
  matrix[N, K] log_x;                // tree-level predictors
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
  L_Omega ~ lkj_corr_cholesky(2);
  to_vector(gamma) ~ normal(0, 2.5);
  for (n in 1:N) {
    log_mu[n] = dot_product(log_x[n], beta[, jj[n]]);
  }
  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] | dot_product(log_x[n], beta[, jj[n]]), sigma);
  }
}
