data {
  int<lower=0> N;                    // num trees
  int<lower=1> K;                    // num of tree-level predictors
  int<lower=1> J;                    // num sp
  int<lower=1> L;                    // num sp-level predictor
  array[N] int<lower=1, upper=J> jj; // sp indicator for each trees
  matrix[N, K] log_x;                // tree-level predictors
  matrix[L, J] u;                    // sp-level predictors
  vector[N] log_y;                   // outcomes
}

parameters {
  matrix[K, L] gamma;
  matrix[K, J] z;
  vector<lower=0>[K] tau;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, J] beta;
  for (k in 1:K) {
    beta[k, ] = gamma[k, ] * u + tau[k] * z[k, ];
  }
}

model {
  vector[N] log_mu;
  to_vector(z) ~ std_normal();
  to_vector(gamma) ~ normal(0, 2.5);
  tau ~ cauchy(0, 2.5);
  for (n in 1:N) {
    log_mu[n] = log_x[n, ] * beta[, jj[n]];
  }
  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] | log_x[n, ] * beta[, jj[n]], sigma);
  }
}
