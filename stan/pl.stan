data {
  int<lower=0> N;                    // num trees
  int<lower=1> K;                    // num of tree-level predictors
  int<lower=1> J;                    // num sp
  array[N] int<lower=1, upper=J> jj; // sp indicator for each trees
  matrix[N, K] log_x;                // tree-level predictors
  vector[N] log_y;                   // outcomes
}

parameters {
  vector[K] gamma;
  matrix[K, J] z;
  vector<lower=0>[K] tau;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, J] beta;
  for (k in 1:K) {
    beta[k, ] = gamma[k] + tau[k] * z[k, ];
  }
}

model {
  vector[N] log_mu;
  sigma ~ normal(0, 1);
  to_vector(z) ~ std_normal();
  tau ~ cauchy(0, 2.5);
  gamma ~ normal(0, 2.5);
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
