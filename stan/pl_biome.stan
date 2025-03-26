data {
  int<lower=0> N;                    // num trees
  int<lower=1> K;                    // num of tree-level predictors
  int<lower=1> J;                    // num sp
  int<lower=1> S;                    // num of biomes
  array[N] int<lower=1, upper=J> jj; // sp indicator for each trees
  array[N] int<lower=1, upper=S> ss; // biome indicator for each tree
  matrix[N, K] log_x;                // tree-level predictors
  vector[N] log_y;                   // outcomes
}

parameters {
  vector[K] gamma;
  matrix[K, J] z_sp;
  matrix[K, S] z_biome;
  vector<lower=0>[K] tau_sp;
  vector<lower=0>[K] tau_biome;
  real<lower=0> sigma;
}

transformed parameters {
  matrix[K, J] beta_sp;
  matrix[K, S] beta_biome;
  // Corrected matrix multiplication and addition
  for (k in 1:K) {
    beta_sp[k, ] = tau_sp[k] * z_sp[k, ];
    beta_biome[k, ] = tau_biome[k] * z_biome[k, ];
  }
}

model {
  vector[N] log_mu;
  sigma ~ normal(0, 1);

  to_vector(z_sp) ~ std_normal();
  tau_sp ~ cauchy(0, 2.5);

  to_vector(z_biome) ~ std_normal();
  tau_biome ~ cauchy(0, 2.5);

  to_vector(gamma) ~ normal(0, 2.5);
  for (n in 1:N) {
    log_mu[n] = log_x[n, ] * (gamma + beta_sp[, jj[n]] + beta_biome[, ss[n]]);
  }
  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] | log_x[n, ] * (gamma + beta_sp[, jj[n]] + beta_biome[, ss[n]]), sigma);
  }
}
