data {
  int<lower=0> N;                    // num trees
  int<lower=1> J;                    // num sp
  int<lower=1> K;                    // num of tree-level predictors
  int<lower=1> S;                    // num of biomes
  array[N] int<lower=1, upper=J> jj; // sp indicator for each tree
  array[N] int<lower=1, upper=S> ss; // biome indicator for each tree
  vector[N] x;                       // DBH (non-log-scale)
  vector[N] log_y;                   // outcomes
}

parameters {
  vector<lower=0>[K] gamma_hat;
  matrix[J, K] z_sp;
  matrix[S, K] z_biome;
  vector<lower=0>[K] tau_sp;
  vector<lower=0>[K] tau_biome;
  real<lower=0> sigma;
}

transformed parameters {
  vector<lower=0>[K] gamma;
  matrix[J, K] beta_sp;
  matrix[S, K] beta_biome;

  gamma[1] = gamma_hat[1] * 10;
  gamma[2] = gamma_hat[2];
  gamma[3] = gamma_hat[3] * 1000;

  // Corrected matrix multiplication and addition
  for (k in 1:K) {
    beta_sp[, k] = tau_sp[k] * z_sp[, k];
    beta_biome[, k] = tau_biome[k] * z_biome[, k];
  }
}

model {
  vector[N] log_mu;
  gamma_hat ~ normal(0, 2.5);
  sigma ~ normal(0, 1);

  to_vector(z_sp) ~ std_normal();
  tau_sp ~ cauchy(0, 2.5);

  to_vector(z_biome) ~ std_normal();
  tau_biome ~ cauchy(0, 2.5);

  log_mu = (gamma[1] + beta_sp[jj, 1] + beta_biome[ss, 1]) +
    gamma[2] + beta_sp[jj, 2] + beta_biome[ss, 2] .* log(x) -
      log(gamma[3] + beta_sp[jj, 3] + beta_biome[ss, 3] +
        pow(x, gamma[2] + beta_sp[jj, 2] + beta_biome[ss, 2]));

  log_y ~ normal(log_mu, sigma);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(log_y[n] | (gamma[1] + beta_sp[jj[n], 1] + beta_biome[ss[n], 1]) +
                (gamma[2] + beta_sp[jj[n], 2] + beta_biome[ss[n], 2]) * log(x[n]) -
                log(gamma[3] + beta_sp[jj[n], 3] + beta_biome[ss[n], 3]) +
                pow(x[n], gamma[2] + beta_sp[jj[n], 2] + beta_biome[ss[n], 2]), sigma);
  }
}
