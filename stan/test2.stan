data {
  int<lower=0> N;
  vector[N] x;
  vector[N] y;
}

parameters {
  vector[2] beta;
  real<lower=0> sigma;
}

model {
  y ~ normal(beta[1] + beta[2] * x, sigma);
  beta ~ normal(0, 100);
}

generated quantities {
  vector[N] log_lik;
  for (n in 1:N) {
    log_lik[n] = normal_lpdf(y[n] | beta[1] + beta[2] * x, sigma);
  }
}
