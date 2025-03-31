data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector[N] y;
}

parameters {
  real<lower=0, upper=1> theta;
  real mu1;
  real<lower=0> sigma1;
  vector[K] beta1;
  real mu2;
  real<lower=0> sigma2;
  vector[K] beta2;
}

model {
  theta ~ beta(1, 1);
  mu1 ~ normal(0, 5);
  sigma1 ~ cauchy(0, 2.5);
  beta1 ~ normal(0, 5);
  mu2 ~ normal(0, 5);
  sigma2 ~ cauchy(0, 2.5);
  beta2 ~ normal(0, 5);
for (n in 1:N) {
    target += log_mix(theta,
      normal_lpdf(y[n] | mu1 + X[n] * beta1, sigma1),
      normal_lpdf(y[n] | mu2 + X[n] * beta2, sigma2));
  }
}

generated quantities {
array[N] int<lower=1, upper=2> z;
  for (n in 1:N) {
    real log_prob1 = log(theta) + normal_lpdf(y[n] | mu1 + X[n] * beta1, sigma1);
    real log_prob2 = log1m(theta) + normal_lpdf(y[n] | mu2 + X[n] * beta2, sigma2);
    z[n] = categorical_rng(softmax([log_prob1, log_prob2]'));
  }
}
