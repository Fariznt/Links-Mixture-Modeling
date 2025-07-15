data {
  int<lower=1> N;
  int<lower=1> K;
  matrix[N, K] X;
  vector<lower=0>[N] y;
  array[N] int<lower=0, upper=1> status;
  int<lower=0, upper=1> use_truncation;
  real<lower=0> truncation_lower;
  real<lower=0> truncation_upper;
}
parameters {
  vector[K] beta1;
  vector[K] beta2;
  real<lower=0> shape1;
  real<lower=0> shape2;
  real<lower=0> scale1;
  real<lower=0> scale2;
  real<lower=0, upper=1> theta;
}
model {
  beta1 ~ normal(0, 2);
  beta2 ~ normal(0, 2);
  shape1 ~ gamma(2, 1);
  shape2 ~ gamma(2, 1);
  scale1 ~ gamma(2, 1);
  scale2 ~ gamma(2, 1);
  theta ~ beta(1, 1);
  
  for (n in 1 : N) {
    real log_p1 = log(theta);
    real log_p2 = log(1 - theta);
    real linpred1 = X[n] * beta1;
    real linpred2 = X[n] * beta2;
    real scale_factor1 = scale1 * exp(-linpred1);
    real scale_factor2 = scale2 * exp(-linpred2);
    
    if (status[n] == 1) {
      log_p1 += weibull_lpdf(y[n] | shape1, scale_factor1);
      log_p2 += weibull_lpdf(y[n] | shape2, scale_factor2);
    } else {
      log_p1 += weibull_lccdf(y[n] | shape1, scale_factor1);
      log_p2 += weibull_lccdf(y[n] | shape2, scale_factor2);
    }
    target += log_sum_exp(log_p1, log_p2);
  }
}
generated quantities {
  array[N] int<lower=1, upper=2> z;
  vector[N] log_lik;
  
  for (n in 1 : N) {
    vector[2] log_p;
    real linpred1 = X[n] * beta1;
    real linpred2 = X[n] * beta2;
    real scale_factor1 = scale1 * exp(-linpred1);
    real scale_factor2 = scale2 * exp(-linpred2);
    
    log_p[1] = log(theta) + weibull_lpdf(y[n] | shape1, scale_factor1);
    log_p[2] = log(1 - theta) + weibull_lpdf(y[n] | shape2, scale_factor2);
    
    z[n] = categorical_rng(softmax(log_p));
    log_lik[n] = log_sum_exp(log_p[1], log_p[2]);
  }
}