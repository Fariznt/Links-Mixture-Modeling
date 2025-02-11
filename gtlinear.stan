data {
  int<lower=1> N; // Number of data points
  int<lower=1> K; // Number of predictors
  matrix[N, K] X; // Predictor matrix
  vector[N] y; // Response vector
}
parameters {
  real<lower=0, upper=1> theta; // Mixture weight for the first component
  real mu1; // Mean of the first component
  real mu2; // Mean of the second component
  real<lower=0> sigma1; // Standard deviation of the first component
  real<lower=0> sigma2; // Standard deviation of the second component
  vector[K] beta1; // Regression coefficients for the first component
  vector[K] beta2; // Regression coefficients for the second component
}
model {
  // Priors
  mu1 ~ normal(0, 5);
  mu2 ~ normal(0, 5);
  sigma1 ~ cauchy(0, 2.5);
  sigma2 ~ cauchy(0, 2.5);
  beta1 ~ normal(0, 5);
  beta2 ~ normal(0, 5);
  theta ~ beta(1, 1);
  
  // Mixture model likelihood
  for (n in 1 : N) {
    target += log_sum_exp(log(theta)
                          + normal_lpdf(y[n] | mu1 + X[n] * beta1, sigma1),
                          log1m(theta)
                          + normal_lpdf(y[n] | mu2 + X[n] * beta2, sigma2));
  }
}
generated quantities {
  array[N] int<lower=1, upper=2> z; // Mixture membership
  for (n in 1 : N) {
    // Calculate unnormalized log probabilities for each component
    real log_prob1 = log(theta)
                     + normal_lpdf(y[n] | mu1 + X[n] * beta1, sigma1);
    real log_prob2 = log1m(theta)
                     + normal_lpdf(y[n] | mu2 + X[n] * beta2, sigma2);
    
    // Normalize probabilities using softmax
    vector[2] prob = softmax([log_prob1, log_prob2]');
    
    // Sample z[n] based on the posterior probabilities
    z[n] = categorical_rng(prob);
  }
}
