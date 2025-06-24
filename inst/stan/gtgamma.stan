data {
  int<lower=1> N; // Number of observations
  int<lower=1> K; // Number of predictors
  vector<lower=0>[N] y; // Response variable (positive values)
  matrix[N, K] X; // Predictor matrix
  
  // prior hyperparameters
  real beta1_loc;
  real beta2_loc;
  real<lower=0> beta1_scale;
  real<lower=0> beta2_scale;
  real<lower=0> theta_alpha;
  real<lower=0> theta_beta;
  real<lower=0> phi1_rate;
  real<lower=0> phi2_rate;
}
parameters {
  real<lower=0, upper=1> theta; // Mixing proportions (must sum to 1)
  vector[K] beta1; // Regression coefficients for component 1
  vector[K] beta2; // Regression coefficients for component 2
  real<lower=0> phi1; // Shape parameter for component 1
  real<lower=0> phi2; // Shape parameter for component 2
}
model {
  vector[N] mu1 = exp(X * beta1); // Mean of Gamma for component 1
  vector[N] mu2 = exp(X * beta2); // Mean of Gamma for component 2
  vector[N] log_lik1;
  vector[N] log_lik2;
  
  // Calculate log-likelihoods for each component
  for (n in 1 : N) {
    log_lik1[n] = gamma_lpdf(y[n] | phi1, phi1 / mu1[n]);
    log_lik2[n] = gamma_lpdf(y[n] | phi2, phi2 / mu2[n]);
    target += log_mix(theta, log_lik1[n], log_lik2[n]);
  }
  
  // Priors for regression coefficients and mix proportion
  beta1 ~ normal(beta1_loc, beta1_scale);
  beta2 ~ normal(beta2_loc, beta2_scale);
  theta ~ beta(theta_alpha, theta_beta);
  // Priors for shape parameters
  phi1 ~ exponential(phi1_rate);
  phi2 ~ exponential(phi2_rate);
}
generated quantities {
  array[N] int<lower=1, upper=2> z; // Mixture membership
  vector[N] mu1 = exp(X * beta1); // Recompute mu1 for generated quantities
  vector[N] mu2 = exp(X * beta2); // Recompute mu2 for generated quantities
  for (n in 1 : N) {
    // Calculate unnormalized log probabilities for each component
    real log_prob1 = log(theta) + gamma_lpdf(y[n] | phi1, phi1 / mu1[n]);
    real log_prob2 = log1m(theta) + gamma_lpdf(y[n] | phi2, phi2 / mu2[n]);
    
    // Normalize probabilities using softmax
    vector[2] prob = softmax([log_prob1, log_prob2]');
    
    // Sample z[n] based on the posterior probabilities
    z[n] = categorical_rng(prob);
  }
}
