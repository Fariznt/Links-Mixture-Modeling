data {
  int<lower=1> N;                             // Number of observations
  vector<lower=0>[N] y;                       // Observed survival times
  array[N] int<lower=0, upper=1> censored;    // Censoring indicator (1 = censored, 0 = observed)
  int<lower=1> K;                             // Number of mixture components
}

parameters {
  simplex[K] theta;                // Mixing proportions
  vector<lower=0>[K] alpha;        // Weibull shape parameters
  vector<lower=0>[K] lambda;       // Weibull scale parameters
}

model {
  for (k in 1:K) {
    alpha[k] ~ gamma(2, 1);        // Prior for shape parameters
    lambda[k] ~ gamma(2, 1);       // Prior for scale parameters
  }
  
  theta ~ dirichlet(rep_vector(1.0, K));

  // Likelihood
  for (n in 1:N) {
    array[K] real lps;             // Log-probabilities for each mixture component
    for (k in 1:K) {
      if (censored[n] == 0) {      // If the observation is not censored
        lps[k] = log(theta[k]) + weibull_lpdf(y[n] | alpha[k], lambda[k]);
      } else {                     // If the observation is censored
        lps[k] = log(theta[k]) + weibull_lccdf(y[n] | alpha[k], lambda[k]);
      }
    }
    target += log_sum_exp(lps);    // Mixture model likelihood
  }
}

generated quantities {
  array[N] int<lower=1, upper=K> z;  

  for (n in 1:N) {
    vector[K] lps;
    for (k in 1:K) {
      if (censored[n] == 0) {
        lps[k] = log(theta[k]) + weibull_lpdf(y[n] | alpha[k], lambda[k]);
      } else {
        lps[k] = log(theta[k]) + weibull_lccdf(y[n] | alpha[k], lambda[k]);
      }
    }
    
    // Normalize log probabilities using softmax
    vector[K] prob = softmax(lps);
    
    // Sample mixture component assignment
    z[n] = categorical_rng(prob);
  }
}
