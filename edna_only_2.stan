data {
  int<lower=1> S;
  int<lower=1> N;
  vector[N] y;
  int<lower=1, upper=S> site_y[N];

}

parameters {
  vector[S] X;

  // eDNA
  real alpha;
  real<lower=0> beta;
  real<lower=0> sigma;  

}

transformed parameters {
  vector[S] mu;

  mu      = alpha + beta * X;
}

model {
  // Latent field
  X ~ normal(0, 1);

  // eDNA priors
  alpha ~ normal(0, 1);
  beta  ~ normal(0, 0.5);
  sigma ~ normal(0, 0.5);  

  // Likelihoods
  y ~ normal(mu[site_y], sigma);
}
