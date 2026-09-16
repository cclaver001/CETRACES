data {
  int<lower=1> S;
  int<lower=1> N;
  vector[N] y;
  int<lower=1, upper=S> site_y[N];

  int<lower=1> M;
  vector[M] z;
  int<lower=1, upper=S> site_z[M];
}

parameters {
  vector[S] X;

  // eDNA
  real alpha;
  real<lower=0> beta;
  real<lower=0> sigma;  

  // Sightings
  real kappa;
  real<lower=0> omega;
  real<lower=0> tau;
}

transformed parameters {
  vector[S] mu;
  vector[S] epsilon;

  mu      = alpha + beta * X;
  epsilon = kappa + omega * X;
}

model {
  // Latent field
  X ~ normal(0, 1);

  // eDNA priors
  alpha ~ normal(0, 1);
  //beta  ~ normal(0, 0.5);
  beta ~ normal(0.1184445, 0.08110615);  // based on the eDNA model
  sigma ~ normal(0, 0.5);  

  // Sightings priors
  kappa ~ normal(0, 1);
  omega ~ normal(0, 0.5);
  tau   ~ exponential(1);

  // Likelihoods
  y ~ normal(mu[site_y], sigma);
  z ~ normal(epsilon[site_z], tau);
}

generated quantities {
  vector[S] ratio_cell;
  vector[N] y_rep;
  vector[M] z_rep;

  for (s in 1:S)
    ratio_cell[s] = exp(epsilon[s]) - 1;

  for (n in 1:N)
    y_rep[n] = normal_rng(mu[site_y[n]], sigma);

  for (m in 1:M)
    z_rep[m] = normal_rng(epsilon[site_z[m]], tau);
}
