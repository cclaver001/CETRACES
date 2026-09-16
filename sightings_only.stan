data {
  int<lower=1> S;              // Number of cells in the grid
 
  // ───── Sightings data ─────
  int<lower=1> M;              
  vector[M] z;             // sightings
  // int<lower=1, upper=S> site_z[M];
  array[M] int<lower=1, upper=S> site_z;

}

parameters {
  vector[S] X;                 // Latent abundance per site
  
  // ───── Sightings parameters ─────
  real kappa;
  real<lower=0> omega;
  real<lower=0> tau;
}

transformed parameters {
  vector[S] epsilon;

  // sightings
  epsilon = kappa + omega * X;
}

model {
  // ───── Priors ─────
  X ~ normal(0, 1);
  //X ~ normal(0, 1);
  
  // Sightings priors
  kappa  ~ normal(0, 1);
  omega   ~ gamma(1,1);
  tau  ~ exponential(1);

  // ───── Likelihoods ─────

  // Sightings likelihood
  z ~ normal(epsilon[site_z], tau);
}

generated quantities {
  vector[S] ratio_cell;

  // Back-transformed sightings metric per cell
  for (s in 1:S) {
    ratio_cell[s] = exp(kappa + omega*X[s]) - 1;
  }
}
