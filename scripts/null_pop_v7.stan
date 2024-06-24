
data {
  int<lower=1> n_sites;
  int<lower=1> n_seasons;
  int<lower=1> nests;
  int<lower=1> site_n[nests];
  int<lower=1> season_n[nests];
  int<lower=1> s[n_sites];
  int reg[n_sites];
  vector[nests] y_n;
  vector[n_sites] y_i;
  matrix[n_sites, 2] ice;
  vector[nests] sigma_n;
  vector[n_sites] sigma_i;
}

parameters {
  real mu_site_raw[n_sites];
  real<lower=0> mu_sigma;
  real<lower=0> tau;
  real<lower=0> sigma_r_raw[n_sites];
  matrix[n_sites, n_seasons] lz;
  real beta0;
  vector[2] beta;
  real<lower=0> sigma_site;
}

transformed parameters {
  real mu_site[n_sites];
  real<lower=0> sigma_r[n_sites];
  
  for (i in 1:n_sites) {
    mu_site[i] = (mu_site_raw[i]*sigma_site) + ice[i,]*beta + beta0;
    sigma_r[i] = (sigma_r_raw[i]*tau) + mu_sigma;
  }
}

model {
  // Priors
  for (i in 1:n_sites) {
    lz[i,s[i]] ~ normal(0, 100);
  }
  
  sigma_site ~ normal(0, 1);
  beta0 ~ normal(0, 1);
  beta ~ normal(0, 1);
  mu_site_raw ~ normal(0,1);
  sigma_r_raw ~ normal(0,1);
  mu_sigma ~ normal(0,1);
  tau ~ normal(0,1);
   
  // Annual population growth
  for (i in 1:n_sites) {
    for (t in s[i]:(n_seasons-1)) {
      lz[i,t+1] ~ normal(lz[i, t] + mu_site[i], sigma_r[i]);
    }
  }
  
  for (i in 1:n_sites) {
    for (t in 1:(s[i] - 1)) {
      lz[i, s[i] - t] ~ normal(lz[i, s[i] - t + 1] - mu_site[i], sigma_r[i]);
    }
  }
  
  // Nest counts at initial year 
  for (i in 1:n_sites) {
    y_i[i] ~ normal(lz[i, s[i]], sigma_n[i]);
  }
    
  // Nest counts after initial year
  for (i in 1:nests) {
    y_n[i] ~ normal(lz[site_n[i], season_n[i]], sigma_n[i]);
  }
}

generated quantities {
  matrix[n_sites, n_seasons-1] r;
  matrix[n_sites, n_seasons-1] r_std;
  real mu_site_pred[n_sites];
  
  for (i in 1:n_sites) {
    for (t in 1:(n_seasons-1)) {
      r[i,t] = lz[i,t+1] - lz[i,t];
      r_std[i,t] = r[i,t] - mu_site[i];
    }
  }
  
  for (i in 1:n_sites) {
    mu_site_pred[i] = ice[i,]*beta + beta0;
  }
}
