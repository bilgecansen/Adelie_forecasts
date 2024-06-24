
library(MCMCvis)
library(foreach)
library(tidyverse)
library(truncnorm)

# Parameter chains 
res_pop <- read_rds("data/results_pop_adpe.rds")

param_chains <- MCMCchains(res_pop, 
                            params = c("beta0", "beta", 
                                       "sigma_site", "mu_sigma"))

idx <- sample(1:nrow(params_chains), 100)

b0 <- param_chains[idx,"beta0"]
b1 <- param_chains[idx,"beta[1]"]
b2 <- param_chains[idx,"beta[2]"]

s_s <- param_chains[idx,"sigma_site"]
s_r <- rtruncnorm(100, a = 0, mean = 0, sd = param_chains[idx,"mu_sigma"])

# Forecast and hindcast
env <- readRDS("data/data_coupled_normalized_adpe.rds")

e_s <- rnorm(100*nrow(env[[1]]), 0, s_s) %>%
  matrix(nrow = nrow(env[[1]]), ncol = 100)

e_r <- rnorm(100*nrow(env[[1]])*ncol(env[[1]]), 0, s_r) %>%
  array(dim = c(nrow(env[[1]]), ncol(env[[1]]), 100))

r <- foreach(i = 1:length(env)) %do% {
  
  x <- env[[i]]
  
  foreach(k = 1:nrow(x)) %:% 
    foreach(h = 1:100, .combine = "rbind") %do% {
      b0[h] + b1[h]*x[k,] + b2[h]*(x[k,]^2) + e_s[k,h] + e_r[k,,h]
    }
}

# Save results in csv
folder <- "growth_projections"
dir.create(folder)

sites <- unique(row.names(env[[1]]))
for (i in 1:length(sites)) {
  
  site_folder <- paste(folder, sites[i], sep = "/")
  dir.create(site_folder)
  
  for (k in 1:50) {
    
    site_file <- paste(site_folder, "/", k, ".csv", sep = "")
    write.csv(r[[k]][[i]], site_file)
  }
}
