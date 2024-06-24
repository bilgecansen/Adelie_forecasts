
library(foreach)
library(tidyverse)
library(abind)

#empe_sites <- read.csv("empe_sites_new.csv")
data_pop <- readRDS("data/data_pop_adpe.rds")
data_stan_null <- readRDS("data/data_stan_null_adpe.rds")


# Load environmental data -------------------------------------------------

# Observed data
data_aice <- readRDS("data/data_forced_finn_500km_adpe.rds") %>%
  filter(season > 1978) %>%
  #filter(site_id %in% sites) %>%
  select(site_id, season, aice) #%>%
  #group_by(site_id) %>%
  #summarise(aice_avg = mean(aice))

# Forecast/hindcast data
data_aice_coupled <- readRDS("data/data_aice_coupled_adpe.rds") %>%
  map(., function(x) arrange(x, site_id))

sites_coupled <- unique(data_aice_coupled[[1]]$site_id)

get_env_mat <- function(data_env, x) {
  select(data_env, site_id, season, x) %>%
    pivot_wider(names_from = season, values_from = x) %>%
    select(-site_id) %>%
    as.matrix()
}

env_mat <- get_env_mat(data_aice, "aice")
env_mat_avg <- apply(env_mat, 1, mean)
x_mean <- mean(env_mat_avg)
x_sd <- sd(env_mat_avg)

env_mat_fore <- foreach(h = 1:length(data_aice_coupled)) %do% 
  get_env_mat(data_aice_coupled[[h]], "aice")

# correct env forecast data
trans_dat <- function(x, y) {
  
  mu_obs <- mean(x)
  sigma_obs <- sd(x)
  sigma_model <- sd(y[80:119])
  mu_model <- mean(y[80:119])
  
  (sigma_obs / sigma_model) * (y - mu_model) + mu_obs
}

env_mat_fore <- foreach(i = 1:length(env_mat_fore)) %:%
  foreach(k = 1:nrow(env_mat_fore[[i]]), .combine = "rbind") %do% {
    trans_dat(env_mat[k,], env_mat_fore[[i]][k,])
  }

for (i in 1:50) {
  env_mat_fore[[i]][env_mat_fore[[i]] < 0] <- 0
  env_mat_fore[[i]][env_mat_fore[[i]] > 1] <- 1
}   

# moving window average
env_mat_fore_avg <- foreach (i = 1:50) %:% 
  foreach(k = 1:length(sites_coupled), .combine = "rbind") %:%
  foreach (h = 1:162, .combine = "c") %do% {
    mean(env_mat_fore[[i]][k,h:(h+39)])
  }

# make sure forecasts are within data bounds
for (i in 1:50) {
  env_mat_fore_avg[[i]][env_mat_fore_avg[[i]] < min(env_mat_avg)] <- 
    min(env_mat_avg)
  env_mat_fore_avg[[i]][env_mat_fore_avg[[i]] > max(env_mat_avg)] <- 
    max(env_mat_avg)
} 

for (i in 1:50) {
  colnames(env_mat_fore_avg[[i]]) <- 1939:2100
  rownames(env_mat_fore_avg[[i]]) <- sites_coupled
}

saveRDS(env_mat_fore_avg, "data/data_coupled_transformed_adpe.rds")

# standardize the data
env_mat_fore_avg_std <- foreach (i = 1:50) %do% {
  (env_mat_fore_avg[[i]] - x_mean)/x_sd
}

saveRDS(env_mat_fore_avg_std, "data/data_coupled_normalized_adpe.rds")
