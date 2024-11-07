
library(MCMCvis)
library(foreach)
library(tidyverse)
library(truncnorm)
library(patchwork)

data_pop <- readRDS("data/data_pop_adpe.rds")

# Parameter chains 
res_pop <- read_rds("data/results_pop_adpe.rds")

param_chains <- MCMCchains(res_pop, 
                            params = c("beta0", "beta", 
                                       "sigma_site", "mu_sigma"))

idx <- sample(1:nrow(param_chains), 100)

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
      b0[h] + b1[h]*x[k,] + b2[h]*(x[k,]^2) + e_r[k,,h] + e_s[k,h]
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


# Forecasts for only observed sites ---------------------------------------

# Mean sea ice across 50 members
env_arr <- array(dim = c(nrow(env[[1]]), ncol(env[[1]]), 50))
for (i in 1:50) {
  env_arr[,,i] <- env[[i]]
}

env_mean <- apply(env_arr, c(1,2), mean)
rownames(env_mean) <- rownames(env[[1]])

env_mean_cp <- env_mean[which(rownames(env_mean) %in% 
                                c("CRZE", "PGEO", "PETE")),]

# Parameter means
params_mean <- MCMCsummary(res_pop, params = c("beta0", "beta", "sigma_r"))

b0_mean <- params_mean[1,1]
b1_mean <- params_mean[2,1]
b2_mean <- params_mean[3,1]

## process variance of CRZE, PGEO, and PETE
sigma_r_mean <- params_mean[c(17,34,35),1]

# Only parameter uncertainty
r_param <- foreach(k = 1:nrow(env_mean_cp)) %:% 
    foreach(h = 1:1000, .combine = "rbind") %do% {
      b0[h] + b1[h]*env_mean_cp[k,] + b2[h]*(env_mean_cp[k,]^2)
    }

# Only annual process variance
e_r_cp <- rnorm(1000*nrow(env_mean_cp)*ncol(env_mean_cp), 
                0, sigma_r_mean) %>%
  array(dim = c(nrow(env_mean_cp), ncol(env_mean_cp), 1000))

r_proc <- foreach(k = 1:nrow(env_mean_cp)) %:% 
  foreach(h = 1:1000, .combine = "rbind") %do% {
    b0_mean + b1_mean*env_mean_cp[k,] + b2_mean*(env_mean_cp[k,]^2) + 
      e_r_cp[k,,h]
  }


# Only site level process variance projected annuaaly
e_s_cp <- rnorm(1000*nrow(env_mean_cp), 0, mean(s_s)) %>%
  matrix(nrow = ncol(env_mean_cp), ncol = 1000)

r_proc_s <- foreach(k = 1:nrow(env_mean_cp)) %:% 
  foreach(h = 1:1000, .combine = "rbind") %do% {
    b0_mean + b1_mean*env_mean_cp[k,] + b2_mean*(env_mean_cp[k,]^2) + 
      e_s_cp[k,h]
  }

# Only climate variance
r_clim <- 
  foreach(k = 1:nrow(env_mean_cp)) %:% 
    foreach(i = 1:length(env), .combine = "rbind") %do% {
      x <- env[[i]]
      x <- x[which(rownames(env_mean) %in% c("CRZE", "PGEO", "PETE")),]
      b0_mean + b1_mean*x[k,] + b2_mean*(x[k,]^2)
    }

# Only climate variance for all colonies
r_clim_all <- 
  foreach(k = 1:nrow(env_mean)) %:% 
  foreach(i = 1:length(env), .combine = "rbind") %do% {
    x <- env[[i]]
    b0_mean + b1_mean*x[k,] + b2_mean*(x[k,]^2)
  }

# Save results in csv
folder <- "growth_projections_climate"
dir.create(folder)

sites <- unique(row.names(env[[1]]))
for (i in 1:length(sites)) {
  
  site_file <- paste(folder, "/", sites[i], ".csv", sep = "")
  write.csv(r_clim_all[[i]], site_file)
  
}


# Plot growth trajectories ------------------------------------------------

theme_set(theme_bw())

plot_traj <- function(r_est, site_name, uncertainty_name) {
  
  y_min <- apply(r_est, 2, quantile, 0.05, na.rm = T)
  y_max <- apply(r_est, 2, quantile, 0.95, na.rm = T)
  y_mean <- apply(r_est, 2, mean, na.rm = T)
  
  data_plot <- data.frame(y_min = y_min,
                          y_max = y_max,
                          y_mean = y_mean,
                          years = as.numeric(colnames(env[[1]])))
  
  ggplot(data_plot) +
    geom_line(aes(y = y_mean, x = years)) +
    geom_ribbon(aes(x = years, ymin = y_min, ymax = y_max), 
                fill = "blue4", alpha = 0.3) +
    labs(y = "Population Growth", 
         title = paste(site_name, uncertainty_name, sep = ", ")) +
    theme(panel.border = element_blank(),
          panel.grid.minor = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(size = 16),
          axis.title = element_text(size = 14),
          axis.text = element_text(size = 12)) 
}

g1 <- plot_traj(r_clim[[3]], "PGEO", "Climate Uncertainty")
g2 <- plot_traj(r_param[[3]], "PGEO", "Parameter Uncertainty")
g3 <- plot_traj(r_proc_s[[3]], "PGEO", "Spatial Process Variance")
g4 <- plot_traj(r_proc[[3]], "PGEO", "Temporal Process Variance")

(g1 | g2) / (g3 | g4)
ggsave("fig_r_traj.jpeg", width = 12, height = 10, units = "in")


