
library(foreach)
library(tidyverse)
library(rstan)
library(MCMCvis)
library(ggrepel)

# Install with devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

source("functions_pe.R")

data_pop <- readRDS("data/data_pop_adpe.rds")
data_stan_null <- readRDS("data/data_stan_null_adpe.rds")


# Calculate long-term ice area average ------------------------------------

sites <- data_pop$site_list$site_id

ice <- readRDS("data/data_forced_finn_500km_adpe.rds") %>%
  filter(season > 1978) %>%
  filter(site_id %in% sites) %>%
  select(site_id, season, aice) %>%
  group_by(site_id) %>%
  summarise(aice_avg = mean(aice))

aice_std <- (ice$aice_avg - mean(ice$aice_avg))/sd(ice$aice_avg)
data_stan_null$ice <- cbind(aice_std, aice_std^2)

# Run Stan model ----------------------------------------------------------

options(mc.cores = parallel::detectCores())

res_null <- stan(file = 'scripts/null_pop_v7.stan', 
                 data = data_stan_null,
                 iter = 4000,
                 control = list(adapt_delta = 0.99, 
                                max_treedepth = 15))

saveRDS(res_null, "data/results_pop_adpe.rds")


# Extract demographic parameters ------------------------------------------

n_sites <- data_stan_null$n_sites
n_seasons <- data_stan_null$n_seasons
sites <- data_pop$site_list$site_id

params_N <- 
  foreach(i = 1:n_sites) %:% 
  foreach(h = 1:n_seasons, .combine = "c") %do%
  paste("lz", "[", i, ",", h, "]", sep = "")

N_all <- foreach(i = 1:n_sites) %do% {
  MCMCsummary(res_null, params = params_N[[i]], ISB = F, exact = T)[,1]
}

season_relative <- foreach(i = 1:n_sites) %do% {
  min_year <- data_pop$abundance_initial$season_relative[i]
  max_year <- max(filter(data_pop$abundance_nests, site_id == sites[i])$season_relative)
  
  c(min_year, max_year)
}

names(season_relative) <- data_pop$site_list$site_id

seasons <-  foreach(i = 1:n_sites) %do% {
  min_year <- data_pop$abundance_initial$season[i]
  max_year <- max(filter(data_pop$abundance_nests, site_id == sites[i])$season)
  
  c(min_year, max_year)
}

names(seasons) <- data_pop$site_list$site_id

# subset model results with years with data
N_all <- foreach(i = 1:n_sites) %do% {
  N_all[[i]][season_relative[[i]][1]:(season_relative[[i]][2])]
}

# Replace abundance estimates with NA for years with no data
for(i in 1:n_sites) {
  x <- filter(data_pop$abundance_nests, site_id == sites[i])
  N_all[[i]][-c(1,x$season_relative)] <- NA
}

names(N_all) <- data_pop$site_list$site_id


# Calculate WPE of sites --------------------------------------------------

wpe_adpe <- foreach(i = 1:length(sites)) %do% {
  tryCatch({test_wpe(N_all[[i]], n_random = 1000)}, error=function(e){NA})
}
names(wpe_adpe) <- sites

idx_np <- which(map_dbl(wpe_adpe, function(x) x[2]) > 7)
wpe_adpe2 <- map_dbl(wpe_adpe, function(x) x[1])[idx_np]
wpe_sig <- map_dbl(wpe_adpe, function(x) x[3])[idx_np]
wpe_sig <- as.factor(as.numeric(wpe_sig < 0.05))

mu_r <- MCMCsummary(res_null, params = "mu_site")[,1]
beta <- MCMCchains(res_null, params = c("beta", "beta0"))
ice_sim <- seq(min(aice_std), max(aice_std), length.out = 100)
beta_idx <- sample(1:nrow(beta), 100)

z <- foreach(i = 1:length(beta_idx), .combine = "rbind") %do% {
  beta[beta_idx[i],3] + beta[beta_idx[i],1]*ice_sim + beta[beta_idx[i],2]*ice_sim^2
}

z2 <- pivot_longer(as.data.frame(z), cols = everything()) %>%
  add_column(iter = rep(1:100, each = 100),
             x = rep(ice_sim*sd(ice$aice_avg) + mean(ice$aice_avg), 100))

idx_wpe <- order(wpe_adpe2, decreasing = T)
sites_wpe <- factor(sites[idx_np][idx_wpe], levels = (sites[idx_np][idx_wpe]))

g_base <- ggplot() +
  geom_col(aes(x = sites_wpe, y = wpe_adpe2[idx_wpe], fill = wpe_sig[idx_wpe])) +
  scale_fill_manual(name = NULL,
                    labels = c("0" = "p > 0.05",
                               "1" = "p < 0.05"),
                    values = c("0" = "darkred", 
                               "1" = "darkblue")) +
  labs(y = "WPE") +
  theme(legend.position = "bottom",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 8, angle = 90),
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

g_quad <- ggplot() +
  geom_line(mapping = aes(x = z2$x,
                          y = z2$value,
                          group = z2$iter), alpha = 0.2) +
  #geom_point(mapping = aes(x = ice$aice_avg[-idx_wpe], 
                           #y = mu_r[-idx_wpe]), size = 4, alpha = 0.9, col = "grey") +
  #geom_point(mapping = aes(x = ice$aice_avg[-idx_wpe], 
                           #y = mu_r[-idx_wpe]), size = 4, shape =  1, col = "black") +
  #geom_point(mapping = aes(x = ice$aice_avg[idx_wpe], 
                           #y = mu_r[idx_wpe], col = wpe_sig), size = 6, alpha = 0.8) +
  #geom_point(mapping = aes(x = ice$aice_avg[idx_wpe], 
                          # y = mu_r[idx_wpe]), size = 6, shape = 1, col = "black") +
  #scale_y_continuous(limits = c(-0.07, 0.05)) +
  #scale_color_manual(name = NULL,
                     #labels = c("0" = "p > 0.05",
                                #"1" = "p < 0.05"),
                     #values = c("0" = "darkred", 
                      #          "1" = "darkblue")) +
  labs(x = "Ice Concentration", y = "Average Growth") +
  theme(legend.position = "bottom",
        title = element_text(size = 12),
        legend.text = element_text(size = 14),
        axis.title.x = element_text(size = 16),
        axis.title.y = element_text(size = 16),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 16),
        panel.border = element_blank(),
        panel.grid.minor = element_blank())


(g_base | g_quad) + 
  plot_annotation(tag_levels = "a")
#ggsave("plot_adpe.jpeg", width = 12, height = 6, units = "in")  

#ggsave("fig_poster13.pdf", width = 8, height = 6, units = "in") 

