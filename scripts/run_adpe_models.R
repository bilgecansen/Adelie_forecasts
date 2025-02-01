
library(foreach)
library(tidyverse)
library(rstan)
library(MCMCvis)
library(ggrepel)

# Install with devtools::install_github("hrbrmstr/ggalt")
library(ggalt)

#source("functions_pe.R")

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


# Plot ice-growth relationship --------------------------------------------

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

ggplot() +
  geom_line(mapping = aes(x = z2$x,
                          y = z2$value,
                          group = z2$iter), alpha = 0.1) +
  geom_line(aes(x = ice_sim*sd(ice$aice_avg) + mean(ice$aice_avg),
                y = z_mean), size = 2, linetype = 2) +
  geom_point(mapping = aes(x = ice$aice_avg, y = mu_r), 
             size = 3, alpha = 0.9, col = "#FF800E") +
  geom_point(mapping = aes(x = ice$aice_avg, y = mu_r), 
             size = 3, shape =  1, col = "black") +
  labs(x = "Ice Concentration", y = "40-Year Average Growth") +
  theme(axis.title = element_text(size = 11),
        axis.text = element_text(size = 10),
        panel.border = element_blank(),
        panel.grid.minor = element_blank())

ggsave("fig_quad.pdf", width = 120, height = 90, units = "mm", dpi = 600)

