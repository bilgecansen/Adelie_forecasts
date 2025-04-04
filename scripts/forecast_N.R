library(MCMCvis)
library(foreach)
library(tidyverse)
library(patchwork)
#devtools::install_github('CCheCastaldo/mapppdr', build_vignettes = TRUE)
library(mapppdr)
library(sf)
library(ggspatial)
library(ggthemes)

sites_adpe <- mapppdr::site_species %>%
  filter(species_id == "ADPE")

sites_adpe_sf <- mapppdr::sites_sf %>%
  filter(site_id %in% sites_adpe$site_id)

data_pop <- readRDS("data/data_pop_adpe.rds")

data_site <- read_csv("data/colony_info_adelie.csv")

# Parameter means
res_pop <- read_rds("data/results_pop_adpe.rds")

params_mean <- MCMCsummary(res_pop, params = c("beta0", "beta"))

b0_mean <- params_mean[1, 1]
b1_mean <- params_mean[2, 1]
b2_mean <- params_mean[3, 1]

# Forecast and hindcast
env <- readRDS("data/data_coupled_normalized_raw_adpe.rds")
idx_s <- which(!rownames(env[[1]]) %in% sites_adpe$site_id)

## remove sites not in mapppd
for (i in 1:length(env)) {
  env[[i]] <- env[[i]][-idx_s, ]
}

# growth projections with only climate uncertainty
r_clim <-
  foreach(k = 1:nrow(env[[1]])) %:%
  foreach(i = 1:length(env), .combine = "rbind") %do%
  {
    x <- env[[i]]

    b0_mean + b1_mean * x[k, ] + b2_mean * (x[k, ]^2)
  }

r_clim_fore <- foreach(i = 1:length(r_clim)) %do%
  {
    idx_start <-
      which(colnames(r_clim[[i]]) == as.character(data_site$median_year[i]))

    r_clim[[i]][, idx_start:ncol(r_clim[[i]])]
  }

r_clim_hind <- foreach(i = 1:length(r_clim)) %do%
  {
    idx_start <-
      which(colnames(r_clim[[i]]) == as.character(data_site$median_year[i]))

    r_clim[[i]][, idx_start:1]
  }

# Abundance projections
N_start <- round(data_site$median)
N_start[is.na(N_start)] <- 1
N_start[N_start == 0] <- 1
N_start <- log(N_start)

## Forecast
N_fore_list <- list()
for (i in 1:length(r_clim_fore)) {
  N_mat <- matrix(
    ncol = ncol(r_clim_fore[[i]]) + 1,
    nrow = nrow(r_clim_fore[[i]])
  )

  N_mat[, 1] <- N_start[i]

  for (k in 1:nrow(r_clim[[1]])) {
    for (h in 1:ncol(r_clim_fore[[i]])) {
      N_mat[k, h + 1] <- N_mat[k, h] + r_clim_fore[[i]][k, h]
    }
  }

  N_mat[N_mat < 0] <- 0
  N_mat <- N_mat[, -ncol(N_mat)]

  colnames(N_mat) <- colnames(r_clim_fore[[i]])

  N_fore_list[[i]] <- N_mat
}

## Hindcast
N_hind_list <- list()
for (i in 1:length(r_clim_hind)) {
  N_mat <- matrix(
    ncol = ncol(r_clim_hind[[i]]) + 1,
    nrow = nrow(r_clim_hind[[i]])
  )

  N_mat[, 1] <- N_start[i]

  for (k in 1:nrow(r_clim[[1]])) {
    for (h in 1:ncol(r_clim_hind[[i]])) {
      N_mat[k, h + 1] <- N_mat[k, h] - r_clim_hind[[i]][k, h]
    }
  }

  N_mat[N_mat < 0] <- 0
  N_mat <- N_mat[, c(-1, -ncol(N_mat))]
  N_mat <- N_mat[, ncol(N_mat):1]

  colnames(N_mat) <- colnames(r_clim_hind[[i]][, (ncol(r_clim_hind[[i]])):2])

  N_hind_list[[i]] <- N_mat
}

N_list <- map2(N_hind_list, N_fore_list, function(x, y) cbind(x, y))

# Weighted average of regional growth
idx1 <- which(
  data_site$longitude <= -53 &
    data_site$longitude >= -64 &
    data_site$latitude <= -62 &
    data_site$latitude >= -67
)

idx4 <- which(
  data_site$longitude <= 142 &
    data_site$longitude >= 136 &
    data_site$latitude <= -65 &
    data_site$latitude >= -67
)

z <- data_site
z$longitude[z$longitude < 0] <- z$longitude[z$longitude < 0] + 360
idx5 <- which(
  z$longitude <= 200 &
    z$longitude >= 165 &
    z$latitude <= -70 &
    z$latitude >= -78
)

idx6 <- which(
  data_site$longitude <= -80 &
    data_site$longitude >= -120 &
    data_site$latitude <= -68 &
    data_site$latitude >= -75
)

idx7 <- which(
  data_site$longitude <= 85 &
    data_site$longitude >= 65 &
    data_site$latitude <= -65 &
    data_site$latitude >= -70
)

calc_weighted_r <- function(idx) {
  foreach(i = 1:50, .combine = "rbind") %do%
    {
      r <- r_clim[idx]
      N <- N_list[idx]

      z <- map_dfr(r, function(x) exp(x[i, ])) %>%
        as.matrix()

      w <- map_dfr(N, function(x) exp(x[i, ])) %>%
        as.matrix()

      log(colSums(z * w) / colSums(w))
    }
}

r1 <- calc_weighted_r(idx1)
r4 <- calc_weighted_r(idx4)
r5 <- calc_weighted_r(idx5)
r6 <- calc_weighted_r(idx6)
r7 <- calc_weighted_r(idx7)
r8 <- calc_weighted_r(1:287)

write.csv(r1, "r_region1_adpe.csv", row.names = F)
write.csv(r4, "r_region4_adpe.csv", row.names = F)
write.csv(r5, "r_region5_adpe.csv", row.names = F)
write.csv(r6, "r_region6_adpe.csv", row.names = F)
write.csv(r7, "r_region7_adpe.csv", row.names = F)
write.csv(r8, "r_cpolar_adpe.csv", row.names = F)

m1 <- data_site[idx1, ] %>%
  select(site_id, longitude, latitude) %>%
  add_column(region = "Region 1")
m4 <- data_site[idx4, ] %>%
  select(site_id, longitude, latitude) %>%
  add_column(region = "Region 4")
m5 <- data_site[idx5, ] %>%
  select(site_id, longitude, latitude) %>%
  add_column(region = "Region 5")
m6 <- data_site[idx6, ] %>%
  select(site_id, longitude, latitude) %>%
  add_column(region = "Region 6")
m7 <- data_site[idx7, ] %>%
  select(site_id, longitude, latitude) %>%
  add_column(region = "Region 7")

site_table <- rbind(m1, m4, m5, m6, m7)
write.csv(site_table, "site_table.csv", row.names = F)


# Plots -------------------------------------------------------------------

# Regional growth trajectories
plot_traj_r <- function(r_est, site_name) {
  theme_set(theme_bw())

  y_min <- apply(r_est, 2, quantile, 0.05, na.rm = T)
  y_max <- apply(r_est, 2, quantile, 0.95, na.rm = T)
  y_mean <- apply(r_est, 2, mean, na.rm = T)

  data_plot <- data.frame(
    y_min = y_min,
    y_max = y_max,
    y_mean = y_mean,
    years = as.numeric(colnames(env[[1]]))
  )

  ggplot(data_plot) +
    geom_ribbon(
      aes(x = years, ymin = y_min, ymax = y_max),
      fill = "#006BA4",
      alpha = 0.6
    ) +
    geom_line(aes(y = y_mean, x = years)) +
    labs(y = "Population Growth", x = "Year", title = site_name) +
    theme(
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      #axis.title.x = element_blank(),
      plot.title = element_text(size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "none"
    )
}

g1 <- plot_traj_r(r1, "Antarctic Peninsula Tip")
g2 <- plot_traj_r(r4, "Terre Adelie")
g3 <- plot_traj_r(r5, "Ross Sea")
g4 <- plot_traj_r(r6, "Amundsen-Bellingshausen")
g5 <- plot_traj_r(r7, "Prydz Bay")
g6 <- plot_traj_r(r8, "Circumpolar")

(g1 + g2 + g3) /
  (g4 + g5 + g6)

ggsave("fig_region.pdf", width = 180, height = 150, units = "mm", dpi = 600)

# Abundance trajectories
p_obs <- mapppdr::penguin_obs %>%
  filter(species_id == "ADPE", type != "chicks")
p_obs$count[p_obs$type == "nests"] <- p_obs$count[p_obs$type == "nests"] * 2

plot_traj_N <- function(N_est, site_name, xlim = NA, ylim) {
  theme_set(theme_bw())

  y_min <- apply(N_est, 2, quantile, 0.05, na.rm = T)
  y_max <- apply(N_est, 2, quantile, 0.95, na.rm = T)
  y_mean <- apply(N_est, 2, mean, na.rm = T)

  data_plot <- data.frame(
    y_min = y_min,
    y_max = y_max,
    y_mean = y_mean,
    years = as.numeric(colnames(env[[1]]))
  )

  g <- ggplot() +
    geom_ribbon(
      data = data_plot,
      mapping = aes(x = years, ymin = y_min, ymax = y_max),
      fill = "#006BA4",
      alpha = 0.6
    ) +
    geom_line(data = data_plot, mapping = aes(y = y_mean, x = years)) +
    geom_line(
      data = filter(p_obs, site_id == site_name),
      mapping = aes(x = season, y = count),
      col = "#FF800E",
      linewidth = 1.1
    ) +
    labs(y = "Colony Abundance", x = "Year", title = site_name) +
    theme(
      panel.border = element_blank(),
      panel.grid.minor = element_blank(),
      #axis.title.x = element_blank(),
      plot.title = element_text(size = 11),
      axis.title = element_text(size = 11),
      axis.text = element_text(size = 10),
      legend.position = "none"
    )

  if (is.na(xlim[1])) {
    g
  } else {
    g +
      scale_x_continuous(limits = xlim) +
      scale_y_continuous(limits = ylim)
  }
}

g6 <- plot_traj_N(
  exp(N_list[[68]]),
  "CRZE"
)

g7 <- plot_traj_N(
  exp(N_list[[151]]),
  "LLAN"
)

g8 <- plot_traj_N(
  exp(N_list[[3]]),
  "ADAR",
  xlim = c(1970, 2020),
  ylim = c(0, 1750000)
)

g9 <- plot_traj_N(
  exp(N_list[[200]]),
  "PGEO",
  xlim = c(1970, 2020),
  ylim = c(0, 150000)
)

(g6 + g7) /
  (g8 + g9)

ggsave("fig_N.pdf", width = 150, height = 150, units = "mm", dpi = 600)
