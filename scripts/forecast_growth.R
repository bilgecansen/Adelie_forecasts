library(MCMCvis)
library(foreach)
library(tidyverse)
library(truncnorm)
library(patchwork)
#devtools::install_github('CCheCastaldo/mapppdr', build_vignettes = TRUE)
library(mapppdr)
library(sf)
library(ggspatial)
library(ggthemes)

sites_adpe <- mapppdr::site_species %>% filter(species_id == "ADPE")
sites_adpe_sf <- mapppdr::sites_sf %>%
  filter(site_id %in% sites_adpe$site_id)

data_pop <- readRDS("data/data_pop_adpe.rds")

# Parameter chains
res_pop <- read_rds("data/results_pop_adpe.rds")

param_chains <- MCMCchains(
  res_pop,
  params = c("beta0", "beta", "sigma_site", "mu_sigma")
)

idx <- sample(1:nrow(param_chains), 100)

b0 <- param_chains[idx, "beta0"]
b1 <- param_chains[idx, "beta[1]"]
b2 <- param_chains[idx, "beta[2]"]

s_s <- param_chains[idx, "sigma_site"]
s_r <- rtruncnorm(100, a = 0, mean = 0, sd = param_chains[idx, "mu_sigma"])

# Forecast and hindcast
env <- readRDS("data/data_coupled_normalized_adpe.rds")
idx_s <- which(!rownames(env[[1]]) %in% sites_adpe$site_id)

## remove sites not in mapppd
for (i in 1:length(env)) {
  env[[i]] <- env[[i]][-idx_s, ]
}

e_s <- rnorm(100 * nrow(env[[1]]), 0, s_s) %>%
  matrix(nrow = nrow(env[[1]]), ncol = 100)

e_r <- rnorm(100 * nrow(env[[1]]) * ncol(env[[1]]), 0, s_r) %>%
  array(dim = c(nrow(env[[1]]), ncol(env[[1]]), 100))

r <- foreach(i = 1:length(env)) %do%
  {
    x <- env[[i]]

    foreach(k = 1:nrow(x)) %:%
      foreach(h = 1:100, .combine = "rbind") %do%
      {
        b0[h] + b1[h] * x[k, ] + b2[h] * (x[k, ]^2) + e_r[k, , h] + e_s[k, h]
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


# Projections with different sources of uncertainty -----------------------

# Mean sea ice across 50 members
env_arr <- array(dim = c(nrow(env[[1]]), ncol(env[[1]]), 50))
for (i in 1:50) {
  env_arr[,, i] <- env[[i]]
}

env_mean <- apply(env_arr, c(1, 2), mean)
rownames(env_mean) <- rownames(env[[1]])

# Parameter means
params_mean <- MCMCsummary(res_pop, params = c("beta0", "beta"))

b0_mean <- params_mean[1, 1]
b1_mean <- params_mean[2, 1]
b2_mean <- params_mean[3, 1]

# Only parameter uncertainty
r_param <- foreach(k = 1:nrow(env_mean)) %:%
  foreach(h = 1:100, .combine = "rbind") %do%
  {
    b0[h] + b1[h] * env_mean[k, ] + b2[h] * (env_mean[k, ]^2)
  }

for (i in 1:length(r_param)) {
  colnames(r_param[[i]]) <- colnames(env[[1]])
}

# Only process uncertainty
e_s2 <- rnorm(100 * nrow(env[[1]]), 0, mean(s_s)) %>%
  matrix(nrow = nrow(env[[1]]), ncol = 100)

r_proc_s <- foreach(k = 1:nrow(env_mean)) %:%
  foreach(h = 1:100, .combine = "rbind") %do%
  {
    b0_mean + b1_mean * env_mean[k, ] + b2_mean * (env_mean[k, ]^2) + e_s2[k, h]
  }

for (i in 1:length(r_proc_s)) {
  colnames(r_proc_s[[i]]) <- colnames(env[[1]])
}

# Only climate uncertainty
r_clim <-
  foreach(k = 1:nrow(env_mean)) %:%
  foreach(i = 1:length(env), .combine = "rbind") %do%
  {
    x <- env[[i]]

    b0_mean + b1_mean * x[k, ] + b2_mean * (x[k, ]^2)
  }

# Climate + site level + parameter
r_all <- foreach(i = 1:nrow(env[[1]])) %do%
  {
    x <- env[[i]]

    foreach(k = 1:nrow(x)) %:%
      foreach(h = 1:100, .combine = "rbind") %do%
      {
        b0[h] + b1[h] * x[k, ] + b2[h] * (x[k, ]^2) + e_s[k, h]
      }
  }
r_all <-
  foreach(k = 1:nrow(env[[1]])) %:%
  foreach(i = 1:length(env), .combine = "rbind") %do%
  {
    r_all[[i]][[k]]
  }


# Save results in csv
folder <- "growth_projections_climate"
dir.create(folder)

sites <- unique(row.names(env[[1]]))
for (i in 1:length(sites)) {
  site_file <- paste(folder, "/", sites[i], ".csv", sep = "")
  write.csv(r_clim[[i]], site_file)
}


# Determine TOE for each colony -------------------------------------------

get_toe <- function(r_est) {
  y_min <- apply(r_est, 2, quantile, 0.05, na.rm = T)
  y_max <- apply(r_est, 2, quantile, 0.95, na.rm = T)
  y_mean <- apply(r_est, 2, mean, na.rm = T)

  TOE_dec <- which(min(y_min[1:32]) > y_max[33:162])[1] %>%
    names() %>%
    as.numeric()

  TOE_inc <- which(max(y_max[1:32]) < y_min[33:162])[1] %>%
    names() %>%
    as.numeric()

  if (y_mean[1] > y_mean[162]) {
    TOE_dec
  } else TOE_inc
}

TOE_param <- map_dbl(r_param, get_toe)
TOE_proc <- map_dbl(r_proc_s, get_toe)
TOE_clim <- map_dbl(r_clim, get_toe)
TOE_all <- map_dbl(r_all, get_toe)

cat_toe <- function(toe) {
  if (is.na(toe)) return(">2100")
  if (toe <= 2025) return("<2025")
  if (toe > 2025 & toe <= 2050) return("2025 - 2050")
  if (toe > 2050 & toe <= 2075) return("2050 - 2075")
  if (toe > 2075 & toe <= 2100) return("2075 - 2100")
}

TOE_param2 <- factor(
  sapply(TOE_param, cat_toe),
  levels = c("<2025", "2025 - 2050", "2050 - 2075", "2075 - 2100", ">2100")
)
TOE_proc2 <- factor(
  sapply(TOE_proc, cat_toe),
  levels = c("<2025", "2025 - 2050", "2050 - 2075", "2075 - 2100", ">2100")
)
TOE_clim2 <- factor(
  sapply(TOE_clim, cat_toe),
  levels = c("<2025", "2025 - 2050", "2050 - 2075", "2075 - 2100", ">2100")
)
TOE_all2 <- factor(
  sapply(TOE_all, cat_toe),
  levels = c("<2025", "2025 - 2050", "2050 - 2075", "2075 - 2100", ">2100")
)


# Plot growth trajectories ------------------------------------------------

plot_traj <- function(r_est, site_name, uncertainty_name, TOE) {
  theme_set(theme_bw())

  y_min <- apply(r_est, 2, quantile, 0.05, na.rm = T)
  y_max <- apply(r_est, 2, quantile, 0.95, na.rm = T)
  y_mean <- apply(r_est, 2, mean, na.rm = T)

  if (y_mean[1] > y_mean[162]) {
    y_thr <- y_min
  } else y_thr <- y_max

  data_plot <- data.frame(
    y_min = y_min,
    y_max = y_max,
    y_mean = y_mean,
    years = as.numeric(colnames(env[[1]])),
    hist = as.factor(
      c(rep("baseline", 32), rep("proj", 130))
    )
  )

  ggplot(data_plot) +
    geom_ribbon(
      aes(x = years, ymin = y_min, ymax = y_max, fill = hist),
      alpha = 0.6
    ) +
    geom_line(aes(y = y_mean, x = years)) +
    labs(
      y = "Population Growth",
      x = "Year",
      title = paste(site_name, uncertainty_name, sep = ", ")
    ) +
    geom_abline(intercept = min(y_thr[1:32]), slope = 0, linetype = 2) +
    geom_vline(xintercept = TOE, linetype = 2) +
    scale_fill_manual(values = c("#006BA4", "#FF800E")) +
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

g1 <- plot_traj(
  r_clim[[172]],
  "MIZU",
  "Climate Uncertainty",
  TOE = TOE_clim[172]
)

g2 <- plot_traj(
  r_proc_s[[172]],
  "MIZU",
  "Spatial Process Uncertainty",
  TOE = TOE_proc[172]
)

g3 <- plot_traj(
  r_param[[172]],
  "MIZU",
  "Parameter Uncertainty",
  TOE = TOE_param[172]
)

g4 <- plot_traj(r_all[[172]], "MIZU", "All Uncertainties", TOE = TOE_all[172])

(((g1 + theme(axis.title.x = element_blank())) /
  (g3 + theme(axis.title = element_blank()))) |
  (g2 + theme(axis.title = element_blank())) /
    (g4 + theme(axis.title.y = element_blank()))) +
  plot_annotation(tag_levels = list(c("a", "c", "b", "d")))

ggsave("fig2.pdf", width = 180, height = 180, units = "mm", dpi = 600)

g5 <- plot_traj(r_clim[[68]], "CRZE", "Climate Uncertainty", TOE = TOE_clim[68])

g6 <- plot_traj(
  r_proc_s[[68]],
  "CRZE",
  "Spatial Process Uncertainty",
  TOE = TOE_proc[68]
)

g7 <- plot_traj(
  r_param[[68]],
  "CRZE",
  "Parameter Uncertainty",
  TOE = TOE_param[68]
)

g8 <- plot_traj(r_all[[68]], "CRZE", "All Uncertainties", TOE = TOE_all[68])

(((g5 + theme(axis.title.x = element_blank())) /
  (g7 + theme(axis.title = element_blank()))) |
  (g6 + theme(axis.title = element_blank())) /
    (g8 + theme(axis.title.y = element_blank()))) +
  plot_annotation(tag_levels = list(c("a", "c", "b", "d")))

ggsave("fig3.pdf", width = 180, height = 180, units = "mm", dpi = 600)

g3
ggsave("fig_traj.pdf", width = 120, height = 90, units = "mm", dpi = 600)


# Plot maps ---------------------------------------------------------------

antarctic <- st_as_sf(maps::map('world', plot = FALSE, fill = TRUE))
idx <- which(antarctic$ID == "Antarctica")
antarctic <- antarctic[idx, ]
antarctic <- st_transform(antarctic, 3031)

sites_adpe_sf <- st_transform(sites_adpe_sf, 3031)

m1 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5) +
  annotation_scale() +
  geom_sf(
    data = sites_adpe_sf,
    aes(col = TOE_clim2),
    alpha = 0.5,
    show.legend = T
  ) +
  scale_color_manual(
    values = c(
      "<2025" = "#009E73",
      "2025 - 2050" = "#56B4E9",
      "2050 - 2075" = "#F0E442",
      "2075 - 2100" = "#CC79A7", #"#E69F00",
      ">2100" = "#D55E00"
    ),
    drop = FALSE
  ) +
  labs(col = "TOE", title = "Climate Uncertainty") +
  coord_sf(label_axes = list(bottom = "E", left = "E")) +
  theme(legend.position = "bottom", plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

m2 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5) +
  annotation_scale() +
  geom_sf(
    data = sites_adpe_sf,
    aes(col = TOE_proc2),
    alpha = 0.5,
    show.legend = T
  ) +
  scale_color_manual(
    values = c(
      "<2025" = "#009E73",
      "2025 - 2050" = "#56B4E9",
      "2050 - 2075" = "#F0E442",
      "2075 - 2100" = "#CC79A7", #"#E69F00",
      ">2100" = "#D55E00"
    ),
    drop = FALSE
  ) +
  labs(col = "TOE", title = "Spatial Process Uncertainty") +
  coord_sf(label_axes = list(bottom = "E", left = "E")) +
  theme(legend.position = "bottom", plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

m3 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5) +
  annotation_scale() +
  geom_sf(
    data = sites_adpe_sf,
    aes(col = TOE_param2),
    alpha = 0.5,
    show.legend = T
  ) +
  scale_color_manual(
    values = c(
      "<2025" = "#009E73",
      "2025 - 2050" = "#56B4E9",
      "2050 - 2075" = "#F0E442",
      "2075 - 2100" = "#CC79A7", #"#E69F00",
      ">2100" = "#D55E00"
    ),
    drop = FALSE
  ) +
  labs(col = "TOE", title = "Parameter Uncertainty") +
  coord_sf(label_axes = list(bottom = "E", left = "E")) +
  theme(legend.position = "bottom", plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

m4 <- ggplot() +
  geom_sf(data = antarctic, alpha = 0.5) +
  annotation_scale() +
  geom_sf(
    data = sites_adpe_sf,
    aes(col = TOE_all2),
    alpha = 0.5,
    show.legend = T
  ) +
  scale_color_manual(
    values = c(
      "<2025" = "#009E73",
      "2025 - 2050" = "#56B4E9",
      "2050 - 2075" = "#F0E442",
      "2075 - 2100" = "#CC79A7", #"#E69F00",
      ">2100" = "#D55E00"
    ),
    drop = FALSE
  ) +
  labs(col = "TOE", title = "All Uncertainties") +
  coord_sf(label_axes = list(bottom = "E", left = "E")) +
  theme(legend.position = "bottom", plot.title = element_text(size = 11)) +
  guides(colour = guide_legend(override.aes = list(size = 3)))

(m1 | m2) /
  (m3 | m4) /
  guide_area() +
  plot_annotation(tag_levels = "a") +
  plot_layout(guides = "collect", heights = c(5, 5, 1))
ggsave("fig4.pdf", width = 180, height = 180, units = "mm", dpi = 600)


# Plot sea ice trajectories -----------------------------------------------

ice <- env <- readRDS("data/data_coupled_transformed_adpe.rds")
ice_mizu <- map_dfr(ice, function(x) x[which(rownames(x) == "MIZU"), ]) %>%
  pivot_longer(everything(), names_to = "year", values_to = "ice") %>%
  add_column(member = rep(1:50, each = 162)) %>%
  mutate(year = as.numeric(year))

ggplot(data = ice_mizu) +
  geom_line(aes(x = year, y = ice, group = member, col = member)) +
  scale_colour_gradient_tableau('Blue') +
  scale_x_continuous(breaks = c(seq(1900, 2100, 50))) +
  labs(x = "Year", y = "Sea Ice Concentration") +
  theme(
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    #axis.title.x = element_blank(),
    plot.title = element_text(size = 11),
    axis.title = element_text(size = 11),
    axis.text = element_text(size = 10),
    legend.position = "none"
  )

ggsave("fig_ice.pdf", width = 120, height = 90, units = "mm", dpi = 600)
