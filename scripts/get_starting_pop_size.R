library(mapppdr)
library(tidyverse)

sites <- read.csv("data/colonies_adelie.csv")

z <- mapppdr::penguin_obs %>%
  filter(species_id == "ADPE")

idx <- which(z$type == "nests")
z2 <- z
z2$count[idx] <- z2$count[idx] * 2

idx2 <- which(z$type == "chicks")
z2$count[idx2] <- NA

z3 <- z2 %>%
  group_by(site_id) %>%
  summarise(
    median = round(median(count, na.rm = T), 2),
    mean = round(mean(count, na.rm = T), 2),
    min = round(min(count, na.rm = T), 2),
    max = round(max(count, na.rm = T), 2),
    sd = round(sd(count, na.rm = T), 2),
    last_year = max(year),
    median_year = round(median(year))
  ) %>%
  as.data.frame()

z4 <- left_join(sites, z3, by = "site_id")

write.csv(z4, "data/colony_info_adelie.csv", row.names = F)
