rm(list = ls())

library(cluster)
library(dplyr)
library(ggplot2)
library(factoextra)
library(vegan)

dat <- read.csv("2410009/data/01_trend_cluster/sen_slope.csv", header = T)
str(dat)

cluster_data <- dat %>% select(-gridID)
scaled_data <- scale(cluster_data)

set.seed(123)
pam_result <- pam(scaled_data, k = 3)

dat_with_cluster <- dat %>%
  mutate(cluster = pam_result$clustering)

write.csv(dat_with_cluster, "2410009/data/01_trend_cluster/clustered_res.csv", row.names = FALSE)