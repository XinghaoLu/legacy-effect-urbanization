rm(list = ls())

library(funspace)
library(RColorBrewer)
library(dplyr)

dat <- read.csv("2410009/data/04_species_lag/species_lag_res.csv", header =T, sep = ',', row.name = 1)
subset_dat <- dat %>%
  select(Hand.Wing.Index, Ave.generation.length, Clutch.size, urtolerance.1, Beak.Length_Culmen, Tail.Length)
subset_dat <- scale(subset_dat)
pca <- princomp(subset_dat)

funtest <- funspace(pca,PCs = c(1, 2),
                    n_divisions = 100)

color_set <- colorRampPalette(c("#FFF9BD","#D3D83D","#6EBC5C", "#279C7F","#31AA77", "#2D7086","#3B3F7C" ,"#382348"))(10)

plot(funtest, pnt = T, quant.plot = TRUE,n_divisions = 1000, pnt.cex = 0.2, arrows = TRUE, 
     xlim = c(-3.6,8.2),
     ylim = c(-3.5,5.2),col = color_set)