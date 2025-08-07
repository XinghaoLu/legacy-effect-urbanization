rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(fields)

## Retrieve the species list
spi.out.red <- read.table("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/Temporally_weigthed_regressions_sp.csv", header = T, sep = ",")
spi.out.red <- na.omit(spi.out.red, cols = "species")

## Model selection
delta.aic <- array(NA, dim = c(dim(spi.out.red)[1], 3))
colnames(delta.aic) <- c('aic', 'aic.null', 'aic.gauss')
models.gauss <- list.files("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/")

for(i in 1:dim(spi.out.red)[1]){ 
  print(i)
  name.sp <- spi.out.red$species[i]
  fil.gauss <- paste("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/",  
                     gsub(" ","_", name.sp), ".rds", sep="")
  if(file.exists(fil.gauss)) {
    mod.gauss <- readRDS(fil.gauss)
    delta.aic[i,1] <- mod.gauss$mod0$aic
    delta.aic[i,2] <- mod.gauss$mod$aic
  }
}

## Residual validation
res <- array(NA, dim = c(dim(spi.out.red)[1],5897))
for(i in 1:dim(spi.out.red)[1]){
  print(i)
  name.sp <- spi.out.red$species[i]
  models <- list.files("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/")
  fil <- paste("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/",  gsub(" ","_", name.sp, models), ".rds", sep="")
  mod <- readRDS(fil)
  res[i,] <- resid(mod$mod, type ="deviance") #
  if(t.test(res[i,])$p.value <= 0.05){print('invalid')} 
}

## Moran validation
China.grid.comp <- terra::vect("/2410009/data/02_timelag/inputs/China/China.comp.shp")
veg.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/veg1990-2020.csv", sep=",", h = T)
wat.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/wat1990-2020.csv", sep=",", h = T)
shdi.comp <- read.table("/2410009/data/02_timelag/Temporal_matrix/shdi1990-2020.csv", sep=",", h = T)
frac.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/frac1990-2020.csv", sep=",", h = T)
contag.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/contag1990-2020.csv", sep=",", h = T)
ntl.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/ntl1990-2020.csv", sep=",", h = T)
tmp.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/tmp1990-2020.csv", sep=",", h = T)
pop.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/pop1990-2020.csv", sep=",", h = T)

centro <- terra::geom(terra::centroids(China.grid.comp[ colSums(is.na(veg.comp))==0 ], inside = T))

## AUC validation
length(which(spi.out.red$AUC.mean < 0.7)) # Remove species with AUC < 0.7
spi.out.red <- filter(spi.out.red, AUC.mean >= 0.7)
mean(spi.out.red$AUC.mean)
sd(spi.out.red$AUC.mean)
mean(spi.out.red$sens.sp)
sd(spi.out.red$sens.sp)
mean(spi.out.red$spe.sp)
sd(spi.out.red$spe.sp)

### Comparison of the estimates of the equilibrium and non-equilibrium models

## Retrieve estimates
coeff0 <- array(NA, dim = c(dim(spi.out.red)[1], 9))
coeff <- array(NA, dim = c(dim(spi.out.red)[1], 9))
for(i in 1:dim(spi.out.red)[1]){
  print(i)
  name.sp <- spi.out.red$species[i]
  models <- list.files(("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/"))
  fil <- paste("/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/",  gsub(" ","_", name.sp, models), ".rds", sep="")
  mod <- readRDS(fil)
  coeff0[i,] <- mod$mod0$coefficients 
  coeff[i,] <- mod$mod$coefficients
}

## Compare the direction and magnitude of the effects
length(which( coeff[,2] > 0 & coeff0[,2] < 0)) # 0%
length(which( coeff[,2] < 0 & coeff0[,2] > 0)) # 0%
length(which( coeff[,3] > 0 & coeff0[,3] < 0)) # 26%
length(which( coeff[,3] < 0 & coeff0[,3] > 0)) # 0%
length(which( coeff[,4] > 0 & coeff0[,4] < 0)) # 0%
length(which( coeff[,4] < 0 & coeff0[,4] > 0)) # 5%
length(which( coeff[,5] > 0 & coeff0[,5] < 0)) # 6%
length(which( coeff[,5] < 0 & coeff0[,5] > 0)) 

t.test(coeff[,2], coeff0[,2], paired = T)
t.test(coeff[,2], coeff0[,2], paired = T, "greater")
t.test(coeff[,3], coeff0[,3], paired = T)
t.test(coeff[,3], coeff0[,3], paired = T, "greater")
t.test(coeff[,4], coeff0[,4], paired = T)
t.test(coeff[,5], coeff0[,5], paired = T)

### Plot the Gaussian weighting functions

################################################################################
par(mfrow=c(3,3))

q95 <- function(lag.param){which(( exp(-0.5*(seq(from=1,to=1000, by=0.1)/(lag.param))^2) )<= 0.05)[1]/10}

## Vegetation
gamma <- spi.out.red$range.veg.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "veg")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Water
gamma <- spi.out.red$range.wat.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "wat")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Shannon
gamma <- spi.out.red$range.shdi.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "shdi")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Fractal dimension
gamma <- spi.out.red$range.frac.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "frac")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Contagion
gamma <- spi.out.red$range.contag.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "contag")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Temperature
gamma <- spi.out.red$range.tmp.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "tmp")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Population
gamma <- spi.out.red$range.pop.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "pop")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

## Night-time light
gamma <- spi.out.red$range.ntl.comp
extent.gamma <- apply(t(t(gamma)),1,q95)
extent.gamma.it = c()

min_val <- min(gamma)/30
max_val <- max(gamma)/30
color_palette <- colorRampPalette(c("#f2f2f2","#f7efba","#dcdf80","#9ac150", "#54a964", "#2f8a79", "#325b7d"))(100)

normalized_gamma <- gamma/30
colors <- color_palette[ceiling(normalized_gamma * 99) + 1]

plot(NA, xlim=c(-30,0), ylim=c(0,1), 
     xlab=expression(paste(Delta,"t from present (years)")), 
     ylab="Estimate weighting", 
     main = "ntl")

for(i in 1:length(gamma)){ 
  a <- gamma[i]
  if(extent.gamma[i] < 1000){ 
    extent.gamma.it <- c(extent.gamma.it, a)
    curve(exp(-0.5*(x/a)^2), from = -30, to = 0, 
          col=colors[i], add=TRUE, lwd=1)
  }
}

### Bar charts with error bars
image.plot(legend.only = TRUE, zlim = c(min_val, max_val), 
           col = color_palette,
           smallplot = c(0.85,0.88,0.2,0.8),
           axis.args = list(cex.axis=0.8),
           legend.args = list(text="legend", side=4, line=2.5, cex=0.8))


range_data <- spi.out.red %>%
  select(starts_with("range.")) %>%
  mutate(across(everything(), ~ ./30))

summary_stats <- range_data %>%
  pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
  group_by(variable) %>%
  summarise(
    mean = mean(value),
    se = sd(value)/sqrt(n())
  ) %>%
  mutate(variable = gsub("range.", "", variable),
         variable = gsub(".comp", "", variable))

ggplot(summary_stats, aes(x = variable, y = mean)) +
  geom_bar(stat = "identity", fill = "#2f8a79", alpha = 0.8) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                width = 0.2, color = "#325b7d", linewidth = 0.8) +
  labs(title = "",
       x = "Environmental Variable",
       y = "") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))


