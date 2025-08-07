rm(list = ls())

set.seed(100)
library(tidyverse)
library(ape)
library(brms)
library(MCMCvis)
library(bayesplot)
library(ggplot2)
library(ggdist)

dat <- read.csv("2410009/data/04_species_lag/species_lag_res.csv", sep = ',', header = TRUE)
dat$species_tree <- gsub(" ", "_", dat$species_tree)
phy.tree <- read.tree("/2410009/data/03_comm_legacy/3.1_phylogenetic_legacy_res/tree.newick")
dat_tree <- dat[match(phy.tree$tip.label, dat$species_tree),]
dat_tree$frac <- dat_tree$γfrac

C <- vcv.phylo(phy.tree, corr = T) 

## Vegetation
mod.veg <- brm(veg ~ 
                 scale(Beak.Length_Culmen) +
                 scale(Tail.Length) +
                 scale(urtolerance.1) +
                 scale(Clutch.size) +
                 scale(Hand.Wing.Index) +
                 scale(Ave.generation.length) +
                 (1|a|gr(species_tree, cov = C)),
               family = Beta(),
               data   = dat_tree,
               data2  = list(C = C),
               control = list(adapt_delta = 0.99),
               chain = 3,
               iter = 3000)


save(mod.veg, file = "/2410009/data/05_life-history_mediation/bayes_veg.rda")

posterior.veg <- as.matrix(mod.veg$fit)
str(posterior.veg )

MCMCplot(posterior.veg, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_veg <- as.data.frame(posterior.veg) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.veg != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_veg, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Water
mod.wat <- brm(wat ~ 
                 scale(Beak.Length_Culmen) +
                 scale(Tail.Length) +
                 scale(urtolerance.1) +
                 scale(Clutch.size) +
                 scale(Hand.Wing.Index) +
                 scale(Ave.generation.length) +
                 (1|a|gr(species_tree, cov = C)),
               family = Beta(),
               data   = dat_tree,
               data2  = list(C = C),
               control = list(adapt_delta = 0.99),
               chain = 3,
               iter = 3000)


save(mod.wat, file = "/2410009/data/05_life-history_mediation/bayes_wat.rda")

posterior.wat <- as.matrix(mod.wat$fit)
str(posterior.wat )

MCMCplot(posterior.wat, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_wat <- as.data.frame(posterior.wat) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.wat != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_wat, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Shannon
mod.shdi <- brm(shdi ~ 
                  scale(Beak.Length_Culmen) +
                  scale(Tail.Length) +
                  scale(urtolerance.1) +
                  scale(Clutch.size) +
                  scale(Hand.Wing.Index) +
                  scale(Ave.generation.length) +
                  (1|a|gr(species_tree, cov = C)),
                family = Beta(),
                data   = dat_tree,
                data2  = list(C = C),
                control = list(adapt_delta = 0.99),
                chain = 3,
                iter = 3000)


save(mod.shdi, file = "/2410009/data/05_life-history_mediation/bayes_shdi.rda")

posterior.shdi <- as.matrix(mod.shdi$fit)
str(posterior.shdi )

MCMCplot(posterior.shdi, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_shdi <- as.data.frame(posterior.shdi) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.shdi != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_shdi, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Contagion
mod.contag <- brm(contag ~ 
                    scale(Beak.Length_Culmen) +
                    scale(Tail.Length) +
                    scale(urtolerance.1) +
                    scale(Clutch.size) +
                    scale(Hand.Wing.Index) +
                    scale(Ave.generation.length) +
                    (1|a|gr(species_tree, cov = C)),
                  family = Beta(),
                  data   = dat_tree,
                  data2  = list(C = C),
                  control = list(adapt_delta = 0.99),
                  chain = 3,
                  iter = 3000)


save(mod.contag, file = "/2410009/data/05_life-history_mediation/bayes_contag.rda")

posterior.contag <- as.matrix(mod.contag$fit)
str(posterior.contag )

MCMCplot(posterior.contag, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_contag <- as.data.frame(posterior.contag) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.contag != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_contag, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Fractal dimension
mod.frac <- brm(frac ~ 
             scale(Beak.Length_Culmen) +
             scale(Tail.Length) +
             scale(urtolerance.1) +
             scale(Clutch.size) +
             scale(Hand.Wing.Index) +
             scale(Ave.generation.length) +
             (1|a|gr(species_tree, cov = C)),
           family = Beta(),
           data   = dat_tree,
           data2  = list(C = C),
           control = list(adapt_delta = 0.99),
           chain = 3,
           iter = 3000)


save(mod.frac, file = "/2410009/data/05_life-history_mediation/bayes_frac.rda")

posterior.frac <- as.matrix(mod.frac$fit)
str(posterior.frac )

MCMCplot(posterior.frac, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_frac <- as.data.frame(posterior.frac) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.frac != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_frac, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Temperature
mod.tmp <- brm(tmp ~ 
                 scale(Beak.Length_Culmen) +
                 scale(Tail.Length) +
                 scale(urtolerance.1) +
                 scale(Clutch.size) +
                 scale(Hand.Wing.Index) +
                 scale(Ave.generation.length) +
                 (1|a|gr(species_tree, cov = C)),
               family = Beta(),
               data   = dat_tree,
               data2  = list(C = C),
               control = list(adapt_delta = 0.99),
               chain = 3,
               iter = 3000)


save(mod.tmp, file = "/2410009/data/05_life-history_mediation/bayes_tmp.rda")

posterior.tmp <- as.matrix(mod.tmp$fit)
str(posterior.tmp )

MCMCplot(posterior.tmp, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_tmp <- as.data.frame(posterior.tmp) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.tmp != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_tmp, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Population
mod.pop <- brm(pop ~ 
                 scale(Beak.Length_Culmen) +
                 scale(Tail.Length) +
                 scale(urtolerance.1) +
                 scale(Clutch.size) +
                 scale(Hand.Wing.Index) +
                 scale(Ave.generation.length) +
                 (1|a|gr(species_tree, cov = C)),
               family = Beta(),
               data   = dat_tree,
               data2  = list(C = C),
               control = list(adapt_delta = 0.99),
               chain = 3,
               iter = 3000)


save(mod.pop, file = "/2410009/data/05_life-history_mediation/bayes_pop.rda")

posterior.pop <- as.matrix(mod.pop$fit)
str(posterior.pop )

MCMCplot(posterior.pop, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_pop <- as.data.frame(posterior.pop) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.pop != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_pop, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )

## Night-time light
mod.ntl <- brm(ntl ~ 
                 scale(Beak.Length_Culmen) +
                 scale(Tail.Length) +
                 scale(urtolerance.1) +
                 scale(Clutch.size) +
                 scale(Hand.Wing.Index) +
                 scale(Ave.generation.length) +
                 (1|a|gr(species_tree, cov = C)),
               family = Beta(),
               data   = dat_tree,
               data2  = list(C = C),
               control = list(adapt_delta = 0.99),
               chain = 3,
               iter = 3000)


save(mod.ntl, file = "/2410009/data/05_life-history_mediation/bayes_ntl.rda")

posterior.ntl <- as.matrix(mod.ntl$fit)
str(posterior.ntl )

MCMCplot(posterior.ntl, ci = c(50, 95),ref_ovl = TRUE,
         params = c("b_scaleurtolerance.1", "b_scaleBeak.Length_Culmen", "b_scaleTail.Length", "b_scaleClutch.size", "b_scaleHand.Wing.Index", "b_scaleAve.generation.length"),
         exact = FALSE, main = 'var')

posterior_ntl <- as.data.frame(posterior.ntl) %>% 
  pivot_longer(
    cols = everything(),
    names_to = "parameter",
    values_to = "value"
  ) %>%
  filter(posterior.ntl != "b_Intercept") %>%
  filter(grepl("^b_", parameter))  

ggplot(posterior_ntl, aes(x = value, y = parameter)) +
  stat_halfeye(
    point_interval = mean_qi, 
    .width = c(0.80, 0.95),
    fill = "#1f77b4",
    slab_alpha = 0.7,
    point_size = 1.5
  ) +
  geom_vline(xintercept=0, linetype ="dashed") +
  labs(
    x = "Parameter Value", 
    y = "Parameter",
    title = "Posterior Distributions (Excluding Intercept)"
  ) +
  theme_classic() +
  theme(
    panel.grid.minor = element_blank(),
    axis.text.y = element_text(hjust = 0),
    plot.title = element_text(face = "bold")
  )