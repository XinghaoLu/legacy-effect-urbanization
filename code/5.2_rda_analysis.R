rm(list = ls())

library(readxl)
library(RColorBrewer)
library(vegan)
library(ks)
library(rdacca.hp)

add.alpha <- function(col, alpha = 1) {
  if(missing(col)) stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) rgb(x[1], x[2], x[3], alpha=alpha))  
}

create_color_set <- function(df, group_var, palette) {
  group_values <- unique(df[[group_var]])
  colors <- brewer.pal(length(group_values), palette)
  names(colors) <- group_values
  return(colors)
}

dat <- read.csv("2410009/data/04_species_lag/species_lag_res.csv", sep = ',', header = TRUE)

x  <- dat[, c("urtolerance.1", "Hand.Wing.Index", "Ave.generation.length", "Clutch.size", "Beak.Length_Culmen", "Tail.Length")]    # 自变量
x <- scale(x)
y1 <- dat[, c("γveg", "γwat", "γshdi", "γfrac", "γcontag", "γntl", "γpop", "γtmp")]

xy_df <- na.omit(data.frame(y1, x))

y1 <- xy_df[, 1:ncol(y1)]
x  <- xy_df[, (ncol(y1)+1) : ncol(xy_df)]

mod <- rda(y1 ~ ., data = x, scale = TRUE)
mod

summary(mod)

mod.hp <- rdacca.hp(y1, x, method = "RDA", type = 'R2', scale = FALSE)
mod.hp

x.fit <- envfit(mod,x,permutations = 999)
r <- as.matrix(x.fit$vectors$r)
p <- as.matrix(x.fit$vectors$pvals)
x.p <- cbind(r,p)
colnames(x.p) <- c("r2", "p-value")
x.imp <- as.data.frame(x.p)
x.imp$x.imp <- x.imp$r2 / sum(x.imp$r2)
x.imp
x.imp$variable <- rownames(x.imp)

ggplot(x.imp, aes(x = reorder(variable, x.imp), y = x.imp, fill = variable)) +
  geom_col(width = 0.7, alpha = 0.8) +
  coord_flip() +
  scale_fill_brewer(palette = "Set2") +
  labs(
    title = "Proportion of Variance Explained (R² Ratio)",
    x = "Variable",
    y = "R² Ratio",
    fill = "Variable"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  ) +
  scale_y_continuous(labels = scales::percent_format())

anova(mod, permutations = 999)
anova(mod, by = "axis", permutations = 999)
anova(mod, by = "terms", permutations = 999)

sc1 <- summary(mod)[["sites"]][, 1:2]

H   <- Hpi(x = sc1)
est <- kde(x = sc1, H = H, compute.cont = TRUE)

minx = min(sc1[,1]) - 0.5
maxx = max(sc1[,1]) + 0.5
miny = min(sc1[,2]) - 0.5
maxy = max(sc1[,2]) + 0.5

color_set <- colorRampPalette(c("white", "#d9dd2b", "#84c253", "#39af77", "#219781", "#267c8a", "#375e86", "#3c4883", "#422769", "#331f41"))(130)
plot(est, 
     cont = seq(1, 100, by = 1), 
     display = "filled.contour2",
     col = color_set,
     add = FALSE, 
     xlab = "RDA1", 
     ylab = "RDA2",
     cex.axis = 1,
     las = 1,
     tck = 0.01,
     xaxs = 'i',
     yaxs = 'i',
     xaxt = 'n',
     yaxt = 'n',
     xlim = c(-1.5, 3.2),
     ylim = c(-1.2, 2.5))

axis(1, at = seq(round(minx), round(maxx)), labels = seq(round(minx), round(maxx)), 
     cex.axis = 1.2, las = 1, tck = -0.012, cex.lab = 3)
axis(2, at = seq(round(miny), round(maxy)), labels = seq(round(miny), round(maxy)),
     cex.axis = 1.2, las = 1, tck = -0.012, cex.lab = 3)

if("fam" %in% colnames(dat)) {
  site_names <- dat$fam[as.numeric(rownames(sc1))]
} else {
  site_names <- rownames(sc1)
}

cl <- contourLevels(est, prob = c(0.9, 0.7, 0.5, 0.3, 0.2, 0.001), approx = TRUE)

plot(est, abs.cont = cl[3], labels = c(0.5),
     add = TRUE, lwd = 0.75, lty = 1, col = "#707070", labcex = 0.75)
plot(est, abs.cont = cl[5], labels = c(0.8),
     add = TRUE, lwd = 0.5, lty = 1, col = "#707070", labcex = 0.75)
plot(est, abs.cont = cl[6], labels = c(0.99),
     add = TRUE, lwd = 0.5, lty = 1, col = "#bababa", labcex = 0.75)
plot(est, abs.cont = cl[1], drawlabels = FALSE,
     add = TRUE, lwd = 0.3, lty = 3, col = "#707070")
plot(est, abs.cont = cl[2], drawlabels = FALSE,
     add = TRUE, lwd = 0.3, lty = 3, col = "#707070")
plot(est, abs.cont = cl[4], drawlabels = FALSE,
     add = TRUE, lwd = 0.3, lty = 3, col = "#707070")

abline(h = 0, lty = 2, lwd = 0.5, col = "grey60")
abline(v = 0, lty = 2, lwd = 0.5, col = "grey60")

dat_xy      <- mod$CCA$v[, 1:2]
dat_xy_text <- dat_xy * 1.7
dat_xy_arr  <- dat_xy * 1.6

for (i in 1:nrow(dat_xy_arr)) {
  arrows(x0 = 0, y0 = 0, 
         x1 = dat_xy_arr[i, 1], y1 = dat_xy_arr[i, 2],
         col = "white", length = 0.05, lwd = 1.2)
}
text(dat_xy_text[, 1], dat_xy_text[, 2], 
     col = "white", cex = 0.8, labels = colnames(y1))

dat_xy1      <- mod$CCA$biplot[, 1:2]
dat_xy1_text <- dat_xy1 * 3.2
dat_xy1_arr  <- dat_xy1 * 3

for (i in 1:nrow(dat_xy1_arr)) {
  arrows(x0 = 0, y0 = 0, 
         x1 = dat_xy1_arr[i, 1], y1 = dat_xy1_arr[i, 2],
         col = 'black', length = 0.05, lwd = 1.2)
}
text(dat_xy1_text[, 1], dat_xy1_text[, 2], 
     col = "black", cex = 0.8, labels = colnames(x))
