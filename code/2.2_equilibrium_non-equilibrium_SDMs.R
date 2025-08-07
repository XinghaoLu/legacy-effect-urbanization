rm(list = ls())

source("/241009/code/2.1_weighting_funct.R")
library(dplyr)
library(scalescape)
library(cvAUC)
library(ecospat)

### Read inputs data
China.grid.comp <- terra::vect("/241009/data/02_timelag/inputs/china/china.comp.shp")
China.grid.comp
length(China.grid.comp)

# Read species data
occ.comp <- read.table("/241009/data/02_timelag/inputs/china/occ.comp.txt", sep=",", header = T)

length(table(occ.comp$scientificName))
species.list <- unique(occ.comp$scientificName)

# Read variables
veg.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/veg1990-2020.csv", sep=",", header = T)
wat.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/wat1990-2020.csv", sep=",", header = T)
shdi.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/shdi1990-2020.csv", sep=",", header = T)
frac.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/frac1990-2020.csv", sep=",", header = T)
contag.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/contag1990-2020.csv", sep=",", header = T)
ntl.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/ntl1990-2020.csv", sep=",", header = T)
tmp.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/tmp1990-2020.csv", sep=",", header = T)
pop.comp <- read.table("/241009/data/02_timelag/inputs/temporal_matrix/pop1990-2020.csv", sep=",", header = T)

colnames(veg.comp)[1] -> "dist"
colnames(wat.comp)[1] -> "dist"
colnames(shdi.comp)[1] -> "dist"
colnames(frac.comp)[1] -> "dist"
colnames(contag.comp)[1] -> "dist"
colnames(ntl.comp)[1] -> "dist"
colnames(tmp.comp)[1] -> "dist"
colnames(pop.comp)[1] -> "dist"

# Normalize variables
veg.comp.norm <- cbind( veg.comp[,1],
                        (veg.comp[,-1] - mean(as.vector(as.matrix(veg.comp[dim(veg.comp)[1],-1]))))
                        / sd(as.vector(as.matrix(veg.comp[dim(veg.comp)[1],-1]))) )

wat.comp.norm <- cbind( wat.comp[,1],
                        (wat.comp[,-1] - mean(as.vector(as.matrix(wat.comp[dim(wat.comp)[1],-1]))))
                        / sd(as.vector(as.matrix(wat.comp[dim(wat.comp)[1],-1]))) )

shdi.comp.norm <- cbind( shdi.comp[,1],
                         (shdi.comp[,-1] - mean(as.vector(as.matrix(shdi.comp[dim(shdi.comp)[1],-1]))))
                         / sd(as.vector(as.matrix(shdi.comp[dim(shdi.comp)[1],-1]))) )

frac.comp.norm <- cbind( frac.comp[,1],
                         (frac.comp[,-1] - mean(as.vector(as.matrix(frac.comp[dim(frac.comp)[1],-1]))))
                         / sd(as.vector(as.matrix(frac.comp[dim(frac.comp)[1],-1]))) )

contag.comp.norm <- cbind( contag.comp[,1],
                           (contag.comp[,-1] - mean(as.vector(as.matrix(contag.comp[dim(contag.comp)[1],-1]))))
                           / sd(as.vector(as.matrix(contag.comp[dim(contag.comp)[1],-1]))) )

ntl.comp.norm <- cbind( ntl.comp[,1],
                        (ntl.comp[,-1] - mean(as.vector(as.matrix(ntl.comp[dim(ntl.comp)[1],-1]))))
                        / sd(as.vector(as.matrix(ntl.comp[dim(ntl.comp)[1],-1]))) )

tmp.comp.norm <- cbind( tmp.comp[,1],
                        (tmp.comp[,-1] - mean(as.vector(as.matrix(tmp.comp[dim(tmp.comp)[1],-1]))))
                        / sd(as.vector(as.matrix(tmp.comp[dim(tmp.comp)[1],-1]))) )

pop.comp.norm <- cbind( pop.comp[,1],
                        (pop.comp[,-1] - mean(as.vector(as.matrix(pop.comp[dim(pop.comp)[1],-1]))))
                        / sd(as.vector(as.matrix(pop.comp[dim(pop.comp)[1],-1]))) )

### Fit temporally-weighted regressions and compare with a null non-temporal model

# Create a list with results to export 
tw.reg.out <- array(NA, dim=c(length(species.list), 23)) 
colnames(tw.reg.out) <- c("species", "Intercept", "veg.comp", "wat.comp", "shdi.comp", "frac.comp", "contag.comp", "ntl.comp", "tmp.comp", "pop.comp",
                          "range.veg.comp", "range.wat.comp", "range.shdi.comp", "range.frac.comp", "range.contag.comp", "range.ntl.comp", "range.tmp.comp", "range.pop.comp",
                          "AUC.mean", "thres.mean", "sens.sp", "spe.sp", "thres.mean.null")  

error_log <- data.frame(
  species = character(),
  error_message = character(),
  stringsAsFactors = FALSE
)

# Iterate following lines for each species
tt <- Sys.time()
for(it in 1:length(species.list)){
  tryCatch({
    print(it)
    
    # Extract occurrence for a given species
    occ.sp <-
      occ.comp %>%
      filter(scientificName == species.list[it])
    occ.sp.grid <- rep(0, length(China.grid.comp$gridID))
    occ.sp.grid[ which(China.grid.comp$gridID %in% occ.sp$gridID == T) ] <- 1
    occ.sp.grid <- occ.sp.grid[colSums(is.na(veg.comp[,-1]))==0 & colSums(is.na(ntl.comp[,-1]))==0]
    
    # Fit the null model without temporally-weighted variables
    dat0 <- data.frame(cbind(occ.sp.grid,
                             as.vector(as.matrix(veg.comp.norm[31,-1])),
                             as.vector(as.matrix(wat.comp.norm[31,-1])),
                             as.vector(as.matrix(shdi.comp.norm[31,-1])),
                             as.vector(as.matrix(frac.comp.norm[31,-1])),
                             as.vector(as.matrix(contag.comp.norm[31,-1])),
                             as.vector(as.matrix(ntl.comp.norm[31,-1])),
                             as.vector(as.matrix(tmp.comp.norm[31,-1])),
                             as.vector(as.matrix(pop.comp.norm[31,-1]))
    ))
    colnames(dat0) <- c("occ.sp.grid", "veg", "wat", "shdi", "frac", "contag", "ntl", "tmp", "pop")
    mod0 <- glm(occ.sp.grid ~ veg + wat + shdi + frac + contag + ntl + tmp + pop, data = dat0, family = binomial(link = "logit"))
    
    # Fit the full model with temporally-weighted variables
    mod.exp <- scalescape::dist_weight(mod0 = mod0,
                                       landscape.vars = list(
                                         veg = as.matrix(veg.comp.norm),
                                         wat = as.matrix(wat.comp.norm),
                                         shdi = as.matrix(shdi.comp.norm),
                                         frac = as.matrix(frac.comp.norm),
                                         contag = as.matrix(contag.comp.norm),
                                         ntl = as.matrix(ntl.comp.norm),
                                         tmp = as.matrix(tmp.comp.norm),
                                         pop = as.matrix(pop.comp.norm)
                                       ), landscape.formula =
                                         'occ.sp.grid ~ veg + wat + shdi + frac + contag + ntl + tmp + pop', lower = 1,
                                       data = dat0, plot.fits = F, weight.fn = "Gaussian")
    
    summary(mod.exp)
    
    # AUC cross-validation & threshold definition
    {
      seq.id=seq(1,dim(dat0)[1],1)
      samples=replicate(100,sample(seq.id,length(seq.id),replace = F))
      dim(samples)
      
      auc = c()
      thres = c()
      sens = c()
      spe = c()
      thres.null = c()
      dat70perc = ceiling(70*length(occ.sp.grid)/100) 
      
      for(i in 1:20){
        
        occ.sp.grid.cv = occ.sp.grid[samples[1:dat70perc,i]] 
        veg.comp.norm.cv = veg.comp.norm[,c(1,samples[1:dat70perc,i])] 
        wat.comp.norm.cv = wat.comp.norm[,c(1,samples[1:dat70perc,i])]
        shdi.comp.norm.cv = shdi.comp.norm[,c(1,samples[1:dat70perc,i])]
        frac.comp.norm.cv = frac.comp.norm[,c(1,samples[1:dat70perc,i])]
        contag.comp.norm.cv = contag.comp.norm[,c(1,samples[1:dat70perc,i])]
        ntl.comp.norm.cv = ntl.comp.norm[,c(1,samples[1:dat70perc,i])]
        tmp.comp.norm.cv = tmp.comp.norm[,c(1,samples[1:dat70perc,i])]
        pop.comp.norm.cv = pop.comp.norm[,c(1,samples[1:dat70perc,i])] 
        
        # Fit the null model without temporally-weighted variables
        dat0.cv <- data.frame(cbind(occ.sp.grid.cv,
                                    as.vector(as.matrix(veg.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(wat.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(shdi.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(frac.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(contag.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(ntl.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(tmp.comp.norm.cv[31,-1])),
                                    as.vector(as.matrix(pop.comp.norm.cv[31,-1]))
        )) 
        colnames(dat0.cv) <- c("occ.sp.grid", "veg", "wat", "shdi", "frac", "contag", "ntl", "tmp", "pop") 
        mod0.cv <- glm(occ.sp.grid ~ veg + wat + shdi + frac + contag + ntl + tmp + pop, data = dat0.cv, family = binomial(link = "logit")) 
        
        # Fit the full model with temporally-weighted variables
        mod.exp.cv <- scalescape::dist_weight(mod0 = mod0.cv,
                                              landscape.vars = list(
                                                veg = as.matrix(veg.comp.norm.cv),
                                                wat = as.matrix(wat.comp.norm.cv),
                                                shdi = as.matrix(shdi.comp.norm.cv),
                                                frac = as.matrix(frac.comp.norm.cv),
                                                contag = as.matrix(contag.comp.norm.cv),
                                                ntl = as.matrix(ntl.comp.norm.cv),
                                                tmp = as.matrix(tmp.comp.norm.cv),
                                                pop = as.matrix(pop.comp.norm.cv)
                                              ), landscape.formula =
                                                'occ.sp.grid ~ veg + wat + shdi + frac + contag + ntl + tmp + pop', lower = 1,
                                              data = dat0.cv, plot.fits = F, weight.fn = "Gaussian")
        
        # Validate the full model & predict
        dat0.pred <- data.frame(cbind(
          occ.sp.grid [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(veg.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(wat.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(shdi.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(frac.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(contag.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(ntl.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(tmp.comp.norm[31,-1])) [samples[(dat70perc+1):length(occ.sp.grid),i]],
          as.vector(as.matrix(pop.comp.norm[31,-1]) [samples[(dat70perc+1):length(occ.sp.grid),i]])
        )) 
        colnames(dat0.pred) <- c("occ.sp.grid", "veg", "wat", "shdi", "frac", "contag", "ntl", "tmp", "pop") 
        
        max.Dist = c(
          max(as.matrix(veg.comp.norm.cv)[,1]),
          max(as.matrix(wat.comp.norm.cv)[,1]),
          max(as.matrix(shdi.comp.norm.cv)[,1]),
          max(as.matrix(frac.comp.norm.cv)[,1]),
          max(as.matrix(contag.comp.norm.cv)[,1]),
          max(as.matrix(ntl.comp.norm.cv)[,1]),
          max(as.matrix(tmp.comp.norm.cv)[,1]),
          max(as.matrix(pop.comp.norm.cv)[,1])
        ) 
        
        data.pred.cv <- weighting_funct(par = mod.exp.cv$opt.range/max.Dist, mod0 = mod0.cv, landscape.formula = mod.exp.cv$landscape.formula,
                                        data = dat0.pred,
                                        landscape.vars = list(
                                          veg = as.matrix(veg.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          wat = as.matrix(wat.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          shdi = as.matrix(shdi.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          frac = as.matrix(frac.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          contag = as.matrix(contag.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          ntl = as.matrix(ntl.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          tmp = as.matrix(tmp.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])]),
                                          pop = as.matrix(pop.comp.norm[,c(1,samples[(dat70perc+1):length(occ.sp.grid),i])])
                                        ),
                                        max.Dist = max.Dist,
                                        weight.fn = "Gaussian", return.coef = T)
        
        pred.valid.cv=predict(mod.exp.cv$mod,newdata = data.pred.cv$data.with.landscape, type = "response") 
        auc <- c(auc,cvAUC::cvAUC(predictions = pred.valid.cv,
                                  labels = occ.sp.grid[samples[(dat70perc+1):length(occ.sp.grid),i]])$cvAUC) 
        thres <- c(thres, ecospat::ecospat.max.tss(pred.valid.cv, occ.sp.grid[samples[(dat70perc+1):length(occ.sp.grid),i]])$max.threshold) 
        
        prd <- rep(0, length(pred.valid.cv))
        prd[which(pred.valid.cv > ecospat::ecospat.max.tss(pred.valid.cv, occ.sp.grid[samples[(dat70perc+1):length(occ.sp.grid),i]])$max.threshold)] <- 1 
        conf <- mda::confusion(prd, occ.sp.grid[samples[(dat70perc+1):length(occ.sp.grid),i]]) 
        sens <- c(sens, conf[4] / (conf[4] + conf[3])) 
        spe <- c(spe, conf[1] / (conf[2] + conf[1])) 
        
        # Define threshold for occurrence prediction with the null model
        pred0.valid.cv = predict(mod0,newdata = dat0.pred, type = "response") 
        thres.null <- c(thres.null, ecospat::ecospat.max.tss(pred0.valid.cv, occ.sp.grid[samples[(dat70perc+1):length(occ.sp.grid),i]])$max.threshold) 
        
      }
      
      AUC.sp <- mean(auc) 
      thres.sp.mean.null <- (mean(thres.null)) 
      sens.sp <- mean(sens) 
      spe.sp <- mean(spe) 
      thres.sp.mean <- (mean(thres)) 
      
    }
    
    # Non-equilibrium predictions
    data.pred <- weighting_funct(par = mod.exp$opt.range/mod.exp$max.Dist, mod0 = mod.exp$mod0, landscape.formula = mod.exp$landscape.formula,
                                 data = mod.exp$mod0$data,
                                 landscape.vars = list(
                                   veg = as.matrix(veg.comp.norm),
                                   wat = as.matrix(wat.comp.norm),
                                   shdi = as.matrix(shdi.comp.norm),
                                   frac = as.matrix(frac.comp.norm),
                                   contag = as.matrix(contag.comp.norm),
                                   ntl = as.matrix(ntl.comp.norm),
                                   tmp = as.matrix(tmp.comp.norm),
                                   pop = as.matrix(pop.comp.norm)
                                 ),
                                 
                                 max.Dist = mod.exp$max.Dist,
                                 weight.fn = "Gaussian", return.coef = T)
    
    pred.sp=predict(mod.exp$mod,newdata = data.pred$data.with.landscape, type = "response", se.fit = T)$fit
    out.fil <- paste("/241009/data/02_timelag/outputs/temporally_weighted_regressions_gaussian/",gsub(" ", "_", species.list[it]),"_predictions.csv", sep="")
    
    # Equilibrium predictions
    data.pred0 <- data.frame(as.vector(as.matrix(veg.comp.norm[(31),-1])),
                             as.vector(as.matrix(wat.comp.norm[(31),-1])),
                             as.vector(as.matrix(shdi.comp.norm[(31),-1])),
                             as.vector(as.matrix(frac.comp.norm[(31),-1])),
                             as.vector(as.matrix(contag.comp.norm[(31),-1])),
                             as.vector(as.matrix(ntl.comp.norm[(31),-1])),
                             as.vector(as.matrix(tmp.comp.norm[(31),-1])),
                             as.vector(as.matrix(pop.comp.norm[(31),-1]))
    ) 
    colnames(data.pred0) <- c("veg", "wat", "shdi", "frac", "contag", "ntl", "tmp", "pop") 
    pred0.sp=predict(mod.exp$mod0,newdata = data.pred0, type = "response", se.fit = T)$fit 
    out.fil0 <- paste("/241009/data/02_timelag/outputs/temporally_weighted_regressions_gaussian/", gsub(" ", "_", species.list[it]), "_predictions0.csv", sep="")   
    
    # Export results
    out.fil.rda <- paste("/241009/data/02_timelag/outputs/temporally_weighted_regressions_gaussian/", species.list[it], ".rds", sep="") 
    out.fil.rda <- gsub(" ", "_", out.fil.rda)
    saveRDS(mod.exp, file = out.fil.rda)
    spi.out <- c(species.list[it], mod.exp$mod$coef, mod.exp$opt.range, AUC.sp,thres.sp.mean, sens.sp,spe.sp,thres.sp.mean.null) 
    tw.reg.out[it,] <- spi.out 
    write.table(tw.reg.out,"/241009/data/02_timelag/outputs/temporally_weighted_regressions_gaussian/temporally_weighted_regressions_sp.csv", col.names = T, sep=";")
    write.table(pred.sp,file = out.fil, sep=";")
    write.table(pred0.sp,file = out.fil0, sep=";")
  }, error = function(e) {
    
    error_log <<- rbind(error_log, data.frame(
      species = species.list[it],
      error_message = as.character(e),
      stringsAsFactors = FALSE
    ))
    
    message(sprintf("Error processing species %s: %s", species.list[it], e))
    
    return(NULL)
  })
}