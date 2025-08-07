rm(list = ls())

library(dplyr)
library(ggplot2)
library(tidyr)
library(fields)

### Map extinction debt and immigration credit

## Load grid & country delineation
China.grid.comp <- terra::vect("/2410009/data/02_timelag/inputs/China/China.comp.shp")
veg.comp <- read.table("/2410009/data/02_timelag/inputs/Temporal_matrix/veg1990-2020.csv", sep=",", header = T) # to handle NA values
delin=geodata::world(resolution=5, level=0, path="/2410009/data/02_timelag/inputs/", version="latest")
delin.3035 <- terra::project(delin, terra::crs(China.grid.comp))

## Predict richness
input_dir <- "/2410009/data/02_timelag/outputs/Temporally_weigthed_regressions_gaussian/"

all_files <- list.files(input_dir, full.names = FALSE)
fi <- all_files[grep("_predictions\\.csv$", all_files)]
fi0 <- all_files[grep("_predictions0\\.csv$", all_files)]

file_species <- gsub("_predictions\\.csv$", "", fi) %>% 
  gsub("_", " ", .) %>% 
  trimws()

valid_idx <- which(file_species %in% spi.out.red$species)
fi_valid <- fi[valid_idx]
fi0_valid <- fi0[valid_idx]
species_valid <- file_species[valid_idx]

if(length(fi_valid) == 0) stop("")
message(sprintf("", length(fi), length(fi_valid)))

first_pred <- read.table(paste0(input_dir, fi_valid[1]), sep = ";")
richness.ts <- array(0, dim = dim(first_pred))
richness0.ts <- array(0, dim = dim(first_pred))

for(i in seq_along(fi_valid)) {
  species_row <- which(spi.out.red$species == species_valid[i])
  pred <- read.table(paste0(input_dir, fi_valid[i]), sep = ";")
  pred0 <- read.table(paste0(input_dir, fi0_valid[i]), sep = ";")
  pred.bin <- (pred > spi.out.red$thres.mean[species_row]) * 1L
  pred0.bin <- (pred0 > spi.out.red$thres.mean.null[species_row]) * 1L
  richness.ts <- richness.ts + pred.bin
  richness0.ts <- richness0.ts + pred0.bin
  if(i %% 10 == 0) {
    message(sprintf("process: %d/%d [%s]", 
                    i, length(fi_valid), 
                    species_valid[i]))
  }
}

write.csv(richness.ts, "/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/richness.ts.csv")
write.csv(richness0.ts, "/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/richness0.ts.csv")

site_ids <- China.grid.comp$gridID

result_df_bin <- data.frame(gridID = site_ids)
result_df_bin0 <- data.frame(gridID = site_ids)

for(i in seq_along(fi_valid)) {
  species_row <- which(spi.out.red$species == species_valid[i])
  species_colname <- gsub(" ", "_", species_valid[i])
  pred_file <- paste0(input_dir, fi_valid[i])
  pred0_file <- paste0(input_dir, fi0_valid[i])
  pred <- read.table(pred_file, sep = ";", header = TRUE)
  pred0 <- read.table(pred0_file, sep = ";", header = TRUE)
  if(nrow(pred) != length(site_ids) | nrow(pred0) != length(site_ids)) {
    stop(sprintf(" ", species_valid[i]))
  }
  result_df_bin[[species_colname]] <- (pred[,1] > spi.out.red$thres.mean[species_row]) * 1L
  result_df_bin0[[species_colname]] <- (pred0[,1] > spi.out.red$thres.mean.null[species_row]) * 1L

  if(i %% 10 == 0) {
    message(sprintf(" ", 
                    i, length(fi_valid), 
                    species_valid[i]))
  }
}

write.csv(result_df_bin, "/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/binary.ts.csv", row.names = FALSE)
write.csv(result_df_bin0, "/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/binary0.ts.csv", row.names = FALSE)