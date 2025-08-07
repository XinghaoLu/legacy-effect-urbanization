rm(list = ls())

library(mFD)
library(dplyr)
library(tidyr)

traits.sp <- read.csv("/2410009/data/03_comm_legacy/3.2_functional_legacy_res/traits.csv", header = TRUE, row.names = 1)
traits.type <- read.csv("/2410009/data/03_comm_legacy/3.2_functional_legacy_res/traits.csv", header = TRUE)

traits.sp$Habitat.Density <- factor(
  traits.sp$Habitat.Density,
  levels = c("open", "semi", "dense"),
  ordered = TRUE
)

traits.sp$Beak.Length_Culmen <- log(traits.sp$Beak.Length_Culmen)
traits.sp$urtolerance.1 <- log(traits.sp$urtolerance.1)
traits.sp$Hand.Wing.Index <- log(traits.sp$Hand.Wing.Index)
traits.sp$Tail.Length <- log(traits.sp$Tail.Length)
traits.sp$Clutch.size <- log(traits.sp$Clutch.size)
traits.sp$Ave.generation.length <- log(traits.sp$Ave.generation.length)

comm <- read.csv("/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/binary.ts.csv", header = TRUE)
comm$gridID <- paste0("grid_", comm$gridID)

comm0 <- read.csv("/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/binary0.ts.csv", header = TRUE)
comm0$gridID <- paste0("grid0_", comm0$gridID)

commall <- rbind(comm, comm0)

comm_matrix <- as.matrix(commall[, -1])
rownames(comm_matrix) <- commall$gridID

zero_cols <- which(colSums(comm_matrix) == 0)
if(length(zero_cols) > 0){
  message(" ", paste(colnames(comm_matrix)[zero_cols], collapse = ", "))
  comm_matrix <- comm_matrix[, -zero_cols]
}

common_species <- intersect(rownames(traits.sp), colnames(comm_matrix))
traits.sp <- traits.sp[common_species, ]
comm_matrix <- comm_matrix[, common_species]

sp_dist <- mFD::funct.dist(
  sp_tr = traits.sp,
  tr_cat = traits.type, 
  metric = "gower",
  scale_euclid = "scale_center", 
  ordinal_var = "classic",
  weight_type = "equal",
  stop_if_NA = TRUE)

fspaces_quality <- mFD::quality.fspaces(sp_dist = sp_dist, 
                                        maxdim_pcoa = 10, 
                                        deviation_weighting = "absolute",
                                        fdist_scaling = FALSE,
                                        fdendro = "average")

round(fspaces_quality$"quality_fspaces", 5)

mFD::quality.fspaces.plot(fspaces_quality = fspaces_quality, 
                          quality_metric = "mad", 
                          fspaces_plot = c("tree_average", "pcoa_1d", "pcoa_2d", "pcoa_3d", "pcoa_4d", "pcoa_5d", "pcoa_6d", "pcoa_7d", "pcoa_8d"), 
                          name_file = NULL, 
                          range_dist = NULL,
                          range_dev = NULL, 
                          range_qdev = NULL,
                          gradient_deviation = c(neg="#5f3a91",nul="grey80",pos="#a04d07"),
                          gradient_deviation_quality = c(low="#FFF9BD",high="#31AA77"),
                          x_lab = "Trait-based distance")

sp_faxes_coord <- fspaces_quality$"details_fspaces"$"sp_pc_coord"

valid_asb <- names(which(rowSums(comm_matrix > 0) > 5))
alpha_fd <- mFD::alpha.fd.multidim(
  sp_faxes_coord = sp_faxes_coord[,c("PC1", "PC2", "PC3", "PC4", "PC5")],
  asb_sp_w = comm_matrix[valid_asb, ],
  ind_vect = c("fric"),
  scaling = TRUE, 
  check_input = TRUE, 
  details_returned = TRUE
)

fd_ind_values<-alpha_fd$"functional_diversity_indices"
fd_ind_values

max_id <- nrow(comm)

fd_ind_values$id <- gsub("grid0?_", "", rownames(fd_ind_values))
fd_ind_values$format <- ifelse(
  grepl("grid0_", rownames(fd_ind_values)), "grid0", "grid"
)

grid_data <- fd_ind_values %>% 
  filter(format == "grid") %>% 
  select(id, fric) %>% 
  rename(fric_grid = fric)

grid0_data <- fd_ind_values %>% 
  filter(format == "grid0") %>% 
  select(id, fric) %>% 
  rename(fric_grid0 = fric)

result <- full_join(grid_data, grid0_data, by = "id") %>% 
  mutate(id = as.numeric(id)) %>%
  arrange(id)

final_result <- result %>%
  complete(id = 1:max_id)


write.csv(final_result, "/2410009/data/03_comm_legacy/3.2_functional_legacy_res/functional.ts.csv")
