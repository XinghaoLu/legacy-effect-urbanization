rm(list = ls())

library(picante)

phy <- read.tree("/2410009/data/03_comm_legacy/3.1_phylogenetic_legacy_res/tree.newick")
summary(phy$edge)
phy <- compute.brlen(phy)

# Load comm0unity data
comm0 <- read.csv("/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/binary.ts.csv",
                  header = TRUE, 
                  sep = ',',
                  row.names = 1)

comm0 <- comm0[rowSums(comm0) > 0, ]

combined <- match.phylo.comm(phy,comm0)
comm0.pd <- pd(comm0,phy,include.root=T)

# Load comm1unity data
comm1 <- read.csv("/2410009/data/03_comm_legacy/3.1_taxonomic_legacy_res/binary.ts.csv",
                  header = TRUE, 
                  sep = ',',
                  row.names = 1)

comm1 <- comm1[rowSums(comm1) > 0, ]

combined <- match.phylo.comm(phy,comm1)
comm1.pd <- pd(comm1,phy,include.root=T)

write.csv(comm0.pd, "/2410009/data/03_comm_legacy/3.1_phylogenetic_legacy_res/phylogenetic0.ts.csv")
write.csv(comm1.pd, "/2410009/data/03_comm_legacy/3.1_phylogenetic_legacy_res/phylogenetic.ts.csv")