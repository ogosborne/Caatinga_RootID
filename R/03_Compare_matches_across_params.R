library(UpSetR)
# load optimal parameter values
optimal <- read.csv("results/parameter_comparison/optimal_params.csv", header = T)
# load RootID matches
load("results/processed_data/match.list.RData")
# all parameter values
US.M.vals <- seq(2,8,2)
CS.n.vals <- c(-1,0,1)
maxmiss.vals <- cbind(c(0,0.2,0.2,0.2,0.4,0.4,0.4),c(0,0,0.1,0.2,0,0.1,0.2))
mindepth.vals <-  c(5,10,20)
maxhaps.vals <- c(2,3,NA)
# optimal parameter values
best.M <- optimal[1, "US_M"]
best.n <- optimal[1, "CS_n"] - optimal[1, "US_M"]
best.mdm <- optimal[1, "FD_max.md.marker"]
best.mdh <- optimal[1, "FD_max.md.hap"]
best.m <- optimal[1, "FD_min.dep"]
best.mh <- optimal[1, "FD_max.haps"]
### UpSet STACKS parameters
dir.create("results/parameter_comparison/STACKS/")
my.list.spp <- list()
my.list.ind <- list()
for(M in US.M.vals){
  for(n in CS.n.vals){
    my.dat  <- match.list[[paste("US_M_", M, sep = "")]][[paste("CS_n_", M+n, sep = "")]][[paste("FD_max.md.marker_", best.mdm, sep = "")]][[paste("FD_max.md.hap_", best.mdh, sep = "")]][[paste("FD_min.dep_", best.m, sep = "")]][[paste("FD_max.haps_", as.character(best.mh), sep = "")]]
    # species
    spp.mlist <- list()
    for(sp in colnames(my.dat$species.diag.markers)[-1]){
      spp.mlist[[sp]] <- paste(sp, my.dat$species.diag.markers[which(my.dat$species.diag.markers[, sp] > 0), "sample"], sep = "_")
    }
    spp.mlist <- do.call(c, spp.mlist)
    my.list.spp[[paste("M = ", M, ";n = ", M+n, sep = "")]] <- spp.mlist
    # individuals
    ind.mlist <- list()
    for(sp in names(my.dat$individual.diag.haplotypes)){
      for(ind in colnames(my.dat$individual.diag.haplotypes[[sp]])[2:(ncol(my.dat$individual.diag.haplotypes[[sp]])-1)]){
        ind.mlist[[ind]] <- paste(sp, my.dat$individual.diag.haplotypes[[sp]][which(my.dat$individual.diag.haplotypes[[sp]][, ind] > 0), "sample"], sep = "_")
      }
    }
    ind.mlist <-  do.call(c, ind.mlist)
    my.list.ind[[paste("M = ", M, ";n = ", M+n, sep = "")]] <- ind.mlist
  }
}
# species pdf
pdf(paste("results/parameter_comparison/STACKS/comp_STACKS_mdm", best.mdm, "_mdh", best.mdh, "_min.dep", best.m, "_max.haps", as.character(best.mh), "_species_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.spp), nintersects = NA, nsets = length(my.list.spp), order.by = "freq", text.scale = 1.5, sets = names(my.list.spp), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()
# individual pdf
pdf(paste("results/parameter_comparison/STACKS/comp_STACKS_mdm", best.mdm, "_mdh", best.mdh, "_min.dep", best.m, "_max.haps", as.character(best.mh), "_individual_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.ind), nintersects = NA, nsets = length(my.list.ind), order.by = "freq", text.scale = 1.5, sets = names(my.list.ind), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()

### UpSet Missing data
dir.create("results/parameter_comparison/missing_data/")
my.list.spp <- list()
my.list.ind <- list()
for(md in 1:nrow(maxmiss.vals)){
  my.dat  <- match.list[[paste("US_M_", best.M, sep = "")]][[paste("CS_n_", best.M+best.n, sep = "")]][[paste("FD_max.md.marker_", maxmiss.vals[md, 1], sep = "")]][[paste("FD_max.md.hap_", maxmiss.vals[md, 2], sep = "")]][[paste("FD_min.dep_", best.m, sep = "")]][[paste("FD_max.haps_", as.character(best.mh), sep = "")]]
  # species
  spp.mlist <- list()
  for(sp in colnames(my.dat$species.diag.markers)[-1]){
    spp.mlist[[sp]] <- paste(sp, my.dat$species.diag.markers[which(my.dat$species.diag.markers[, sp] > 0), "sample"], sep = "_")
  }
  spp.mlist <- do.call(c, spp.mlist)
  my.list.spp[[paste("marker = ", maxmiss.vals[md, 1], ";haplotype = ", maxmiss.vals[md, 2], sep = "")]] <- spp.mlist
  # individuals
  ind.mlist <- list()
  for(sp in names(my.dat$individual.diag.haplotypes)){
    for(ind in colnames(my.dat$individual.diag.haplotypes[[sp]])[2:(ncol(my.dat$individual.diag.haplotypes[[sp]])-1)]){
      ind.mlist[[ind]] <- paste(sp, my.dat$individual.diag.haplotypes[[sp]][which(my.dat$individual.diag.haplotypes[[sp]][, ind] > 0), "sample"], sep = "_")
    }
  }
  ind.mlist <-  do.call(c, ind.mlist)
  my.list.ind[[paste("marker = ", maxmiss.vals[md, 1], ";haplotype = ", maxmiss.vals[md, 2], sep = "")]] <- ind.mlist
}

# species pdf
pdf(paste("results/parameter_comparison/missing_data/comp_missing_data_M", best.M, "_n", best.M+best.n, "_min.dep", best.m, "_max.haps", as.character(best.mh), "_species_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.spp), nintersects = NA, nsets = length(my.list.spp), order.by = "freq", text.scale = 1.5, sets = names(my.list.spp), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()
# individual pdf
pdf(paste("results/parameter_comparison/missing_data/comp_missing_data_M", best.M, "_n", best.M+best.n, "_min.dep", best.m, "_max.haps", as.character(best.mh), "_individual_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.ind), nintersects = NA, nsets = length(my.list.ind), order.by = "freq", text.scale = 1.5, sets = names(my.list.ind), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()


### UpSet Min depth
dir.create("results/parameter_comparison/min_depth/")
my.list.spp <- list()
my.list.ind <- list()
for(m in mindepth.vals){
  my.dat  <- match.list[[paste("US_M_", best.M, sep = "")]][[paste("CS_n_", best.M+best.n, sep = "")]][[paste("FD_max.md.marker_", best.mdm, sep = "")]][[paste("FD_max.md.hap_", best.mdh, sep = "")]][[paste("FD_min.dep_", m, sep = "")]][[paste("FD_max.haps_", as.character(best.mh), sep = "")]]
  # species
  spp.mlist <- list()
  for(sp in colnames(my.dat$species.diag.markers)[-1]){
    spp.mlist[[sp]] <- paste(sp, my.dat$species.diag.markers[which(my.dat$species.diag.markers[, sp] > 0), "sample"], sep = "_")
  }
  spp.mlist <- do.call(c, spp.mlist)
  my.list.spp[[paste("min.dep = ", m, sep = "")]] <- spp.mlist
  # individuals
  ind.mlist <- list()
  for(sp in names(my.dat$individual.diag.haplotypes)){
    for(ind in colnames(my.dat$individual.diag.haplotypes[[sp]])[2:(ncol(my.dat$individual.diag.haplotypes[[sp]])-1)]){
      ind.mlist[[ind]] <- paste(sp, my.dat$individual.diag.haplotypes[[sp]][which(my.dat$individual.diag.haplotypes[[sp]][, ind] > 0), "sample"], sep = "_")
    }
  }
  ind.mlist <-  do.call(c, ind.mlist)
  my.list.ind[[paste("min.dep = ", m, sep = "")]] <- ind.mlist
}

# species pdf
pdf(paste("results/parameter_comparison/min_depth/comp_min_depth_M", best.M, "_n", best.M+best.n, "_mdm", best.mdm, "_mdh", best.mdh, "_max.haps", as.character(best.mh), "_species_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.spp), nintersects = NA, nsets = length(my.list.spp), order.by = "freq", text.scale = 1.5, sets = names(my.list.spp), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()
# individual pdf
pdf(paste("results/parameter_comparison/min_depth/comp_min_depth_M", best.M, "_n", best.M+best.n, "_mdm", best.mdm, "_mdh", best.mdh, "_max.haps", as.character(best.mh), "_individual_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.ind), nintersects = NA, nsets = length(my.list.ind), order.by = "freq", text.scale = 1.5, sets = names(my.list.ind), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()

### UpSet Max haps
dir.create("results/parameter_comparison/max_haps/")
my.list.spp <- list()
my.list.ind <- list()
for(mh in maxhaps.vals){
  my.dat  <- match.list[[paste("US_M_", best.M, sep = "")]][[paste("CS_n_", best.M+best.n, sep = "")]][[paste("FD_max.md.marker_", best.mdm, sep = "")]][[paste("FD_max.md.hap_", best.mdh, sep = "")]][[paste("FD_min.dep_", best.m, sep = "")]][[paste("FD_max.haps_", as.character(mh), sep = "")]]
  # species
  spp.mlist <- list()
  for(sp in colnames(my.dat$species.diag.markers)[-1]){
    spp.mlist[[sp]] <- paste(sp, my.dat$species.diag.markers[which(my.dat$species.diag.markers[, sp] > 0), "sample"], sep = "_")
  }
  spp.mlist <- do.call(c, spp.mlist)
  my.list.spp[[paste("max.haps = ", mh, sep = "")]] <- spp.mlist
  # individuals
  ind.mlist <- list()
  for(sp in names(my.dat$individual.diag.haplotypes)){
    for(ind in colnames(my.dat$individual.diag.haplotypes[[sp]])[2:(ncol(my.dat$individual.diag.haplotypes[[sp]])-1)]){
      ind.mlist[[ind]] <- paste(sp, my.dat$individual.diag.haplotypes[[sp]][which(my.dat$individual.diag.haplotypes[[sp]][, ind] > 0), "sample"], sep = "_")
    }
  }
  ind.mlist <-  do.call(c, ind.mlist)
  my.list.ind[[paste("max.haps = ", mh, sep = "")]] <- ind.mlist
}

# species pdf
pdf(paste("results/parameter_comparison/max_haps/comp_max_haps_M", best.M, "_n", best.M+best.n, "_mdm", best.mdm, "_mdh", best.mdh, "_min.dep", best.m, "_species_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.spp), nintersects = NA, nsets = length(my.list.spp), order.by = "freq", text.scale = 1.5, sets = names(my.list.spp), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()
# individual pdf
pdf(paste("results/parameter_comparison/max_haps/comp_max_haps_M", best.M, "_n", best.M+best.n, "_mdm", best.mdm, "_mdh", best.mdh, "_min.dep", best.m, "_individual_Upset.pdf", sep = ""))
# upset
print(upset(fromList(my.list.ind), nintersects = NA, nsets = length(my.list.ind), order.by = "freq", text.scale = 1.5, sets = names(my.list.ind), mb.ratio = c(0.5,  0.5), mainbar.y.label = "N shared matches", sets.x.label = "N matches", keep.order = T))
# close pdf
dev.off()

