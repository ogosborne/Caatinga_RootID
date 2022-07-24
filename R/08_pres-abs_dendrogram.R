library(ape)
library(RootID)
# load data
load("results/processed_data/best.stacks.dat.RData")
# load metadata
metadata <- read.table("data/metadata.csv",sep=",",header=T,na.strings = "NA",stringsAsFactors = F)
# make list of which individuals are from which species
ind.list <- list()
for(i in sort(unique(metadata$species))){
  ind.list[[i]] <- metadata[which(metadata$species == i),"sample"]
}
rm(i)
# make output directory
dir.create("results/dendrogram/")
# make dendrogram
ind.dend <- shared.marker.tree(stacks.dat = best.stacks.dat, root = "Neoglaziovia_variegata", ind.list = ind.list, collapse.to.species = F)
# plot dendrogram
pdf("results/dendrogram/individual_level_UPGMA.pdf")
plot.phylo(ind.dend, no.margin = T, cex=0.9)
dev.off()
