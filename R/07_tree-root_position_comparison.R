# make output directory
dir.create("results/tree_vs_root_pos/")
# load data
load("results/processed_data/best.match.md0.RData")
load("results/processed_data/best.match.md3.RData")
metadata <- read.table("data/metadata.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
# get root coordinates
root.pos <- data.frame(matrix(NA, nrow = (nrow(best.matches$species.diag.markers)), ncol = 4, dimnames = list(NULL, c("sample", "x", "y", "depth"))))
root.pos$sample <- best.matches$species.diag.markers[, "sample"]
pos <- seq(3, 11, 2)
names(pos) <- c("A", "B", "C", "D", "E")
root.pos$x <- pos[sapply(X = root.pos$sample, FUN = function(x) strsplit(strsplit(x, split = "_")[[1]][[2]], split = "")[[1]][[2]])]
root.pos$y <- pos[as.numeric(sapply(X = root.pos$sample, FUN = function(x) strsplit(strsplit(x, split = "_")[[1]][[2]], split = "")[[1]][[1]]))]
root.pos$depth <- as.numeric(sapply(X = root.pos$sample, FUN = function(x) strsplit(x, split = "_")[[1]][[3]]))
rm(pos)
# save root.pos
write.table(root.pos, file = "results/tree_vs_root_pos/root.pos.csv", sep = ",", row.names = F, quote = F)

### distance from tree to root sample for spp 
# function to get the horizontal distances between each root sample and the closest members of each species, dividing them into those that are present in the root sample and those which aren't.
make.spp.dist.list <- function(root.pos, matches, metadata){
  # make list to hold results
  root.dist.list <- list()
  root.dist.list$p <- list()
  root.dist.list$a <- list()
  # for each root sample
  for (r in 1:nrow(root.pos)){
    # make a list of present and absent species for each root sample
    r.name <- root.pos$sample[[r]]
    root.dist.list$p[[r.name]] <- list()
    root.dist.list$a[[r.name]] <- list()
    # get column numbers for the species that are present/absent
    p.n <- which(matches$species.diag.markers[r, 2:ncol(matches$species.diag.markers)] > 0) + 1
    a.n <- which(matches$species.diag.markers[r, 2:ncol(matches$species.diag.markers)] == 0) + 1
    # get the names of the species that are present/absent
    p.names <- gsub(".N.reads","",colnames(best.matches$species.diag.markers[,p.n]),fixed=T)
    a.names <- gsub(".N.reads","",colnames(best.matches$species.diag.markers[,a.n]),fixed=T)
    # get the coordinates of the root sample
    r.coords <- root.pos[r, c("x", "y")]
    # for each species which is present
    for(n in 1:length(p.n)){
      # get coordinates for all individuals of the species
      p.coords <- metadata[which(metadata$species == p.names[[n]]), c("tree.pos.x", "tree.pos.y")]
      # if coordinates are present for each tree
      if(all(!is.na(p.coords))){
        # get distances from all trees
        p.dist <- proxy::dist(r.coords, p.coords, method = "Euclidean")
        # add smallest distance to list
        root.dist.list[["p"]][[r.name]][[p.names[[n]]]] <- min(p.dist)
      }
    }
    # flatten the sublist
    root.dist.list[["p"]][[r.name]] <- unlist(root.dist.list[["p"]][[r.name]])
    # for each species which isn't present
    for(n in 1:length(a.n)){
      # get coordinates for all individuals of the species
      a.coords <- metadata[which(metadata$species == a.names[[n]]), c("tree.pos.x", "tree.pos.y")]
      # if coordinates are present for each tree
      if(all(!is.na(a.coords))){
        # get distances from all trees
        a.dist <- proxy::dist(r.coords, a.coords, method = "Euclidean")
        # add smallest distance to list
        root.dist.list[["a"]][[r.name]][[a.names[[n]]]] <- min(a.dist)
      }
      # flatten the sublist
      root.dist.list[["a"]][[r.name]] <- unlist(root.dist.list[["a"]][[r.name]])
    }
  }
  root.dist.list
}
# get distance lists (with and without a match depth filter)
sp.root.dist.list.1 <- make.spp.dist.list(root.pos = root.pos, matches = best.matches, metadata = metadata)
sp.root.dist.list.3 <- make.spp.dist.list(root.pos = root.pos, matches = best.matches.3.3, metadata = metadata)
# make boxplot of distance between root samples and present vs absent species
# min depth = 1
pdf("results/tree_vs_root_pos/root.tree.dist.species.boxplot.md1.pdf")
boxplot(list("Species detected" = do.call("c", sp.root.dist.list.1$p), "Species undetected" = do.call("c", sp.root.dist.list.1$a)), ylab = "Root-to-tree distance (m)")
dev.off()
# min depth = 3
pdf("results/tree_vs_root_pos/root.tree.dist.species.boxplot.md3.pdf")
boxplot(list("Species detected" = do.call("c", sp.root.dist.list.3$p), "Species undetected" = do.call("c", sp.root.dist.list.3$a)), ylab = "Root-to-tree distance (m)")
dev.off()
# test difference with Mann-Whitney U tests
# sink to text file
sink("results/tree_vs_root_pos/root.tree.dist.species.tests.txt")
cat("Distance from root samples to present vs absent species, min read depth = 1:","\n")
print(wilcox.test(do.call("c", sp.root.dist.list.1$p), do.call("c", sp.root.dist.list.1$a), paired = F))
cat("Distance from root samples to present vs absent species, min read depth = 3:","\n")
print(wilcox.test(do.call("c", sp.root.dist.list.3$p), do.call("c", sp.root.dist.list.3$a), paired = F))
sink()

### distance from tree to root sample for individuals
# function to get the horizontal distances between each root sample and each individual, dividing them into those that are present in the root sample and those which aren't.
make.ind.dist.list <- function(root.pos, matches, metadata){
  # concatenate individual matches for each species and remove non-diagnostic and sample names
  per.root.ind.match <- do.call("cbind", matches$individual.diag.haplotypes)
  per.root.ind.match <- per.root.ind.match[, grep(".non.diagnostic|.sample", colnames(per.root.ind.match), fixed = F, invert = T)]
  # make list to hold results
  root.ind.dist.list <- list()
  root.ind.dist.list$p <- list()
  root.ind.dist.list$a <- list()
  # for each root sample
  for (r in 1:nrow(root.pos)){
    # make a list of present and absent individuals for each root sample
    r.name <- root.pos$sample[[r]]
    root.ind.dist.list$p[[r.name]] <- list()
    root.ind.dist.list$a[[r.name]] <- list()
    # get the coordinates of the root sample
    r.coords <- root.pos[r, c("x","y")]
    # get column numbers for the individuals that are present/absent
    p.n <- which(per.root.ind.match[r, ] > 0)
    a.n <- which(per.root.ind.match[r, ] == 0)
    # if there are any individuals present
    if(length(p.n) > 0){
      # get the names of the individuals that are present
      p.names <- do.call("rbind", strsplit(colnames(per.root.ind.match)[p.n], split = ".", fixed = T))[, 2]
      # for each present individual
      for(n in 1:length(p.n)){
        # get coordinates of the tree
        p.coords <- metadata[which(metadata$sample == p.names[[n]]), c("tree.pos.x", "tree.pos.y")]
        # get the distance
        p.dist <- proxy::dist(r.coords, p.coords, method = "Euclidean")
        # add to list
        root.ind.dist.list[["p"]][[r.name]][[p.names[[n]]]] <- min(p.dist)
      }
    }
    # flatten the sublist
    root.ind.dist.list[["p"]][[r.name]] <- unlist(root.ind.dist.list[["p"]][[r.name]])
    # get the names of the individuals that aren't present
    a.names <- do.call("rbind", strsplit(colnames(per.root.ind.match)[a.n], split = ".", fixed=T))[, 2]
    # for each absent individual
    for(n in 1:length(a.n)){
      # get coordinates of the tree
      a.coords <- metadata[which(metadata$sample == a.names[[n]]), c("tree.pos.x", "tree.pos.y")]
      # get the distance
      a.dist <- proxy::dist(r.coords, a.coords, method = "Euclidean")
      # add to list
      root.ind.dist.list[["a"]][[r.name]][[a.names[[n]]]] <- min(a.dist)
    }
    # flatten the sublist
    root.ind.dist.list[["a"]][[r.name]] <- unlist(root.ind.dist.list[["a"]][[r.name]])
  }
  root.ind.dist.list
}
# get distance lists (with and without a match depth filter)
in.root.dist.list.1 <- make.ind.dist.list(root.pos = root.pos, matches = best.matches, metadata = metadata)
in.root.dist.list.3 <- make.ind.dist.list(root.pos = root.pos, matches = best.matches.3.3, metadata = metadata)
# make boxplot of distance between root samples and present vs absent species
# min depth = 1
pdf("results/tree_vs_root_pos/root.tree.dist.individual.boxplot.md1.pdf")
boxplot(list("Individuals detected" = do.call("c", in.root.dist.list.1$p), "Individuals undetected" = do.call("c", in.root.dist.list.1$a)), ylab = "Root-to-tree distance (m)")
dev.off()
# min depth = 3
pdf("results/tree_vs_root_pos/root.tree.dist.individual.boxplot.md3.pdf")
boxplot(list("Individuals detected" = do.call("c", in.root.dist.list.3$p), "Individuals undetected" = do.call("c", in.root.dist.list.3$a)), ylab = "Root-to-tree distance (m)")
dev.off()
# test difference with Mann-Whitney U tests
# sink to text file
sink("results/tree_vs_root_pos/root.tree.dist.individual.tests.txt")
cat("Distance from root samples to present vs absent individuals, min read depth = 1:","\n")
print(wilcox.test(do.call("c", in.root.dist.list.1$p), do.call("c", in.root.dist.list.1$a), paired = F))
cat("Distance from root samples to present vs absent individuals, min read depth = 3:","\n")
print(wilcox.test(do.call("c", in.root.dist.list.3$p), do.call("c", in.root.dist.list.3$a), paired = F))
sink()


