library(viridis)
library(ppcor)
# load data
load("results/processed_data/best.match.md0.RData")
load("results/processed_data/best.match.md3.RData")
metadata <- read.table("data/metadata.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
root.pos <- read.table("results/tree_vs_root_pos/root.pos.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
# make list of which individuals are from which species
ind.list <- list()
for(i in sort(unique(metadata$species))){
  ind.list[[i]] <- metadata[which(metadata$species == i),"sample"]
}
rm(i)
# make results dir
dir.create("results/ag_bg_correlations/")
# add N.ind per species
metadata$N.ind <- apply(metadata, 1, function(x) length(which(metadata$species == x["species"] )))
# add canopy lateral spread (maximum canopy radius)
metadata$canopy.LS <- apply(metadata[, c("canopy_diameter_NS","canopy_diameter_EW")], 1, function(x) max(x)/2)
# add canopy area
metadata$canopy.area <- (metadata$canopy_diameter_NS / 2) * (metadata$canopy_diameter_EW / 2) * pi
# add canopy volume
metadata$canopy.vol <- 4 / 3 * pi * (metadata$canopy_diameter_NS/2) * (metadata$canopy_diameter_EW/2) * ((metadata$h_tree - metadata$h_crown)/2)

# function to get root lateral spread (i.e. root radius) for each individual
get_root_LS <- function(matches, root.pos, metadata, ind.list, min.marker.n = 1, min.hap.prop = .Machine$double.xmin){
  # get species names
  multi.ind <- names(which(lengths(ind.list) > 1))
  single.ind <- setdiff(names(ind.list), multi.ind)
  # initiate data frame
  my.df <- data.frame(matrix(NA, ncol=3, nrow = sum(lengths(ind.list)), dimnames = list(NULL, c("species", "ind", "est.root.radius"))))
  # fill species and ind columns
  my.df$species <- gsub("[0-9]", "", names(unlist(ind.list)))
  my.df$ind <- unlist(ind.list)
  # get matches
  for (n in 1:nrow(my.df)){
    my.sp <- my.df[n,"species"]
    my.ind <- my.df[n, "ind"]
    if(my.sp %in% single.ind){
      my.matches <- matches$species.diag.markers[which(matches$species.diag.markers[, paste(my.sp, "N.reads", sep = ".")] >= min.marker.n), "sample"]
    } else {
      my.matches <- matches$individual.diag.haplotypes[[my.sp]][which(matches$individual.diag.haplotypes[[my.sp]][ ,paste(my.ind, "N.reads", sep = ".")] >= min.hap.prop), "sample"]
    }
    if(length(my.matches) > 0){
      # get tree pos
      tree.coords <- metadata[which(metadata$sample == my.ind), c("tree.pos.x", "tree.pos.y")]
      match.coords <- root.pos[which(root.pos$sample %in% my.matches), c("x", "y")]
      my.df[n,"est.root.radius"] <- max(proxy::dist(tree.coords, match.coords, method = "Euclidean"))
      }
    my.df[which(my.df$est.root.radius == -Inf),"est.root.radius"] <- NA
    
  }
  my.df
}
# get root lateral spread
root.LS <- get_root_LS(matches = best.matches, root.pos = root.pos, metadata = metadata, ind.list = ind.list)
root.LS.3 <- get_root_LS(matches = best.matches.3.3, root.pos = root.pos, metadata = metadata, ind.list = ind.list)
metadata$root.LS <- root.LS[match(x = metadata$sample, table = root.LS$ind),"est.root.radius"]
metadata$root.LS.3 <- root.LS.3[match(x = metadata$sample, table = root.LS.3$ind),"est.root.radius"]
# function to get the number of root matches per individual
get_N_root_match <-  function(matches, ind.list, min.marker.n=NA, min.hap.prop = NA){
  # get species names
  multi.ind <- names(which(lengths(ind.list) > 1))
  single.ind <- setdiff(names(ind.list), multi.ind)
  # initiate data frame
  my.df <- data.frame(matrix(NA, ncol = 3, nrow = sum(lengths(ind.list[multi.ind])) + length(single.ind), dimnames = list(NULL, c("species","ind","N.root.matches"))))
  # fill species and ind columns
  spp.names <- list()
  ind.names <- list()
  for(i in names(ind.list)){
    spp.names[[i]] <- rep(i, length(ind.list[[i]]))
    ind.names[[i]] <- ind.list[[i]]
  }
  my.df$species <- unname(unlist(spp.names))
  my.df$ind <- unname(unlist(ind.names))
  # get N root matches
  for (n in 1:nrow(my.df)){
    my.sp <- my.df[n, "species"]
    if(my.sp %in% single.ind){
      if(is.na(min.marker.n)){
        my.df[n, "N.root.matches"] <- length(which(matches$species.diag.markers[, paste(my.sp, "N.reads", sep = ".")] > 0))
      } else {
        my.df[n, "N.root.matches"] <- length(which(matches$species.diag.markers[, paste(my.sp, "N.reads", sep = ".")] >= min.marker.n))
      }
    } else {
      my.ind <- my.df[n, "ind"]
      if(is.na(min.hap.prop)){
        my.df[n, "N.root.matches"] <- length(which(matches$individual.diag.haplotypes[[my.sp]][, paste(my.ind, "N.reads", sep = ".")] > 0))
      } else {
        my.df[n, "N.root.matches"] <- length(which(matches$individual.diag.haplotypes[[my.sp]][, paste(my.ind, "N.reads", sep = ".")] >= min.hap.prop))
      }
    }
  }
  return(my.df)
}
# get N root matches
N.root <- get_N_root_match(matches = best.matches, ind.list = ind.list)
N.root.3 <- get_N_root_match(matches = best.matches.3.3, ind.list = ind.list)
metadata$N.root <- N.root[match(x = metadata$sample, table = N.root$ind),"N.root.matches"]
metadata$N.root.3 <- N.root.3[match(x = metadata$sample, table = N.root.3$ind),"N.root.matches"]

# function to get partial correlations between each aboveground metric and each belowground mtric controlling for a counfounding variable (z), P-values and data with NA removed between aboveground and belowground metrics
get_ag_bg_pcor <- function(ag.metrics, bg.metrics, z, metadata){
  # set up output
  cor.mat <- matrix(NA, ncol = length(ag.metrics), nrow = length(bg.metrics), dimnames = list(bg.metrics, ag.metrics))
  p.mat <- cor.mat
  reduced.dat <- list()
  cor.tests <- list()
  # for each combination
  for(b in bg.metrics){
    for(a in ag.metrics){
      # get reduced data
      dat <- metadata[, c("species",a, b, z)]
      dat <- dat[complete.cases(dat),]
      # get corr
      my.cor <- ppcor::pcor.test(x = dat[, a], y = dat[, b], z = dat[, z], method = "spearman")
      # add corr to matrix
      cor.mat[b, a] <- my.cor$estimate
      # add P value to matrix
      p.mat[b, a] <- my.cor$p.value
      reduced.dat[[paste(a,b,sep="_vs_")]] <- dat
      cor.tests[[paste(a,b,sep="_vs_")]] <- my.cor
    }
  }
  list(rho = cor.mat, p = p.mat, cor = cor.tests, data = reduced.dat)
}
# get corrlations
ag.metrics = c("h_tree", "canopy.LS", "canopy.area", "canopy.vol", "h_crown")
bg.metrics = c("root.LS", "N.root")
bg.metrics.3 = c("root.LS.3", "N.root.3")
ag_bg_pcor <- get_ag_bg_pcor(ag.metrics = ag.metrics, bg.metrics = bg.metrics, z = "N.ind", metadata = metadata)
ag_bg_pcor.3 <- get_ag_bg_pcor(ag.metrics = ag.metrics, bg.metrics = bg.metrics.3, z = "N.ind", metadata = metadata)
# print correlations to files
# as table
out.tab <- data.frame(ag = rep(colnames(ag_bg_pcor$rho), each = 2), bg = rep(rownames(ag_bg_pcor$rho),5), rho.1 = as.vector(ag_bg_pcor$rho), p.1 = as.vector(ag_bg_pcor$p), rho.3 = as.vector(ag_bg_pcor.3$rho), p.3 = as.vector(ag_bg_pcor.3$p))
write.csv(out.tab, file = "results/ag_bg_correlations/ppcor.tab.csv", quote = F)
# as text
sink("results/ag_bg_correlations/partial.cor.txt")
cat("Min match depth = 1", "\n\n")
for(n in names(ag_bg_pcor$cor)){
  cat(n, "\n")
  print(ag_bg_pcor$cor[[n]])
}
cat("\n\n", "Min match depth = 3", "\n\n")
for(n in names(ag_bg_pcor.3$cor)){
  cat(n, "\n")
  print(ag_bg_pcor.3$cor[[n]])
}
sink()
# plot and output correlations
# nice names
nice.names <- list(h_tree = 'Tree height (m)',
                   h_crown = 'Crown base height (m)',
                   canopy.LS = 'Canopy radius (m)',
                   canopy.area = expression('Canopy area (m'^2~')'),
                   canopy.vol = expression('Canopy volume (m'^3~')'),
                   root.LS = 'Root radius (m)',
                   root.LS.3 = 'Root radius (m)',
                   N.root = 'N root matches',
                   N.root.3 = 'N root matches')
# plot function
plot_ag.bg <- function(ag.metrics, bg.metrics, ag_bg_pcor, nice.names, prefix, cols){
  for(b in bg.metrics){
    for(a in ag.metrics){
      pdf(paste(prefix, paste(a,b,sep="_vs_"), "pdf", sep = "."),paper = "a4r")
      par(mfrow=c(2,2),mar=c(5.1, 5.1, 4.1, 2.1))
      my.cols <- cols[which(names(cols) %in% ag_bg_pcor$data[[paste(a, b, sep = "_vs_")]][, "species"])]
      plot(ag_bg_pcor$data[[paste(a, b, sep = "_vs_")]][, a], ag_bg_pcor$data[[paste(a, b, sep = "_vs_")]][,b], xlab = nice.names[[a]], ylab = nice.names[[b]], pch = 21, cex = 1.2, bg = my.cols[ag_bg_pcor$data[[paste(a, b, sep = "_vs_")]][, "species"]], cex.lab = 1.3)
      abline( lm(ag_bg_pcor$data[[paste(a, b, sep = "_vs_")]][,b] ~ ag_bg_pcor$data[[paste(a, b, sep = "_vs_")]][, a]))
      par(mar=c(0,0,0,2.1))
      plot.new()
      legend("left", legend = gsub("_"," ", names(my.cols)), pt.bg = my.cols, text.font = 3, pch = 21, cex = 1, bty = "n", ncol = 1, pt.cex = 2)
      dev.off()
    }
  }
}
# get cols
my.cols = viridis_pal(option = "C")(length(unique(ag_bg_pcor$data$h_tree_vs_N.root$species)))
names(my.cols) <- sort(unique(unique(ag_bg_pcor$data$h_tree_vs_N.root$species)))
#plot
plot_ag.bg(ag.metrics = ag.metrics, bg.metrics = bg.metrics, ag_bg_pcor = ag_bg_pcor, nice.names = nice.names, prefix = "./results/ag_bg_correlations/cor.md1", cols = my.cols)
plot_ag.bg(ag.metrics = ag.metrics, bg.metrics = bg.metrics.3, ag_bg_pcor = ag_bg_pcor.3, nice.names = nice.names, prefix = "./results/ag_bg_correlations/cor.md3", cols = my.cols)

