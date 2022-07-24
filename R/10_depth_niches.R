library(coin)
library(gtools)
library(RColorBrewer)
# load data
load("results/processed_data/best.match.md0.RData")
load("results/processed_data/best.match.md3.RData")
metadata <- read.table("data/metadata.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
root.pos <- read.table("results/tree_vs_root_pos/root.pos.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
# create output directory
dir.create("results/depth_distribution/")
# define test root depth bias function
test_depth_bias <- function(matches, root.pos){
  # get depths
  depths <- factor(sort(unique(root.pos$depth)), ordered = T)
  # species 
  species <- gsub(".N.reads","",colnames(matches$species.diag.markers[-1]),fixed = T)
  # dataframe
  my.df <- data.frame(matrix(NA,ncol=length(depths)+4,nrow = length(species),dimnames = list(NULL,c("species",paste("N.depth.",depths,sep=""),"lbl.stat","lbl.P","lbl.Padj"))))
  # for each species
  for(n in 1:length(species)){
    # add species name
    my.df[n,"species"] <- species[[n]]
    # initiate matrix
    my.mat <- matrix(NA,nrow=length(depths),ncol=2,dimnames = list(depths,c("pres","abs")))
    # get matches
    my.matches <- matches$species.diag.markers[which(matches$species.diag.markers[,paste(species[[n]],"N.reads",sep=".")]>0),"sample"]
    # get depths
    for(d in depths){
      my.mat[d,"pres"] <- length(which(root.pos[which(root.pos$sample %in% my.matches),"depth"] == d))
      my.mat[d,"abs"] <- length(which(root.pos$depth == d)) - my.mat[d,"pres"]
    }
    # add depths to data.frame
    my.df[n,2:(length(depths)+1)] <- my.mat[,"pres"]
    # do lbl test
    my.lbl <- coin::lbl_test(as.table(my.mat))
    my.df[n,"lbl.stat"] <- my.lbl@statistic@teststatistic
    my.df[n,"lbl.P"] <- pvalue(my.lbl)
  }
  # multiple test correction
  my.df$lbl.Padj <- p.adjust(my.df$lbl.P,method = "fdr")
  return(my.df)
}
# get depth bias
my.depth.bias <- test_depth_bias(matches = best.matches,root.pos = root.pos)
my.depth.bias.3 <- test_depth_bias(matches = best.matches.3.3,root.pos = root.pos)
# list of which are trees
tree_shrub <- metadata[match(my.depth.bias$species, metadata$species), "type"]
# remove underscore from species names
my.names <- gsub("_", " ", my.depth.bias$species, fixed = T)
# calculate spaces for adding text
my.spaces <- 0.2 * 1:length(my.names)
# make pdf
pdf("./results/depth_distribution/Depth.barplot.md1.pdf")
par(mar=c(15,6,4.1,2.1))
# bar plot
barplot(t(as.matrix(my.depth.bias[, 2:5])),
        beside = F,
        col = brewer.pal(4, name = "YlGnBu"),
        names.arg = rep("", length(my.names)),
        las = 2,
        ylab = "N root samples",
        cex.names = 1.2,
        cex.axis = 1.2,
        font = 1,
        col.sub = "blue",
        font.axis = 1,
        cex.lab = 1.5)
# stars for P vals
text(1:19 * 1.2 - 0.5, rowSums(my.depth.bias[, 2:5]) + 3, labels = stars.pval(my.depth.bias$lbl.Padj), cex = 2)
# x axis with bold for trees
axis(1, at = which(tree_shrub == "tree") + my.spaces[which(tree_shrub == "tree")] - 0.45, labels = my.names[which(tree_shrub == "tree")], las = 2, font = 4)
# x axis with non-bold for non-trees
axis(1, at = which(tree_shrub == "shrub") + my.spaces[which(tree_shrub == "shrub")] - 0.45, labels = my.names[which(tree_shrub == "shrub")], las = 2, font = 3)
# legend
legend("topleft", legend = c("0 - 5 cm", "5 - 10 cm", "15 - 20 cm", "45 - 50 cm"), fill = brewer.pal(4, name = "YlGnBu"), bty = "n", cex = 0.8, title = "Depth")
dev.off()
# save stat tables
write.csv(my.depth.bias, file = "results/depth_distribution/depth.bias.tests.md1.csv", row.names = F, quote = F)
write.csv(my.depth.bias.3, file = "results/depth_distribution/depth.bias.tests.md3.csv", row.names = F, quote = F)

