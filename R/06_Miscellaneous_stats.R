# set seed
set.seed(1234)
# load data
load("results/processed_data/best.diag.RData")
load("results/processed_data/best.match.md0.RData")
load("results/processed_data/best.match.md3.RData")
metadata <- read.table("data/metadata.csv", sep = ",", header = T, na.strings = "NA", stringsAsFactors = F)
# make list of which individuals are from which species
ind.list <- list()
for(i in sort(unique(metadata$species))){
  ind.list[[i]] <- metadata[which(metadata$species == i),"sample"]
}
rm(i)
# make output directory
dir.create("results/misc_stats")

### Number of diagnostic haplotypes and markers identified
# table of N diagnostic haplotypes
all.hap.summary <- do.call(rbind, best.diag$summary$individual.diagnostic.haplotypes)
all.hap.summary$species <- gsub("~.+","",gsub(".","~",row.names(all.hap.summary),fixed=T))
rownames(all.hap.summary) <- NULL
# sink to text file
sink("results/misc_stats/range.n.diagnostic.txt")
cat("N species diagnostic markers (range):","\n")
cat(range(best.diag$summary$species.diagnostic.markers$N.diag.markers),"\n")
cat("N individual diagnostic haplotypes (range):","\n")
cat(range(all.hap.summary[which(all.hap.summary$individual != "non.diagnostic"), "N.diag.haplotypes"]),"\n")
sink()
# Table of N diagnostic markers per species 
tab.sp.diag <- best.diag$summary$species.diagnostic.markers
colnames(tab.sp.diag) <- c("Species", "N diagnostic markers")
write.table(tab.sp.diag, file = "results/misc_stats/N.diagnostic.markers.csv", sep = ",", row.names = F, quote = F)
# Table of N diagnostic haplotypes per individual
tab.ind.diag <- all.hap.summary[,c("species", "individual", "N.diag.haplotypes")]
colnames(tab.ind.diag) <- c("Species", "Individual", "N diagnostic haplotypes")
write.table(tab.ind.diag, file = "results/misc_stats/N.diagnostic.haplotypes.csv", sep = ",", row.names = F, quote = F)

### Matches per root sample
# Some of the info was already calculated in the 05_Rarefaction_analysis.R script, so load the output of that first
N.matches <- read.table("results/rarefaction/rarefaction.slopes.csv", sep=",", header=T, na.strings = "NA", stringsAsFactors = F)
# Concatenate diagnostic haplotype matches
per.root.ind.match.1 <- do.call("cbind",best.matches$individual.diag.haplotypes)
per.root.ind.match.3 <- do.call("cbind",best.matches.3.3$individual.diag.haplotypes)
# Remove sample and non-diagostic columns
per.root.ind.match.1 <- per.root.ind.match.1[,grep(".non.diagnostic|.sample",colnames(per.root.ind.match.1),fixed=F,invert = T)]
per.root.ind.match.3 <- per.root.ind.match.3[,grep(".non.diagnostic|.sample",colnames(per.root.ind.match.3),fixed=F,invert = T)]
# convert to counts (matches are reported as proportions of species diagnostic marker reads that have diagnostic haplotypes for each individual)
for(n in 1:ncol(per.root.ind.match.1)){
  sp <- strsplit(colnames(per.root.ind.match.1)[[n]],".",fixed=T)[[1]][1]
  per.root.ind.match.1[, n] <- per.root.ind.match.1[, n] * best.matches$species.diag.markers[, paste(sp, "N.reads", sep = ".")]
  per.root.ind.match.3[, n] <- per.root.ind.match.3[, n] * best.matches.3.3$species.diag.markers[, paste(sp, "N.reads", sep = ".")]
}
# get totals
per.root.ind.match.1$total <- rowSums(per.root.ind.match.1)
per.root.ind.match.3$total <- rowSums(per.root.ind.match.3)
# make table
match.counts <- data.frame(sample = best.matches$species.diag.markers$sample,
                           Total.reads = N.matches$n.reads,
                           Total.cat.match.reads = N.matches$n.cat.reads,
                           Total.spp.diag.reads = rowSums(best.matches$species.diag.markers[,2:ncol(best.matches$species.diag.markers)]),
                           N.ind.diag.reads = per.root.ind.match.1$total
)
match.counts$N.unmatched.reads <- match.counts$Total.reads - match.counts$Total.cat.match.reads
match.counts$N.nondiagnostic.cat.reads <- match.counts$Total.cat.match.reads - match.counts$Total.spp.diag.reads
match.counts$N.only.spp.diag.reads <- match.counts$Total.spp.diag.reads - match.counts$N.ind.diag.reads
match.counts <- match.counts[,c("sample", "Total.reads", "Total.cat.match.reads", "Total.spp.diag.reads", "N.unmatched.reads", "N.nondiagnostic.cat.reads", "N.only.spp.diag.reads", "N.ind.diag.reads")]
# save table
tab.match.counts <- match.counts
colnames(tab.match.counts) <- c("Sample", "N total reads", "N total catalogue-matching reads", "N total species-diagnostic reads", "N unmatched reads", "N nondiagnostic catalogue-matching reads", "N only species-diagnostic reads", "N individual-diagnostic reads")
write.csv(tab.match.counts, file = "results/misc_stats/root.read.match.counts.csv", row.names = F)
# sink to text file
sink("results/misc_stats/matches.per.root.sample.txt")
cat("N root reads matching catalogue (range):","\n")
cat(range(match.counts$Total.cat.match.reads), "\n")
cat("Percent of catalogue-matching reads that match species-diagnostic markers (range):","\n")
cat(round(range(match.counts$Total.spp.diag.reads / match.counts$Total.cat.match.reads * 100), 2), "\n")
cat("Percent of catalogue-matching reads that match individual-diagnostic haplotypes (range):","\n")
cat(round(range(match.counts$N.ind.diag.reads / match.counts$Total.cat.match.reads * 100), 2), "\n")
cat("N root samples with at least one matching species-diagnostic marker:","\n")
cat(length(which(match.counts$Total.spp.diag.reads > 0)), "\n")
cat("N root samples with at least one matching individual-diagnostic haplotype:","\n")
cat(length(which(match.counts$N.ind.diag.reads > 0)), "\n")
sink()

### Matches per species/individual
# number of root samples that each species was detected in
detected.spp.1 <- data.frame(species = sort(unique(metadata$species)), n.roots = NA, n.ind = NA)
detected.spp.3 <- detected.spp.1
for(n in 1:nrow(detected.spp.1)){
  detected.spp.1[n,"n.roots"] <- length(which(best.matches$species.diag.markers[,paste(detected.spp.1[n,"species"],".N.reads",sep="")] > 0))
  detected.spp.3[n,"n.roots"] <- length(which(best.matches.3.3$species.diag.markers[,paste(detected.spp.3[n,"species"],".N.reads",sep="")] > 0))
  n.ind <- length(ind.list[[detected.spp.1[n,"species"]]])
  detected.spp.1[n,"n.ind"] <- n.ind
  detected.spp.3[n,"n.ind"] <- n.ind
}
detected.spp.1$tree <- detected.spp.1$species %in% metadata[which(metadata$type == "tree"), "species"]
detected.spp.3$tree <- detected.spp.1$tree
# number of root samples that each individual was detected in
detected.ind.1 <- unname(colSums(per.root.ind.match.1 / per.root.ind.match.1, na.rm = T))
detected.ind.3 <- unname(colSums(per.root.ind.match.3 / per.root.ind.match.3, na.rm = T))
root.ind.names <- as.data.frame(do.call("rbind", strsplit(colnames(per.root.ind.match.1), split = ".", fixed = T))[, 1:2])
detected.ind.1 <- cbind(root.ind.names, detected.ind.1)
detected.ind.3 <- cbind(root.ind.names, detected.ind.3)
colnames(detected.ind.1) <- c("species", "ind", "n.roots")
colnames(detected.ind.3) <- c("species", "ind", "n.roots")
detected.ind.1 <- detected.ind.1[which(detected.ind.1$species != "total"), ]
detected.ind.3 <- detected.ind.3[which(detected.ind.3$species != "total"), ]
detected.ind.1$n.ind <- NA
for (n in 1:nrow(detected.ind.1)) {
  detected.ind.1[n, "n.ind"] <- length(ind.list[[detected.ind.1[n, "species"]]])
}
detected.ind.3$n.ind <- detected.ind.1$n.ind
# sink to text file
sink("results/misc_stats/matches.per.sp+ind.txt")
cat("N root samples per tree/shrub species, min read depth = 1 (range):","\n")
cat(range(detected.spp.1[which(detected.spp.1$tree == TRUE), "n.roots"]), "\n")
cat("N root samples per tree/shrub species, min read depth = 3 (range):","\n")
cat(range(detected.spp.3[which(detected.spp.3$tree == TRUE), "n.roots"]), "\n")
cat("N root samples per subshrub/herb species, min read depth = 1 (range):","\n")
cat(range(detected.spp.1[which(detected.spp.1$tree == FALSE), "n.roots"]), "\n")
cat("N root samples per subshrub/herb species, min read depth = 3 (range):","\n")
cat(range(detected.spp.3[which(detected.spp.3$tree == FALSE), "n.roots"]), "\n")
cat("N tree/shrub species detected, min read depth = 1:","\n")
cat(length(which(detected.spp.1$n.roots > 0 & detected.spp.1$tree == TRUE)),"\n")
cat("N tree/shrub species detected, min read depth = 3:","\n")
cat(length(which(detected.spp.3$n.roots > 0 & detected.spp.3$tree == TRUE)),"\n")
cat("N subshrub/herb species detected, min read depth = 1:","\n")
cat(length(which(detected.spp.1$n.roots > 0 & detected.spp.1$tree == FALSE)),"\n")
cat("N subshrub/herb species detected, min read depth = 3:","\n")
cat(length(which(detected.spp.3$n.roots > 0 & detected.spp.3$tree == FALSE)),"\n")
cat("Median N root samples per individual, min read depth = 1:","\n")
cat(median(detected.ind.1$n.roots), "\n")
cat("Median N root samples per individual, min read depth = 3:","\n")
cat(median(detected.ind.3$n.roots), "\n")
cat("N individuals detected, min read depth = 1:","\n")
cat(length(which(detected.ind.1$n.roots > 0)),"\n")
cat("N individuals detected, min read depth = 3:","\n")
cat(length(which(detected.ind.3$n.roots > 0)),"\n")
sink()

### Effect of N individuals per species
# sink to text file
sink("results/misc_stats/Nind.cor.txt")
cat("N individuals per species vs N species-diagnostic markers correlation","\n")
print(cor.test(best.diag$summary$species.diagnostic.markers$N.diag.markers, lengths(ind.list[best.diag$summary$species.diagnostic.markers$species]), method = "s", exact = F))
cat("N individuals per species vs N individual-diagnostic haplotypes correlation","\n")
print(cor.test(all.hap.summary[which(all.hap.summary$individual != "non.diagnostic"),"N.diag.haplotypes"], lengths(ind.list[all.hap.summary[which(all.hap.summary$individual != "non.diagnostic"),"species"]]), method = "s", exact = F))
cat("N individuals per species vs N species matches correlation, min read depth = 1","\n")
print(cor.test(detected.spp.1$n.ind, detected.spp.1$n.roots, method = "s", exact = F))
cat("N individuals per species vs N species matches correlation, min read depth = 3","\n")
print(cor.test(detected.spp.3$n.ind, detected.spp.3$n.roots, method = "s", exact = F))
cat("N individuals per species vs N individual matches correlation, min read depth = 1","\n")
print(cor.test(detected.ind.1$n.ind, detected.ind.1$n.roots, method = "s", exact = F))
cat("N individuals per species vs N individual matches correlation, min read depth = 3","\n")
print(cor.test(detected.ind.3$n.ind, detected.ind.3$n.roots, method = "s", exact = F))
sink()
