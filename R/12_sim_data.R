library(RootID)
library(wesanderson)

# make list of which individuals are from which species
ind.lists <- list()
spp = c("Ao", "At", "Eg", "Ha", "Ls", "Pt")
nreads <- as.integer(seq(100000, 1000000, 100000))
for(n in nreads){
  ind.lists[[paste("N", n, sep = "")]]
  for(i in spp){
    ind.lists[[paste("N", n, sep = "")]][[i]] <- paste(i, 0:9, n, sep = "_" )
  }
}
rm(list = c("i","n"))
# genome sizes
genome.size <- c(Ao = 1187539004, At = 119667750, Eg = 691346258, Ha = 3027963057, Ls = 2399509484, Pt = 392162179)
# full names
full.names <- c(Ao = "Asparagus officinalis", At  = "Arabidopsis thaliana", Eg = "Eucalyptus grandis", Ha = "Helianthus annuus", Ls = "Lactuca sativa", Pt = "Populus trichocarpa")
# read stacks output
indir <- "SIM_data/sstacks/"
stacks.dat.list <- list()
diag.list <- list()
for(n in nreads){
  stacks.dat.list[[paste("N", n, sep = "")]] <- read.stacks(data.dir=paste(indir, "N", n, sep=""), min.dep = 5, verbose = T)
}
# find diagnostic markers and haplotypes
diag.list <- list()
for(n in nreads){
  diag.list[[paste("N", n, sep = "")]] <- find.diag(ind.list = ind.lists[[paste("N", n, sep = "")]],
                                                    stacks.dat = stacks.dat.list[[paste("N", n, sep = "")]],
                                                    max.md.marker = 0.4,
                                                    max.md.hap = 0.2,
                                                    min.dep = 5,
                                                    max.haps = NA)
}
# get N species diagnostic markers
N.sp.mar <- data.frame(matrix(NA, ncol = 7, nrow = 10, dimnames = list(NULL, c("N.reads", spp))))
N.sp.mar$N.reads <- nreads
for(n in 1:length(nreads)){
  N.sp.mar[n,2:7] <- diag.list[[paste("N", nreads[[n]], sep = "")]]$summary$species.diagnostic.markers$N.diag.markers
}
# get mean N individual-diagnostic haplotypes 
N.ind.hap <- list()
N.ind.hap$mean <- data.frame(matrix(NA, ncol = 7, nrow = 10, dimnames = list(NULL, c("N.reads", spp))))
N.ind.hap$mean$N.reads <- nreads
for(sp in spp){
  for(n in 1:10){
    N.ind.hap$mean[n,sp] <-  mean(diag.list[[paste("N", nreads[[n]], sep = "")]]$summary$individual.diagnostic.haplotypes[[sp]][1:10,"N.diag.haplotypes"])
  }
}
# get N individual-diagnostic haplotypes per individual
for(sp in spp){
  N.ind.hap[[sp]] <- data.frame(matrix(NA, ncol = 11, nrow = 10, dimnames = list(NULL, c("N.reads", paste(sp, 0:9, sep = "_")))))
  N.ind.hap[[sp]]$N.reads <- nreads
  for(ind in 0:9){
    for(n in 1:10){
      N.ind.hap[[sp]][n, paste(sp, ind, sep = "_")] <- diag.list[[paste("N", nreads[[n]], sep = "")]]$summary$individual.diagnostic.haplotypes[[sp]][(ind+1), "N.diag.haplotypes"]
    }
  }
}
#### plot
# create output directory
dir.create("results/sim_data")
## layout
pdf("results/sim_data/Nreads.vs.genome.size.pdf")
layout(matrix(c(1,2,5,3,4,5), ncol = 2))
cols <- wes_palette("Zissou1", 6, type = "continuous")
names(cols) <- names(sort(genome.size))
## plot 1: N reads vs N species diagnostic markers
plot(N.sp.mar$N.reads, N.sp.mar$Ao, type = "b", ylim = c(0, max(N.sp.mar[,2:7])), col = cols[["Ao"]], lwd = 4, xlab = "N reads", ylab = "N species diagnostic markers", main = "(a)")
for(sp in spp[2:6]){
  points(N.sp.mar$N.reads, N.sp.mar[[sp]], type = "b", col = cols[[sp]], lwd = 4)
}
## plot 2: N reads vs N individual diagnostic haplotypes
plot(N.ind.hap$mean$N.reads, N.ind.hap$mean$Ao, type = "b", ylim = c(0, 490), col = cols[["Ao"]], lwd = 4, xlab = "N reads", ylab = "N individual diagnostic haplotypes", main = "(c)")
# plot per individual
for(sp in spp){
  for(ind in 0:9){
    points(N.ind.hap$mean$N.reads, N.ind.hap[[sp]][,paste(sp, ind, sep = "_")], type = "b", col = cols[[sp]], lwd = 1, lty = 2)
  }
}
# plot means
for(sp in spp){
  points(N.ind.hap$mean$N.reads, N.ind.hap$mean[[sp]], type = "b", col = cols[[sp]], lwd = 4)
}
## plot 3: Genome size vs N species diagnostic markers
plot(genome.size/1000000, N.sp.mar[10,2:7], col = cols[names(genome.size)], xlab = "Genome size (Mbp)", ylab = "N species diagnostic markers", pch = 19, cex = 2, main = "(b)")
## plot 4: Genome size vs N individual diagnostic haplotypes
all.ind <- do.call(cbind, N.ind.hap[c(spp)])
all.ind <- all.ind[,grep("N.reads", colnames(all.ind), invert = T)]
sp.vec <- gsub("\\..*", "", colnames(all.ind))
df <- data.frame(gs = genome.size[sp.vec]/1000000, nh = unlist(all.ind[10,]))
plot(df$gs, df$nh, col = cols[sp.vec], xlab = "Genome size (Mbp)", ylab = "N individual diagnostic haplotypes", pch = 1, cex = 2, main = "(d)")
abline(lm(nh ~ gs, data = df)) 
# legend
plot.new()
legend("top", legend = full.names[names(sort(genome.size))], fill = cols[names(sort(genome.size))])
dev.off()
# test
sink("results/sim_data/Genome_size.cor.txt")
cat("Genome size vs N species diagnostic markers","\n")
print(cor.test(genome.size, unlist(N.sp.mar[10,2:7]), method = "s", exact = F))
cat("Genome size vs mean N individual diagnostic haplotypes","\n")
print(cor.test(genome.size, unlist(N.ind.hap$mean[10,2:7]), method = "s", exact = F))
sink()
