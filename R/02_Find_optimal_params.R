# load metadata
metadata <- read.table("data/metadata.csv", sep=",", header=T, na.strings = "NA", stringsAsFactors = F)
# load RootID runs
load("results/processed_data/stacks.dat.list.RData")
load("results/processed_data/diag.list.RData")
load("results/processed_data/match.list.RData")
# make output directory
dir.create("results/parameter_comparison")
# get names of species
spp.names <- sort(unique(metadata$species)) 
# get names of multi-individual species
Mspp.names <- names(table(metadata$species)[table(metadata$species)>1])
# get names of individuals from multi-individual species
ind.names <- sort(paste(metadata[which(metadata$species %in% Mspp.names),"species"],metadata[which(metadata$species %in% Mspp.names),"sample"],sep="_"))
# all parameter values to test
US.M.vals <- seq(2,8,2)
CS.n.vals <- c(-1,0,1)
maxmiss.vals <- cbind(c(0,0.2,0.2,0.2,0.4,0.4,0.4),c(0,0,0.1,0.2,0,0.1,0.2))
mindepth.vals <-  c(5,10,20)
maxhaps.vals <- c(2,3,NA)
# initiate data frame to hold numbers of diagnostic markers/haplotypes
n.diag <- data.frame(matrix(NA,
                            nrow = (length(US.M.vals) * length(CS.n.vals) * nrow(maxmiss.vals) * length(mindepth.vals) * length(maxhaps.vals)),
                            ncol =length(ind.names) + 6 + length(spp.names),
                            dimnames = list(NULL, c("US_M", "CS_n", "FD_max.md.marker", "FD_max.md.hap", "FD_min.dep", "FD_max.haps", ind.names, spp.names))))
# get number of diagnostic haplotypes and markers per individual or species
N=1
for(M in US.M.vals){
  for(n in CS.n.vals){
    for(md in 1:nrow(maxmiss.vals)){
      for(m in mindepth.vals){
        for(mh in maxhaps.vals){
          n.diag[N,1] <- M
          n.diag[N,2] <- as.integer(M+n)
          n.diag[N,3:4] <- c(maxmiss.vals[md,1],maxmiss.vals[md,2])
          n.diag[N,5] <- m
          n.diag[N,6] <- mh
          for(i in Mspp.names){
            for (j in metadata[which(metadata$species == i),"sample"]){
              n.diag[N,paste(i,j,sep="_")] <- length(diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]][[paste("FD_max.haps_",as.character(mh),sep="")]][["individual.diag.haplotypes"]][[i]][[j]])
            }
          }
          for(s in spp.names){
            n.diag[N,s] <- length(diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]][[paste("FD_max.haps_",as.character(mh),sep="")]][["species.diag.markers"]][[s]])
          }
          N = N+1
        }
      }
    }
  }
}
# rank rows for each column
rank.n.diag <- n.diag
for(i in 7:ncol(n.diag)){
  rank.n.diag[,i] <- rank(rank.n.diag[,i],ties.method = "min")
}
# get mean rank
rank.n.diag$mean.rank <- rowMeans(rank.n.diag[,7:ncol(n.diag)])
# choose column with best mean rank
optimal <- rank.n.diag[which.max(rank.n.diag$mean.rank),1:6]
optimal
#####    US_M CS_n FD_max.md.marker FD_max.md.hap FD_min.dep FD_max.haps
#####     6    7       0.4              0.2          5          NA
# add mean rank to n.diag
n.diag$mean.rank <- max(rank.n.diag$mean.rank)+1-rank.n.diag$mean.rank
# get best
best.stacks.dat <- stacks.dat.list$US_M_6$CS_n_7
best.diag <- diag.list$US_M_6$CS_n_7$FD_max.md.marker_0.4$FD_max.md.hap_0.2$FD_min.dep_5$FD_max.haps_NA
best.matches <- match.list$US_M_6$CS_n_7$FD_max.md.marker_0.4$FD_max.md.hap_0.2$FD_min.dep_5$FD_max.haps_NA
# save results
save(best.stacks.dat, file = "results/processed_data/best.stacks.dat.RData")
save(best.diag, file = "results/processed_data/best.diag.RData")
save(best.matches, file = "results/processed_data/best.match.md0.RData")
write.table(optimal,file="results/parameter_comparison/optimal_params.csv",sep=",",row.names = F,quote = F)
write.table(n.diag,file="results/parameter_comparison/N.diag.across.params.csv",sep=",",row.names = F,quote = F)
