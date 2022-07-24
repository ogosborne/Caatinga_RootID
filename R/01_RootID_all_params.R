library(RootID)
# load metadata
metadata <- read.table("data/metadata.csv",sep=",",header=T,na.strings = "NA",stringsAsFactors = F)
# make list of which individuals are from which species
ind.list <- list()
for(i in sort(unique(metadata$species))){
  ind.list[[i]] <- metadata[which(metadata$species == i),"sample"]
}
rm(i)
## lists to hold outputs
stacks.dat.list <- list()
diag.list <- list()
match.list <- list()
# parameter values to test
US.M.vals <- seq(2,8,2)
CS.n.vals <- c(-1,0,1)
maxmiss.vals <- cbind(c(0,0.2,0.2,0.2,0.4,0.4,0.4),c(0,0,0.1,0.2,0,0.1,0.2))
mindepth.vals <-  c(5,10,20)
maxhaps.vals <- c(2,3,NA)
# initiate data frame to record timing
timing <- data.frame(matrix(NA,nrow = (length(US.M.vals) * length(CS.n.vals) * nrow(maxmiss.vals) * length(mindepth.vals) * length(maxhaps.vals)),ncol = 10,dimnames = list(NULL,c("US_M","CS_n","FD_max.md.marker","FD_max.md.hap","FD_min.dep","FD_max.haps","time_RS","time_FD","time_MD","time_total"))))
# make output directories
dir.create("results")
dir.create("results/processed_data")
dir.create("results/timing")
#### run RootID on all parameter combinations, save timings
row.num <- 1
# USTACKS M
for(M in US.M.vals){
  # msg
  cat("US M =",M,"\n")
  # lists
  diag.list[[paste("US_M_",M,sep="")]] <- list()
  match.list[[paste("US_M_",M,sep="")]] <- list()
  stacks.dat.list[[paste("US_M_",M,sep="")]] <- list()
  # CSTACKS n
  for (n in CS.n.vals){
    # lists
    diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]] <- list()
    match.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]] <- list()
    # msg
    cat("CS n =",M+n,"\n")
    # timing
    time1 <- Sys.time()
    # read.stacks
    stacks.dat.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]] <- read.stacks(data.dir=paste("STACKS_data/sstacks/US_M",M,"/CS_n",M+n,"/leaf",sep=""),
                                                                                           min.dep=5,
                                                                                           verbose=F)
    # timing
    time2 <- Sys.time()
    td12 <- difftime(time2,time1,units = "secs")
    # lists
    for(i in unique(maxmiss.vals[,1])){
      diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",i,sep = "")]] <- list()
      match.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",i,sep = "")]] <- list()
    }
    # max missing data
    for(md in 1:nrow(maxmiss.vals)){
      # lists
      diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]] <- list()
      match.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]] <- list()
      # min depth 
      for(m in mindepth.vals){
        #lists
        diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]] <- list()
        match.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]] <- list()
        # max haplotypes
        for(mh in maxhaps.vals){
          # timing
          time3 <- Sys.time()
          # find.diag
          diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]][[paste("FD_max.haps_",as.character(mh),sep="")]] <- find.diag(ind.list = ind.list,
                                                                                                                                                                                                                                                                                            stacks.dat = stacks.dat.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]],
                                                                                                                                                                                                                                                                                            max.md.marker = maxmiss.vals[md,1],
                                                                                                                                                                                                                                                                                            max.md.hap = maxmiss.vals[md,2],
                                                                                                                                                                                                                                                                                            min.dep = m,
                                                                                                                                                                                                                                                                                            max.haps = mh)
          # timing
          time4 <- Sys.time()
          td34 <- difftime(time4,time3,units = "secs")
          time5 <- Sys.time()
          # match.diag
          match.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]][[paste("FD_max.haps_",as.character(mh),sep="")]] <- 
            match.diag(data.dir = paste("STACKS_data/sstacks/US_M",M,"/CS_n",M+n,"/root",sep=""),
                       diag = diag.list[[paste("US_M_",M,sep="")]][[paste("CS_n_",M+n,sep="")]][[paste("FD_max.md.marker_",maxmiss.vals[md,1],sep = "")]][[paste("FD_max.md.hap_",maxmiss.vals[md,2],sep = "")]][[paste("FD_min.dep_",m,sep="")]][[paste("FD_max.haps_",as.character(mh),sep="")]])
          # timing
          time6 <- Sys.time()
          td56 <- difftime(time6,time5,units = "secs")
          # write timing
          timing[row.num,1] <- M
          timing[row.num,2] <- M+n
          timing[row.num,3:6] <- c(maxmiss.vals[md,1],maxmiss.vals[md,2],m,mh)
          timing[row.num,7:9] <- c(td12,td34,td56)
          timing[row.num,10] <- sum(as.numeric(timing[row.num,7:9]))
          # iterator
          row.num <- row.num +1
        }
      }
    }
  }
}
# max and min timing 
sink("results/timing/timing.range.txt")
cat("Timing range", "\n", ":")
cat(range(timing$time_total))
sink()
# save results
save(stacks.dat.list,file = "results/processed_data/stacks.dat.list.RData")
save(diag.list,file = "results/processed_data/diag.list.RData")
save(match.list,file = "results/processed_data/match.list.RData")
write.table(timing,file="results/timing/timing.csv",sep=",",row.names = F,quote = F)

