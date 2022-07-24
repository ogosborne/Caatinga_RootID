library(RootID)
# make output directory
dir.create("results/rarefaction/")
# load best diagnostic markers
load("results/processed_data/best.diag.RData")
# Number of subsamples
N=50
# run rarefaction
set.seed(1234)
rarefaction.dat <- sample.rarefaction(data.dir = "STACKS_data/US_M6/CS_n7/root/",
                                      n.subsamples = N,
                                      n.reps = 100,
                                      diag = best.diag,
                                      verbose = T)
# save result
save(rarefaction.dat,file="results/rarefaction/rarefaction.dat.RData")
# initialise data frame to hold slopes
slopes <- data.frame(sample = names(rarefaction.dat$species$mean), m.markers = NA, m.haplotypes = NA, m.species = NA, m.individuals = NA, n.markers=NA,n.haplotypes=NA,n.species=NA,n.individuals=NA,n.cat.reads=NA)
# add N reads per sample
nreads <- read.table("data/nreads.tsv",header=T,sep=" ",stringsAsFactors = F)
slopes$n.reads <- nreads[match(slopes$sample, nreads$sample),"N.reads"]
rm(nreads)
# proportion of the final curve to calculate slope for
pslope=0.1
# calculate final slope
for(n in 1:nrow(slopes)){
  for (i in c("species","individuals","haplotypes","markers")){
    # slope 
    N90 <- rarefaction.dat[[i]]$mean[(round(N-pslope*N)),slopes[n,"sample"]]
    N100 <- rarefaction.dat[[i]]$mean[N,slopes[n,"sample"]]
    P90 <-  N90 / N100
    if(is.na(P90)) P90 <- 1
    m <- (1 - P90) / pslope
    slopes[n,paste("m",i,sep=".")] <- m
    # N
    slopes[n,paste("n",i,sep=".")] <- rarefaction.dat[[i]]$mean[N,slopes[n,"sample"]]
  }
  # n.cat reads
  slopes[n,"n.cat.reads"] <- rarefaction.dat$n.reads[[slopes[n,"sample"]]][[N]]
}
# save slopes table
write.table(slopes,file="results/rarefaction/rarefaction.slopes.csv",sep=",",row.names = F,quote = F)
# plot N reads vs slope 
pdf("results/rarefaction/n.read.slope.correlations.pdf")
par(mfrow=c(2,2),mar=c(4,4,2,1))
for (i in c("markers","haplotypes","species","individuals")){
  plot(slopes$n.reads, slopes[,paste("m", i, sep = ".")], xlab = expression(italic("N")~"reads"), ylab = bquote(italic("m")[.(i)]))
}
dev.off()
# sink stats to text file
sink("results/rarefaction/slope.stats.txt")
# cor tests
cat(" Spearman's correlation tests","\n")
for (i in c("markers","haplotypes","species","individuals")){
  cat(i, ":","\n")
  print(cor.test(slopes$n.reads, slopes[, paste("m", i, sep = ".")], method = "s", exact = FALSE))
}
# median
cat("\n","Median final slopes:","\n")
for (i in c("markers","haplotypes","species","individuals")){
  cat(i, ":", median(slopes[, paste("m", i, sep = ".")]), "\n")
}
# range
cat("\n","Range of final slopes:","\n")
for (i in c("markers","haplotypes","species","individuals")){
  cat(i, ":", range(slopes[, paste("m", i, sep = ".")]), "\n")
}
# N = 0
cat("\n","N final slopes == 0:","\n")
for (i in c("markers","haplotypes","species","individuals")){
  cat(i, ":",length(which(slopes[, paste("m", i, sep = ".")] == 0)), "\n")
}
# N < 0.05
cat("\n","N final slopes < 0.05:","\n")
for (i in c("markers","haplotypes","species","individuals")){
  cat(i, ":",length(which(slopes[, paste("m", i, sep = ".")] < 0.05)), "\n")
}
# stop sinking to text file
sink()
# plot slopes for markers, haplotypes, species and individuals
options(scipen=999)
for(i in c("species","individuals","haplotypes","markers")){
  layout.mat <- matrix(c(rep(26,7),28,1:5,27,28,6:10,27,28,11:15,27,28,16:20,27,28,21:25,27,rep(29,7)),ncol=7)
  for(z in c("05","10","20","50")){
    pdf(paste("results/rarefaction/",i,"_rarefaction_ordered_",z,".pdf",sep=""))
    layout(layout.mat)
    par(mar=rep(0,4))
    for(x in c("A","B","C","D","E")){
      for(y in rev(c("1","2","3","4","5"))){
        # initial plot
        plot(c(0,rarefaction.dat$n.reads[,paste("R_",y,x,"_",z,sep="")]),
             c(0,rarefaction.dat[[i]]$mean[,paste("R_",y,x,"_",z,sep="")]),
             main="",
             xlab="",
             ylab="",
             xaxt='n',
             yaxt='n',
             pch=19,
             cex=0.4,
             xlim=c(0-10*min(rarefaction.dat$n.reads[,paste("R_",y,x,"_",z,sep="")]),max(rarefaction.dat$n.reads[,paste("R_",y,x,"_",z,sep="")])),
             ylim=c(0,max(rarefaction.dat[[i]]$mean[,paste("R_",y,x,"_",z,sep="")]))) 
        # Y axis
        my.at <- 0:ceiling(max(rarefaction.dat[[i]]$mean[,paste("R_",y,x,"_",z,sep="")])) 
        if(max(my.at)>10){
          my.at <- my.at[which(my.at %% 10 == 0)]
        } 
        if(max(my.at)>100){
          my.at <- my.at[which(my.at %% 100 == 0)]
        }
        if(max(my.at)>1000){
          my.at <- my.at[which(my.at %% 1000 == 0)]
        }
        if(max(my.at)>10000){
          my.at <- my.at[which(my.at %% 10000 == 0)]
        }
        axis(side = 2,
             at = my.at,
             tck=0.05,
             las=0,
             cex=0.5,
             padj = 4)
        #X axis if on bottom row
        if(y=="1"){
          axis(side=1,
               at=rarefaction.dat$n.reads[which(1:nrow(rarefaction.dat$n.reads) %% 10 == 0),paste("R_",y,x,"_",z,sep="")],
               labels = 1/nrow(rarefaction.dat$n.reads)*which(1:nrow(rarefaction.dat$n.reads) %% 10 == 0))
        }
        # add subplot letters on top edge
        if(y == "5"){
          mtext(text=x,
                side=3,
                line=1)
        }
        # add subplot numbers on right edge
        if(x == "E"){
          mtext(text=y,
                side=4,
                las=2,
                line=1)
        }
        # plot SD
        polygon(c(rarefaction.dat$n.reads[,paste("R_",y,x,"_",z,sep="")],
                  rev(rarefaction.dat$n.reads[,paste("R_",y,x,"_",z,sep="")])),
                c(rarefaction.dat[[i]]$upper95CI[,paste("R_",y,x,"_",z,sep="")],
                  rev(rarefaction.dat[[i]]$lower95CI[,paste("R_",y,x,"_",z,sep="")])),
                col="gray75",
                border=NA)
        # replot points
        points(rarefaction.dat$n.reads[,paste("R_",y,x,"_",z,sep="")],
               rarefaction.dat[[i]]$mean[,paste("R_",y,x,"_",z,sep="")],pch=19,cex=0.4)
        # get N reads and slope
        N <- slopes[which(slopes$sample == paste("R_",y,x,"_",z,sep="")),"n.cat.reads"]
        N <- format(N, big.mark=",",scientific=F)
        m <- slopes[which(slopes$sample == paste("R_",y,x,"_",z,sep="")),paste("m",i,sep=".")]
        if(m == 1){
          m <- 1
        } else {
          m <- format(round(m,2),nsmall = 2,scientific = F)
        }
        # write 10% final slope and number of individuals in full dataset
        mtext(bquote(italic("m") == .(m)),side = 1,line=-1,cex=0.7,col="firebrick",adj=1)
        mtext(bquote(italic("N") == .(N)),side = 1,line=-2,cex=0.7,col="cornflowerblue",adj=1)
      }
    }
    # add Y axis label
    plot(c(0,1),c(0,1),axes=F,col="white",xlab="",ylab="")
    text(0.75,0.5,paste("Number of ",i," detected"),cex=2,srt=90)
    # add X axis label
    plot(c(0,1),c(0,1),axes=F,col="white",xlab="",ylab="")
    text(0.5,0.5,"Proportion of total reads subsampled",cex=2)
    dev.off()
  }
}

