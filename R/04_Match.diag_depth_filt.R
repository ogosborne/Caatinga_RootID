library(RootID)
# load best diagnostic markers
load("results/processed_data/best.diag.RData")
# run match.diag with a depth filter of 3
best.matches.3.3 <- match.diag(data.dir = "STACKS_data//US_M6/CS_n7/root", diag = best.diag, min.reads.mar = 3, min.reads.hap = 3)
# save result
save(best.matches.3.3, file = "results/processed_data/best.match.md3.RData")

