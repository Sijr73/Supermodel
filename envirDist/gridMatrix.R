#!/usr/bin/Rscript

library(methods)
library(parallel)
library(dplyr)
library(ggplot2)
library(Matrix)

node <- Sys.info()["nodename"]
if(grepl("^jedi|sith$", node)){
	cat(">>> config for jedi|sith\n")
	options(mc.cores=28L)
}

print(load("~/DATA/mediaBiGGwithRandom.Rdata"))
mediaBiGG <- mediaBiGG[grep("^R\\|", names(mediaBiGG))]
essentials <- names(which(table(unlist(mediaBiGG)) == length(mediaBiGG)))

mediaBiGG <- mclapply(mediaBiGG, function(x) setdiff(x, essentials))

cat("building gridMatrix\n")
mList <- mclapply(mediaBiGG, function(x){
ret <- sapply(mediaBiGG, function(y){
		length(intersect(x, y)) == 0
	})
})
gridMatrix <- do.call(rbind, mList)
rownames(gridMatrix) <- colnames(gridMatrix) <- names(mediaBiGG)



save(gridMatrix, file="~/DATA/gridMatrix.Rdata")





