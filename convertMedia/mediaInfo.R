#!/usr/bin/Rscript

library(methods)
library(parallel)
library(dplyr)
library(cplexAPI)

print(load("mediaBiGGwithRandom.Rdata"))


mediaInfo <- data.frame(mediumName=names(mediaBiGG))

mediaInfo$mediumType <- "normal"
mediaInfo$mediumType[grepl("^[CNPS]_", mediaInfo$mediumName)] <- "minimal"
mediaInfo$mediumType[grepl("^R\\|", mediaInfo$mediumName)] <- "random"

mediaInfo$mediumSize <- sapply(mediaBiGG, length)

mediaInfo$hasGlc <- sapply(mediaBiGG, function(x) "glc__D" %in% x)



















