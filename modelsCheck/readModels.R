#!/usr/bin/Rscript

library(methods)
##use sybilSBML version 2 from its archive
library(sybilSBML)
library(parallel)

f <- dir("~/sourceData/models", pattern="*.xml", full.names=T)
##for E.coli super model
#f <- dir("~/Ecolisupermodel/xml_models", pattern="*.xml", full.names=T)

m <- mclapply(f, function(x){
                tryCatch(readSBMLmod(x), error=function(e) NULL )
        }
)

names(m) <- gsub("\\.xml", "", gsub(".*/models//", "", f))

cat("Problem with following models:\n")
err <- !sapply(m, is, "modelorg")
print(f[err])

m <- m[sapply(m, is, "modelorg")]


save(m, file="~/DATA/models.Rdata")
