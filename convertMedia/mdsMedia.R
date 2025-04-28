#!/usr/bin/Rscript
library(parallel)



print(load("mediaBiGGwithRandom.Rdata"))

jacDist <- function(x, y){
	if(length(x) == 0 && length(y)==0){
		return(0)
	}
	1 - (length(intersect(x, y)) / length(union(x, y)))
}

mediaTypes <- list()
mediaTypes$minimal <- grep("^[CNPS]_cpd", names(mediaBiGG), value=T)
mediaTypes$random <- grep("^R\\|", names(mediaBiGG), value=T)
mediaTypes$seed <- setdiff(names(mediaBiGG), unlist(mediaTypes))

mds <- list()


for(type in names(mediaTypes)){
	m <- mclapply(mediaTypes[[type]], function(i){
		ji <- mediaBiGG[[i]]
		sapply(mediaTypes[[type]], function(j){
			return(jacDist(ji, mediaBiGG[[j]]))
		})
	})
	distMatrix <- do.call(rbind, m)
	
	
	mds[[type]] <- cmdscale(distMatrix, eig=TRUE, k=2)
	
	
}

mdsPoints <- lapply(names(mds), function(n){
	m <- mds[[n]]
	df <- m$points
	colnames(df) <- c("x", "y")
	rownames(df) <- mediaTypes[[n]]
	return(as.data.frame(df))
})
names(mdsPoints) <- names(mds)


save(mds, mdsPoints, file="mdsMedia.Rdata")


































