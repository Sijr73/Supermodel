#!/usr/bin/Rscript
library(methods)
library(sybil)
library(cplexAPI)
library(sybilSWITCH)
library(parallel)
source("initData.R")
source("../helper/applyMediaBoundsOnModel.R")

if(!exists("modelNumber")){
	modelNumber <- 1
	mn <- names(modelReactMap)[modelNumber]
}

erg <- lapply(names(df), function(e){
	bof <- switch(e,
				normal=modelBiomassMapSelectionById[mn], 
				general=which(react_id(uni) == "GENERAL_BOF"),
				energy=which(react_id(uni) == "ENERGY_BOF")
	)
	euni <- changeObjFunc(uni, bof)
	lapply(names(df[[e]]), function(m){
		print(c(e, m))
		optErg <- df[[c(e, m)]]
		
		mod <- rmReact(euni, setdiff(1:react_num(euni), union(modelReactMapById[[mn]], optErg$solution)))
		mod <- applyMediaBoundsOnModel(mod, mediaBiGG[[m]])
		
		opt <- optimizeProb(mod, solverParm=sp, retOptSol=F)
		
		data.frame(mn, e, m, obj=opt$obj, stat=opt$stat)
	})
})

erg <- do.call(rbind, unlist(erg, recursive=F))


if(any(erg$obj < bioThreshold)){
	erg %>% filter(obj < bioThreshold) %>% print()
	stop("some optimizations failed!")
}















