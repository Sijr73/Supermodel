#!/usr/bin/Rscript
library(methods)
library(sybil)
library(cplexAPI)
#library(sybilSWITCH)
library(parallel)
source("~/envir/runFBA.R")
source("~/EGC/initData1.R")
source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")

SYBIL_SETTINGS("SOLVER", "cplexAPI")
stopifnot(initData())
uni=uni2
cpuNum <- 8
cpuNum <- 4
maxMem <- 60
chunk <- function(x, n) split(x, factor(sort(rank(x)%%n)))

if(!interactive()){
		args <- commandArgs(trailingOnly = TRUE)
		modelNumber <- as.numeric(args[1])
}else{
	cpuNum <- 4
	if(!exists("modelNumber")){
		modelNumber <- 1
	}
}
options(mc.cores=as.integer(cpuNum))

sp <- list()
sp$CPX_PARAM_WORKMEM <- floor((maxMem/cpuNum) * 1024)
sp$CPX_PARAM_TRELIM <- floor((maxMem/cpuNum) * 1024)
sp$CPX_PARAM_THREADS <- 1L
sp$CPX_PARAM_PARALLELMODE <- CPX_PARALLEL_OPPORTUNISTIC
sp$CPX_PARAM_EPINT <- 1e-9
sp$CPX_PARAM_EPRHS <- 1e-9
#sp$CPX_PARAM_EPOPT <- 1e-9
#sp$CPX_PARAM_NUMERICALEMPHASIS <- CPX_ON



stopifnot(initData())
mn <- names(modelReactMap)[modelNumber]
filename <- paste0("~/EGC/FBAdata/", mn, ".Rdata")

if(!file.exists(filename) || interactive()){
	df <- list()
	for(e in c("normal", "general", "energy")){
		print(c(mn, e))
		bof <- switch(e,
			normal=modelBiomassMapSelectionById[mn], 
			general=which(react_id(uni) == "GENERAL_BOF"),
			energy=which(react_id(uni) == "ENERGY_BOF")
		)
		if(is.na(bof)){
			df[[e]] <- data.frame(medium=NULL, run=NULL, model=NULL, growth=NULL, status=NULL)
		}else{
			df[[e]] <- runFBA(model=uni, mn=mn, bof=bof, cores=cpuNum)
			df[[e]]$bof <- e
		}
	}
	df <- do.call(rbind, df)
	if(!interactive()){
		save(df, file=filename)
	}else{
		cat(">>> result not saved.\n")
	}
}else{
	cat(">>> file already existing.\n")
	print(filename)
}
cat(">>> EOF\n")



































