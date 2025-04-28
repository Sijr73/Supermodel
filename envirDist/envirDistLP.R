#!/usr/bin/Rscript

library(methods)
library(sybil)
library(cplexAPI)
#library(sybilSWITCH)
library(parallel)
source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")
source("~/envirDist/runARM.R")
source("~/envirDist/initData.R")

SYBIL_SETTINGS("SOLVER", "cplexAPI")
SYBIL_SETTINGS("TOLERANCE", 1e-10)

cpuNum <- 4
maxMem <- 80
optRounds <- 3
randomAfter <- 4
chunk <- function(x, n) split(x, factor(sort(rank(x)%%n)))

sp <- list()
cplexThreads <- 1L
sp$CPX_PARAM_WORKMEM <- floor((maxMem/cpuNum/cplexThreads) * 1024)
sp$CPX_PARAM_TRELIM <- floor((maxMem/cpuNum/cplexThreads) * 1024)
sp$CPX_PARAM_THREADS <- cplexThreads
sp$CPX_PARAM_EPRHS <- 1e-9
sp$CPX_PARAM_EPINT <- 1e-9
sp$CPX_PARAM_EPOPT <- 1e-9
sp$CPX_PARAM_NUMERICALEMPHASIS <- CPX_ON

if(!interactive()){
		args <- commandArgs(trailingOnly = TRUE)
		modelNumber <- as.integer(args[1])
}else{
	cpuNum <- 4
	if(!exists("modelNumber")){
		modelNumber <- 1
	}
}

options(mc.cores=as.integer(cpuNum))

stopifnot(initData())
bioThreshold <- 1e-2


mn <- names(modelReactMap)[modelNumber]
if(!file.exists("envirdata")){
	dir.create("~/")
}

filename <- paste0("~/", mn, ".Rdata")

stop()

rm(uni2, modelReactMap, modelRmEnergyCycles)


if(!file.exists(filename) || interactive()){
	
	print(load("~/DATA/compareDf.Rdata"))
	# only NG cases and currently used model.
	compareDf <- compareDf[compareDf$model==mn & compareDf$NG, ]
	#compareDf <-compareDf[compareDf$model==mn & compareDf$GG, ]
	mediaToUse <- sapply(split(compareDf$medium, compareDf$bof), as.character)
#	mediaToUse <- lapply(mediaToUse, function(x) sample(x, size=40)) #TODO remove before flight
	
	mediaBiGG <- mediaBiGG[unique(unlist(mediaToUse))]
	mediaLowBounds <- mediaLowBounds[unique(unlist(mediaToUse))]
	
	df <- list()
	for(e in c("normal", "general", "energy")){
		print(c(mn, e))
		splitfile <- paste0("~/split/", paste(mn, e, "Rdata", sep="."))
		print(splitfile)
		
		if(file.exists(splitfile) && !interactive()){
			cat(">>> restoring calculations from file:\n")
			print(splitfile)
			load(splitfile)
			df[[e]] <- erg
		}
		else{
			bof <- switch(e,
				normal=modelBiomassMapSelectionById[mn], 
				general=which(react_id(uni) == "GENERAL_BOF"),
				energy=which(react_id(uni) == "ENERGY_BOF")
			)
			if(!is.na(bof)){
				erg <- runARM(model=uni, mn=mn, bof=bof, media=mediaToUse[[e]], cores=cpuNum)
			}else{
				erg <- list()
			}
			df[[e]] <- erg
			cat(">>> saving results in file:\n")
			print(splitfile)
			if(interactive()){
				cat(">>> splitfile not saved\n")
			}else{
				save(erg, file=splitfile)
			}
		}
	}
	if(!interactive()){
		save(df, file=filename)
	}else{
		cat(">>> result not saved.\n")
	}
}else{
	cat(">>> file already existing.\n")
	print(filename)
}
save(df, file=filename)
save(erg, file=splitfile)

cat(">>> EOF\n")
#ggplot(df, aes(growth, fill=model)) + geom_histogram(binwidth=0.2) + facet_grid( bof ~ run)



































