rm(list=ls(all=TRUE))
for(i in 1:55){
  modelNumber <- i

#!/usr/bin/Rscript
library(methods)
library(dplyr)
library(sybil)
library(cplexAPI)
#library(sybilSWITCH)
library(parallel)
source("New analysis bacillus/networkComplexityBigg-master/envirDist/runARMMILP.R")
#source("New analysis bacillus/networkComplexityBigg-master/envirDist/initData.R")
source("New analysis bacillus/DATA/E.coli/initData1.R")
source("New analysis bacillus/DATA/E.coli/sysBiolAlg_armClass.R")
#source("New analysis bacillus/DATA/E.coli/algorithms.R")

SYBIL_SETTINGS("SOLVER", "cplexAPI")
SYBIL_SETTINGS("TOLERANCE", 1e-6)

cpuNum <- 1
maxMem <- 100
optRounds <- 3
randomAfter <- 4
chunk <- function(x, n) split(x, factor(sort(rank(x)%%n)))

sp <- list()
cplexThreads <- 8L
sp$CPX_PARAM_WORKMEM <- floor((maxMem/cpuNum/cplexThreads) * 1024)
sp$CPX_PARAM_TRELIM <- floor((maxMem/cpuNum/cplexThreads) * 1024)
sp$CPX_PARAM_THREADS <- cplexThreads
sp$CPX_PARAM_EPRHS <- 1e-9
sp$CPX_PARAM_EPINT <- 1e-9
sp$CPX_PARAM_EPOPT <- 1e-9
sp$CPX_PARAM_NUMERICALEMPHASIS <- CPX_ON
sp$CPX_PARAM_TILIM <- 4*60 # timelimit 2 min

poCmd = list(c("getMIPrelGapCPLEX", "LP_PROB@oobj@env", "LP_PROB@oobj@lp"),
                         c("getBestObjValCPLEX", "LP_PROB@oobj@env", "LP_PROB@oobj@lp"),
                         c("getObjValCPLEX", "LP_PROB@oobj@env", "LP_PROB@oobj@lp")
)

if(!interactive()){
		args <- commandArgs(trailingOnly = TRUE)
		modelNumber <- 73
}else{
	cpuNum <- 1
	if(!exists("modelNumber")){
		modelNumber <- 1
	}
}
modelNumber <- i

options(mc.cores=as.integer(cpuNum))

stopifnot(initData())
mn <- names(modelReactMap)[modelNumber]
#if(!file.exists("New analysis bacillus/DATA/envirdata1")){
#	dir.create("New analysis bacillus/DATA/envirdata1")
#}

#filename <- paste0("envirDistMILP.", mn, ".Rdata")
filename <- paste0("New analysis bacillus/DATA/E.coli/eni/Normal/", mn, ".Rdata")


rm(uni2, modelReactMap, modelRmEnergyCycles)

if(!file.exists(filename) || interactive()){
	
#	print(load("New analysis bacillus/DATA/compareDf.Rdata"))
  print(load("New analysis bacillus/DATA/E.coli/compareDf.Rdata"))
  
	# only NG cases and currently used model.
	#compareDf <- compareDf %>% filter(model==mn, fullG, bof=="normal")
	compareDf <- compareDf[compareDf$model==mn & compareDf$NG & compareDf$mediumType=="seed", ]
#	compareDf <- compareDf %>% filter(medium %in% c("N_cpd00043", "R|C|dcyt|N|tyrp|P|acgam1p|S|dms", "R|C|gal|N|dtmp|P|tyrp|S|h2s", 
#												"R|C|mbdg|N|asp__L|P|pser__L|S|dms", "R|C|metsox_R__L|N|spmd|P|g3pc|S|taur/taur_V3"
#										))
	
	mediaToUse <- lapply(split(compareDf$medium, compareDf$bof), as.character)
#	mediaToUse <- lapply(mediaToUse, function(x) sample(x, size=40)) #TODO remove before flight
	
	mediaBiGG <- mediaBiGG[unique(unlist(mediaToUse))]
	mediaLowBounds <- mediaLowBounds[unique(unlist(mediaToUse))]
	df <- list()
	for(e in c( "normal")){
		print(c(mn, e))
		splitfile <- paste0("New analysis bacillus/DATA/E.coli/eni/split/", paste(mn, e, "Rdata", sep="."))
		print(splitfile)
		
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
	}
#	if(!interactive()){
		save(df, file=filename)
	#}else
	{	cat(">>> result not saved.\n")
	}
}else{
	cat(">>> file already existing.\n")
	print(filename)
}
}
cat(">>> EOF\n")
#ggplot(df, aes(growth, fill=model)) + geom_histogram(binwidth=0.2) + facet_grid( bof ~ run)



































