#!/usr/bin/Rscript
library(methods)
library(sybil)
library(cplexAPI)
#library(sybilSWITCH)
library(parallel)
source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")
source("New analysis bacillus/DATA/E.coli/sysBiolAlg_armClass.R")
source("~/envir/runARM.R")
source("New analysis bacillus/DATA/E.coli/runARM.R")
source("~/EGCecoli/initData1.R")
source("New analysis bacillus/DATA/E.coli/initData1.R")

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
#modelNumber <- SJ

options(mc.cores=as.integer(cpuNum))

stopifnot(initData())
bioThreshold <- 1e-2


mn <- names(modelReactMap)[modelNumber]
#if(!file.exists("envirdata")){
#	dir.create("New analysis bacillus/envirdata")
#}

#filename <- paste0("~/EGCecoli/envirdata/", mn, ".Rdata")
filename <- paste0("New analysis bacillus/DATA/E.coli/eni/", mn, ".Rdata")

#stop()

rm(uni2, modelReactMap, modelRmEnergyCycles)


if(!file.exists(filename) || interactive()){
	
	#print(load("~/EGCecoli/compareDf.Rdata"))
  print(load("New analysis bacillus/DATA/E.coli/compareDf.Rdata"))
  
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
		splitfile <- paste0("New analysis bacillus/DATA/E.coli/eni/split/", paste(mn, e, "Rdata", sep="."))
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
#data=df
print(load("New analysis bacillus/DATA/E.coli/universalBiGG.ver2.Rdata"))
files <- dir("New analysis bacillus/DATA/E.coli/eni/", pattern="*.Rdata", full.names=T)
data <- mclapply(files, function(f){
  get(load(f))
})
names(data) <- gsub("New analysis bacillus/DATA/E.coli/eni//", "", gsub(".Rdata", "", files, fixed=T), fixed=T)
addedReactions <- lapply(names(data), function(mn){
  e <- data[[mn]]
  erg <- lapply(names(e), function(etype){
    lapply(e[[etype]], function(x){
      react_id(uni2)[x$solution]
    })
  })
  names(erg) <- names(e)
  return(erg)
})


save(df, file=filename)
save(erg, file=splitfile)

cat(">>> EOF\n")
#ggplot(df, aes(growth, fill=model)) + geom_histogram(binwidth=0.2) + facet_grid( bof ~ run)



































