#!/usr/bin/Rscript
library(methods)
library(sybil)
library(cplexAPI)
library(parallel)
library(dplyr)
if(!interactive()){
		args <- commandArgs(trailingOnly = TRUE)
		modelNumber <- as.integer(args[1])
		cpuNum <- 8L
		maxMem <- 30
}else{
	cpuNum <- 4
	maxMem <- 40
	if(!exists("modelNumber")){
		modelNumber <- 1
	}
}
#Rprof(filename = paste0("Rprof.", modelNumber, ".out"), append = FALSE, interval = 0.02,
#		memory.profiling = T, gc.profiling = T, 
#		line.profiling = T)

source("~/envirDist/runCrossFBA.R")
source("~/envirDist/initData.R")
print(load("~/DATA/mergedResults.Rdata"))
sessionInfo()
print(R.version)

uptakeLowBnd <- -10
optRounds <- 3
bioThreshold <- 1e-2

print(load("~/DATA/universalBiGG.ver2.Rdata"))
uni <- uni2

print(load("~/DATA/mediaBiGGwithRandom.Rdata"))

SYBIL_SETTINGS("SOLVER", "cplexAPI")

options(mc.cores=as.integer(cpuNum))
chunk <- function(x, n) split(x, factor(sort(rank(x)%%n)))

sp <- list()
sp$CPX_PARAM_WORKMEM <- floor((maxMem/cpuNum) * 1024)
sp$CPX_PARAM_TRELIM <- floor((maxMem/cpuNum) * 1024)
sp$CPX_PARAM_THREADS <- 1L
sp$CPX_PARAM_PARALLELMODE <- CPX_PARALLEL_OPPORTUNISTIC
sp$CPX_PARAM_EPINT <- 1e-9
sp$CPX_PARAM_EPRHS <- 1e-9
sp$CPX_PARAM_NUMERICALEMPHASIS <- CPX_ON
#sp$CPX_PARAM_MIPEMPHASIS <- CPX_MIPEMPHASIS_FEASIBILITY


stopifnot(initData())
mn <- names(modelReactMap)[modelNumber]
filename <- paste0("~/DATA/crossdata/", mn, ".Rdata")

modelReactions <- union(modelReactMap[[mn]], react_id(ex)) # allow all exchange reactions.
addedReactionsModel <- addedReactions[[mn]]

lowbnd(uni)[react_id(uni) %in% modelRmEnergyCycles[[mn]]$bwd] <- 0
uppbnd(uni)[react_id(uni) %in% modelRmEnergyCycles[[mn]]$fwd] <- 0


scratchdirpre <- "~/DATA/scratch_gs/tmp/cross/"
if(interactive()){
	scratchdirpre <- "~/DATA/tmpdir/"
}

scratchdirpre <- paste0(scratchdirpre, mn, "/")

if(!file.exists(scratchdirpre)){
	stopifnot(dir.create(scratchdirpre, recursive=T))
}

result <- result %>% filter(model== mn, NG)
tasks <- split(result, result$bof)
#tasks <- lapply(tasks, head) #TODO remove before flight

#if(!interactive()){
	cat(">>> removing unused objects\n")
	rm(uni2, resultWithModel, modelReactMap, modelReactMapById, 
		modelRmEnergyCycles, modelRmEnergyCyclesById, mediaLowBounds, 
		addedReactions, result, addSummary)
#}
gc()

if(!exists("cl")){
#cl <- makeCluster(cpuNum)
	cl<-makeCluster(cpuNum,type="FORK")###
}
clInit <- function(){
	library(sybil)
	library(cplexAPI)
	SYBIL_SETTINGS("SOLVER", "cplexAPI")
	return(T)
}
stopifnot(all(unlist(clusterCall(cl, clInit))))
setDefaultCluster(cl)


if(!file.exists(filename) || interactive()){
	for(e in names(tasks)){
		print(c(mn, e))
		bof <- switch(e,
			normal=modelBiomassMapSelection[mn], 
			general="GENERAL_BOF",
			energy="ENERGY_BOF"
		)
		if(is.na(bof) || length(tasks[[e]]$medium) == 0){
			df[[e]] <- data.frame(addedForMedium=NULL, medium=NULL, mediumType=NULL, model=NULL, bof=NULL, growth=NULL, status=NULL)
			next
		}
		
#		df[[e]] <- list()
		
		mediaToTestByType <- split(as.character(tasks[[e]]$medium), tasks[[e]]$mediumType)
		
		for(mediumType in names(mediaToTestByType)){
			mediaToTest <- mediaToTestByType[[mediumType]]
			for(addedForMedium in mediaToTest){
				print(addedForMedium)
				mediumfile <- paste0(scratchdirpre, e, ".", gsub("\\W", "_", addedForMedium), ".Rdata")
				print(mediumfile)
				if(file.exists(mediumfile) & !interactive()){
					#load(mediumfile)
					#cat("loaded file from scratch_gs\n")
					cat("results found on hdd:\n")
					print(mediumfile)
				}else{
					add <- addedReactionsModel[[c(e, addedForMedium)]]
					stopifnot(!is.null(add))
		
					model <- rmReact(changeObjFunc(uni, bof), setdiff(react_id(uni), union(modelReactions, add)))
					
					out <- runFBA(model=model, mn=mn, addedForMedium=addedForMedium, mediaToTest=mediaToTest, bof=e, cores=cpuNum)
#					browser(expr=(e=="general" & addedForMedium == "R|C|thymd|N|23cgmp|P|23camp|S|h2s"))
					out$mediumType <- mediumType
					save(out, file=mediumfile)
				}
#				df[[c(e, addedForMedium)]] <- out
			
			}
		}
#		df[[e]] <- do.call(rbind, df[[e]])
		gc()
	}
	
	cat(">>> collecting data now\n")
	df <- list()
	for(e in names(tasks)){
		df[[e]] <- list()
		mediaToTestByType <- split(as.character(tasks[[e]]$medium), tasks[[e]]$mediumType)
		for(mediumType in names(mediaToTestByType)){
			mediaToTest <- mediaToTestByType[[mediumType]]
			for(addedForMedium in mediaToTest){
				print(addedForMedium)
				mediumfile <- paste0(scratchdirpre, e, ".", gsub("\\W", "_", addedForMedium), ".Rdata")
				load(mediumfile)
				df[[c(e, addedForMedium)]] <- out
			}
		}
		df[[e]] <- do.call(rbind, df[[e]])
	}
	df <- do.call(rbind, df)
	rownames(df) <- NULL
	if(!interactive()){
		save(df, file=filename)
	}else{
		cat(">>> result not saved.\n")
	}
}else{
	cat(">>> file already existing.\n")
	print(filename)
}

stopCluster(cl)

cat(">>> EOF\n")
#ggplot(df, aes(growth, fill=model)) + geom_histogram(binwidth=0.2) + facet_grid( bof ~ run)



































