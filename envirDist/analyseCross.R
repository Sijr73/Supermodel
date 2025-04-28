#!/usr/bin/Rscript

library(methods)
library(parallel)
library(dplyr)
library(ggplot2)
library(cplexAPI)

print(load("~/DATA/universalBiGG.ver2.Rdata"))

print(load("~/DATA/gridMatrix.Rdata"))

bioThreshold <- 1e-2

files <- dir("~/DATA/crossdata/", pattern="*.Rdata", full.names=T)

node <- Sys.info()["nodename"]
if(grepl("^jedi|sith$", node)){
	cat(">>> config for jedi|sith\n")
	options(mc.cores=28L)
}else{
	options(mc.cores=4)
}


filesToDelete <- character(0)
jobsToRestart <- integer(0)
scratchdirpre <- "~/DATA/scratch_gs/tmp/cross/"

filename <- "~/DATA/innovationIndex1.Rdata"
writePDF <- T

if(! file.exists(filename) || ! interactive()){
	dflist <- mclapply(files, function(f){
		print(f)
		df <- get(load(f))
		dim(df)
	
		df <- df %>% mutate(growth = ifelse(status==CPX_STAT_INFEASIBLE, 0, growth))
	
		if(nrow(df %>% filter(status == CPX_STAT_OPTIMAL_INFEAS, growth >= 1e-8))> 0){
			probDf <- df %>% filter(status == CPX_STAT_OPTIMAL_INFEAS, growth >= 1e-8) %>% print()
		
			mn <- probDf$model[1]
			mediumfiles <- paste0(paste0(scratchdirpre, mn, "/"), probDf$bof, ".", gsub("\\W", "_", probDf$addedForMedium), ".Rdata")
			filesToDelete <<- append(filesToDelete, mediumfiles)
			jobsToRestart <<- unique(c(jobsToRestart, which(names(modelReactMap)==mn)))
			warning("problem with some solutions")
		}
	
		stopifnot(all(df$status %in% c(CPX_STAT_OPTIMAL_INFEAS, CPX_STAT_INFEASIBLE, CPX_STAT_OPTIMAL_INFEAS, CPX_STAT_OPTIMAL)))
		df <- df %>% mutate(isViable=growth>=bioThreshold)
	
		if((df %>% filter(as.character(medium) == as.character(addedForMedium), !isViable) %>% select(isViable) %>% nrow()) > 0){
			probDf <- df %>% filter(as.character(medium) == as.character(addedForMedium), !isViable) %>% select(isViable)
			mn <- probDf$model[1]
			mediumfiles <- paste0(paste0(scratchdirpre, mn, "/"), probDf$bof, ".", gsub("\\W", "_", probDf$addedForMedium), ".Rdata")
			filesToDelete <<- append(filesToDelete, mediumfiles)
			jobsToRestart <<- unique(c(jobsToRestart, which(names(modelReactMap)==mn)))
			warning("problem with some solutions")
		}
		
		df$distinctNutrients <- NA
		df[df$mediumType == "random", "distinctNutrients"] <- gridMatrix[as.matrix(df %>% filter(mediumType == "random") %>% select(addedForMedium, medium))]
		
#		browser()
		
		retDf <- df %>% group_by(model, bof, mediumType, addedForMedium) %>% 
							summarize(innoIndex = sum(isViable),
										innoIndexRel = sum(isViable)/n(),
										innoCount = n(),
										innoIndexDistinct = sum(isViable & distinctNutrients),
										innoIndexDistinctRel = innoIndexDistinct / sum(distinctNutrients),
										innoCountDistinct = sum(distinctNutrients))
		
		return(retDf)
	})

	innovationIndex <- do.call(rbind, dflist)
	innovationIndex <- innovationIndex %>% ungroup() %>% mutate_if(is.character, factor)
	gc()

	stopifnot(length(filesToDelete) == 0)

	stopifnot(setequal(names(modelReactMap), innovationIndex %>% select(model) %>% unique() %>% unlist()))

	#  We define the innovation index, IGlucose, of a network to be the number of 
	#  additional carbon sources on which each network is viable.
	
	print(load("New analysis bacillus/mergedResults.Rdata"))
	innovationIndex <- innovationIndex %>% left_join(addSummary %>% select(model, bof, mediumType, NN, NG, GN, GG, count), by=c("model", "bof", "mediumType"))
#	adaptCoef := innoIndex / (GG+1)
#	innovationIndex <- innovationIndex %>% mutate(adaptCoef=innoIndex/(GG+1))
	
	save(innovationIndex, file=filename)
}else{
	load(filename)
}

saveAsPDF <- function(filename=null, plot=g){
	stopifnot(!is.null(filename))
	if(!interactive() || writePDF){
			ggsave(plot, filename=filename, path="plots", width=297, height=210, units="mm")
	}else{
			#print(plot)
	}
}

# violins of relative innovation index
g <- ggplot(innovationIndex, aes(x=factor(model), y=innoIndexRel, colour=model)) +
	geom_violin() +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
#	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "none",
		axis.text.x = element_text(angle = 90, hjust = 1, size=7)) +
	ggtitle("Innovation Index")

saveAsPDF("innovationIndexByModel.pdf", g)

# density of relative innovation index
#g <- ggplot(innovationIndex, aes(fill=factor(model), x=innoIndexRel, colour=model)) +
#	geom_density(alpha = 0.1, trim=T) +
#	guides(colour=guide_legend(title.position = "top", ncol=9), 
#			shape=guide_legend(title.position = "top", nrow=1)) +
##	scale_y_log10() +
#	facet_grid(mediumType ~ bof) +
#	theme(legend.position = "none",
#		axis.text.x = element_text(angle = 90, hjust = 1, size=7)) +
#	ggtitle("Innovation Index")

#saveAsPDF("innovationIndexByModel.pdf", g)





































