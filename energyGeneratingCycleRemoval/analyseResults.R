#!/usr/bin/Rscript

library(methods)
library(sybil)
library(cplexAPI)
#library(sybilSWITCH)
library(ggplot2)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")
threadCount <- 7


print(load("~/DATA/universalBiGG.ver1.2.Rdata"))
#print(load("~/DATA/E.coli/universalBiGG.ver1.2.Rdata"))
uni <- uni1.2
ex <- findExchReact(uni)
lowbnd(uni)[react_pos(ex)] <- 0
uni <- changeObjFunc(uni, grep("^FCC_", react_id(uni)))

modelReactMapId <- lapply(modelReactMap, match, table=react_id(uni))


cycleResults <- lapply(names(modelReactMap), function(n){
	filename <- paste0("~/energyGeneratingCycleRemoval/resultsfinal/egcResults", n, ".Rdata")
	if(file.exists(filename)){
		(load(filename))
		return(list(cycles=cycles, removedBwd=removedBwd, removedFwd=removedFwd))
	}else{
		return(NULL)
	}
})
names(cycleResults) <- names(modelReactMap)

stopifnot(all(!sapply(cycleResults, is.null)))


egcTest <- lapply(names(cycleResults), function(n){
	fwd <- unlist(cycleResults[[c(n, "removedFwd")]])
	bwd <- unlist(cycleResults[[c(n, "removedBwd")]])
	
	ub <- uppbnd(uni)
	lb <- lowbnd(uni)
	
	ub[fwd] <- 0
	lb[bwd] <- 0
	
	optimizeProb(uni, react=1:react_num(uni), ub=ub, lb=lb, retOptSol=F)
})
names(egcTest) <- names(cycleResults)
stopifnot(length(checkSolStat(sapply(egcTest, "[[", "stat")))==0)


stopifnot(all(round(sapply(egcTest, "[[", "obj"), digits=6) == 0))


if(T){
	biomassTest <- lapply(names(cycleResults), function(n){
		fwd <- unlist(cycleResults[[c(n, "removedFwd")]])
		bwd <- unlist(cycleResults[[c(n, "removedBwd")]])
	
		ub <- uppbnd(uni1.2)
		lb <- lowbnd(uni1.2)
	
		ub[fwd] <- 0
		lb[bwd] <- 0
	
		bof <- modelBiomassMapSelection[n]
		if(is.na(bof)){
			m <- uni1.2
		}else{
			m <- changeObjFunc(uni1.2, bof)
		}
	
		optimizeProb(m, react=1:react_num(uni), ub=ub, lb=lb, retOptSol=F)
	})
	names(biomassTest) <- names(cycleResults)
	stopifnot(length(checkSolStat(sapply(biomassTest, "[[", "stat")))==0)


	if(!interactive()){
		pdf("cycleRemovalReport.biomassProduction.pdf", width=7, paper="a4r")
		barplot(sapply(biomassTest, "[[", "obj"))
		dev.off()
	}else{
		barplot(sapply(biomassTest, "[[", "obj"))
	}
}



cycleCounts <- do.call(rbind, lapply(names(modelReactMap), function(i){
	if(is.null(cycleResults[[i]])){
		return(NULL)
	}
	counts <- sapply(cycleResults[[c(i, "cycles", i)]], length)
	data.frame(baseModel=i, addedModels=names(counts), count=counts)
}))




(p <- ggplot(cycleCounts, aes(baseModel, addedModels)) +
		geom_tile(aes(fill = count), colour = "white") +
		scale_fill_gradient(low = "white", high = "steelblue") +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)))
		
		
if(!interactive()){
	pdf("cycleRemovalReport.pdf", width=7, paper="a4r")
	print(p)
	dev.off()
}else{
	print(p)
}


modelRmEnergyCycles <- lapply(names(modelReactMap), function(n){
	fwd <- unlist(cycleResults[[c(n, "removedFwd")]])
	bwd <- unlist(cycleResults[[c(n, "removedBwd")]])
	
	
	list(fwd=react_id(uni)[fwd], bwd=react_id(uni)[bwd])
})
names(modelRmEnergyCycles) <- names(modelReactMap)


modelReactMap <- lapply(modelReactMap, function(x) union(x, "ENERGY_BOF"))

uni2 <- uni1.2
save(uni2, modelReactMap, modelBiomassMap, modelBiomassMapSelection, modelRmEnergyCycles, file="~/DATA/universalBiGG.ver2.Rdata")
#save(uni2, modelReactMap, modelBiomassMap, modelBiomassMapSelection, modelRmEnergyCycles, file="~/DATA/E.coli/universalBiGG.ver2.Rdata")












