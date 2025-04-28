initData <- function(){
	uptakeLowBnd <- -10
	bioThreshold <- 1e-2
	#bioThreshold <- 10
	
	#print(load("New analysis bacillus/networkComplexityBigg-master/models/universalBiGG.ver2.Rdata"))
	#print(load("~/EGCecoli/universalBiGG.ver2.Rdata"))
	print(load("New analysis bacillus/DATA/E.coli/universalBiGG.ver2.Rdata"))

	uni <- uni2

#	print(load("~/envir/mediaBiGGwithRandom.Rdata"))
	print(load("New analysis bacillus/DATA/E.coli/mediaBiGGwithRandom.Rdata"))

	#print(load("~/envir/mediaBiGG.Rdata"))
	

	ex <- findExchReact(uni)
	exOnly <- ex[grep("^EX_", react_id(ex))]
	exMap <- react_pos(exOnly)
	names(exMap) <- gsub("\\[\\w\\]$", "", met_id(exOnly))
	
	# setting media constraints
	mediaLowBounds <- lapply(mediaBiGG, function(m){
		lb <- lowbnd(uni)
		lb[exMap] <- 0
		lb[exMap[names(exMap) %in% m]] <- uptakeLowBnd
		return(lb)
	})

	modelReactMapById <- lapply(modelReactMap, function(x){
		union(match(x, react_id(uni)), exMap)# ... adding all exchange reactions to the models.
	})

	modelBiomassMapSelectionById <- match(modelBiomassMapSelection, react_id(uni))
	names(modelBiomassMapSelectionById) <- names(modelBiomassMapSelection)
	
	modelRmEnergyCyclesById <- lapply(modelRmEnergyCycles, function(x){
		lapply(x, function(y) match(y, react_id(uni)))
	})
	
	for(i in ls()) assign(x=i, value=get(i), envir=.GlobalEnv)
	return(T)
}

