initData <- function(){
	uptakeLowBnd <- -10
	bioThreshold <- 1e-2
	#bioThreshold <- 10
	
	print(load("~/DATA/universalBiGG.ver2.Rdata"))
	uni <- uni2

	print(load("~/DATA/mediaBiGGwithRandom.Rdata"))
modelReactMap[["iIS312"]]=append(modelReactMap[["iIS312"]],c("UMPK_B", "ADKd", "DADNt2_LR", "DGSNt2_LR","NTD8_B", 
                                                             "TYRabc", "TRPabc", "AMPt6", "CMPt6", "DTMPt6", "DURIK1_1"
	                                                             ,"DURIt2_B", "MANpts_B" ))
modelReactMap[["iIS312_Amastigote"]]=append(modelReactMap[["iIS312_Amastigote"]],c("DADK_B", "DGSNt","DTTPt","DCTPD2","DADNt2_LR",
                                                                        "TYRabc","TRPabc", "MAN6Pt6_2", "AMPt6","ADADir_LR",
                                                             "CMPt6","DCYTt2_B","DGNSK_1"  ))	
	modelReactMap[["iIS312_Epimastigote"]]=append(modelReactMap[["iIS312_Epimastigote"]],c("UMPK_B", "ADKd","DTTPt","DGSNt2_LR",
	"NTD8_B","TYRabc","TRPabc","ADADir_B","AMPt6","CMPt6", "MAN1Pt6","DADNt2_B","DCYTt2_B"))	
	modelReactMap[["iIS312_Trypomastigote"]]=append(modelReactMap[["iIS312_Trypomastigote"]],c("ADKd","UMPK_copy2_LR",
	"DGSNt2_LR","TYRabc","TRPabc","F6Pt6_2","AMPt6","CMPt6","DTMPt6","DADNt2_B","DGNSK_1","DURIK1_1","DURIt2_B"))
	modelReactMap[["iLJ478"]]=append(modelReactMap[["iLJ478"]],"ASNt2r_LR")
modelReactMap[["iSM199"]]=append(modelReactMap[["iSM199"]],c("LEUTAi","ILEabc","VALt2r_LR","ILETA2_LR","PHEabc"))
	
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
		union(match(x, react_id(uni)), exMap)#  adding all exchange reactions to the models.
	})

	modelBiomassMapSelectionById <- match(modelBiomassMapSelection, react_id(uni))
	names(modelBiomassMapSelectionById) <- names(modelBiomassMapSelection)
	
	modelRmEnergyCyclesById <- lapply(modelRmEnergyCycles, function(x){
		lapply(x, function(y) match(y, react_id(uni)))
	})
	
	for(i in ls()) assign(x=i, value=get(i), envir=.GlobalEnv)
	return(T)
}

