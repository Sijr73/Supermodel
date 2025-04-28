applyMediaBoundsOnModel <- function(model, medium, uptakeLowBnd=-10){
	
	ex <- findExchReact(model)
	
	lowbnd(model)[react_pos(ex)] <- 0
	
	exOnly <- ex[grep("^EX_", react_id(ex))]
	exMap <- react_pos(exOnly)
	names(exMap) <- gsub("\\[\\w\\]$", "", met_id(exOnly))
	
	lb <- rep(0, length(exMap))
	lb[names(exMap) %in% medium] <- uptakeLowBnd
	
	lowbnd(model)[exMap] <- lb
	
	return(model)
}
