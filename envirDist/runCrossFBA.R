
#objectSizes <- function(env=as.environment(-1L)){
#	sizes <- sapply(ls(envir=env), function(x) object.size(get(x)))
#	if(identical(globalenv(), env)){
#		return(sizes)
#	}else{
#		return(append(objectSizes(parent.env(env)), sizes))
#	}
#}


runFBA <- function(model=NULL, mn=NULL, addedForMedium=NULL, mediaToTest=NULL, bof=NULL, cores=1L){
	if(is.null(model)){
		stop("no model given")
	}
	
	if(is.null(bof)){
		stop("bof mustnt be NULL")
	}
	
	stopifnot(!is.null(addedForMedium))
	stopifnot(!is.null(mediaToTest))
	
	mediaBiGG <- mediaBiGG[mediaToTest]
	
	ex <- findExchReact(model)
	lowbnd(model)[react_pos(ex)] <- 0
	exOnly <- ex[grep("^EX_", react_id(ex))]
	exMap <- react_pos(exOnly)
	names(exMap) <- gsub("\\[\\w\\]$", "", met_id(exOnly))
	
	# setting media constraints
	mediaLowBounds <- lapply(mediaBiGG, function(m){
		lb <- rep(0, length(exMap))
		lb[names(exMap) %in% m] <- uptakeLowBnd
		return(lb)
	})
	
	# making list objects for optimizer function
	rid <- rep(list(exMap), length(mediaLowBounds))
	ub <- rep(list(uppbnd(model)[exMap]), length(mediaLowBounds))
	lb <- mediaLowBounds
	
	tasks <- chunk(seq(along=mediaBiGG), n=cores)
#	objectSizesValues <<- objectSizes()
#	dput(objectSizesValues)
	
	stopifnot(all(clusterExport(cl=NULL, varlist=c("model", "rid", "lb", "ub", "sp"), envir=environment())))
	erg <- clusterApply(cl=NULL, tasks, function(t){
#	erg <- mclapply(tasks, mc.cores=cores, mc.cleanup=T, mc.preschedule=T, function(t){
#	erg <- lapply(tasks, function(t){
		return(list(mediaId=t, opt=optimizer(model, react=rid[t], lb=lb[t], ub=ub[t], verboseMode = 1, solverParm = sp)))
	})
#	browser()
	
	makeDf <- function(e){
		media <- names(mediaBiGG)[e$mediaId]
		value <- round(e$opt$obj, digits=6)
		return(data.frame(addedForMedium=addedForMedium, medium=media, model=mn, bof=bof, growth=value, status=e$opt$stat))
	}
	
	df <- do.call(rbind, lapply(erg, makeDf))
	rownames(df) <- NULL
	return(df)

}









