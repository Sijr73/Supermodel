

runFBA <- function(model=NULL, mn=NULL, bof=NULL, cores=1L){
	if(is.null(model)){
		if(exists("uni")){
			model <- uni
		}
		else{
			stop("no model given")
		}
	}
	
	if(is.null(bof)){
		stop("bof mustnt be NULL")
	}
	
	#if no bof return empty data.frame
	if(is.na(bof) || length(bof) == 0){
		stop("no bof given!")
	}
	
	# setting bof and making list for optimizer
	obj <- rep(0, react_num(model))
	obj[bof] <- 1
	obj <- rep(list(obj), length(mediaLowBounds))

	# making list objects for optimizer function
	rid <- rep(list(1:react_num(model)), length(mediaLowBounds))
	ub <- rep(list(uppbnd(model)), length(mediaLowBounds))
	lb <- mediaLowBounds
	
	# setting cycles to off
	ub <- lapply(ub, function(x) {
		x[modelRmEnergyCyclesById[[c(mn, "fwd")]]] <- 0
		x
	})
	lb <- lapply(lb, function(x) {
		x[modelRmEnergyCyclesById[[c(mn, "bwd")]]] <- 0
		x
	})
	
	tasks <- chunk(seq(along=mediaBiGG), n=cores)
	ergFull <- mclapply(tasks, mc.cores=cores, function(t){
		return(list(mediaId=t, opt=optimizer(model, react=rid[t], lb=lb[t], ub=ub[t], obj_coef=obj[t], verboseMode = 1, solverParm = sp)))
	})
	
#	ergFull <- optimizer(model, react=rid, lb=lb, ub=ub, obj_coef=obj)
#	stopifnot(length(checkSolStat(ergFull$stat))==0)
#	checkFull <- checkSolStat(ergFull$stat)
	
	offid <- setdiff(1:react_num(model), modelReactMapById[[mn]])

	ub <- mclapply(ub, function(x){
		x[offid] <- 0
		x
	})
	lb <- mclapply(lb, function(x){
		x[offid] <- 0
		x
	})
	
	erg <- mclapply(tasks, mc.cores=cores, function(t){
		return(list(mediaId=t, opt=optimizer(model, react=rid[t], lb=lb[t], ub=ub[t], obj_coef=obj[t], verboseMode = 1, solverParm = sp)))
	})
	
#	erg <- optimizer(model, react=rid, lb=lb, ub=ub, obj_coef=obj)
#	stopifnot(length(checkSolStat(erg$stat))==0)
#	check <- checkSolStat(erg$stat)
	
	
	makeDf <- function(e, run=""){
		media <- names(mediaBiGG)[e$mediaId]
		value <- round(e$opt$obj, digits=6)
		return(data.frame(medium=media, run=run, model=mn, growth=value, status=e$opt$stat))
	}
	
	dfF <- do.call(rbind, lapply(ergFull, makeDf, run="full"))
	dfM <- do.call(rbind, lapply(erg, makeDf, run="sub"))
	
	return(rbind(dfF, dfM))

}





