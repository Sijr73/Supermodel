
# runARM(model=uni, mn=mn, bof=bof, media= m, cores=1)





runARM <- function(model=NULL, mn=NULL, bof=NULL, media=NULL, cores=1L){
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
	
	if(is.null(media)){
		return(list())
	}
	
	
	model <- changeObjFunc(model, bof)
	
	# making list objects for optimizer function
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
	
	ub <- ub[match(media, names(mediaBiGG))]
	lb <- lb[match(media, names(mediaBiGG))]
	
	addid <- setdiff(1:react_num(model), modelReactMapById[[mn]])
	
	tasks <- chunk(seq(along=media), n=cores)
#	tasks <- lapply(tasks, head)
	
	erg <- mclapply(tasks, mc.cores=cores, function(t){
#	erg <- lapply(tasks, function(t){
		sba <- NULL
		res <- list()
		
		for(i in t){
			opt <- optimizeProbArmLP(model=model, sysBiolAlg=sba,
						additionalReact=addid,
						react=1:react_num(model),
						ub=ub[[i]],
						lb=lb[[i]],
						biomassThreshold= bioThreshold,
						randomAfterOptimization= randomAfter,
						maxOptimizationRounds= optRounds*randomAfter,
						returnSBA = ifelse(is.null(sba), T, F),
						threads=1,
						printDebug=interactive(),
						solverParm=sp
			)
			
			browser()
			
			res[[media[i]]] <- list(	solution=intersect(opt$solution, addid),
								solStat=opt$solStat)
			
			if(is.null(sba)){
				sba <- opt$sba
			}
		}
		return(res)
	})
	names(erg) <- NULL
	erg <- unlist(erg, recursive=F)
	return(erg)
}
