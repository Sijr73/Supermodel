
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
	
	sba <- sysBiolAlg(model, algorithm = "arm", additionalReact=addid, biomassThreshold= bioThreshold, solverParm=sp)
	res <- list()
	
	if(file.exists(splitfile)){
		cat("loading saved split results\n")
		load(splitfile)
	}
	count <- 0
	for(i in seq(along=media)){
		if(!is.null(res[[media[i]]])){
			next
		}
		opt <- optimizeProb(sba,
					react=1:react_num(model),
					ub=ub[[i]],
					lb=lb[[i]],
					poCmd=poCmd
		)
		
		res[[media[i]]] <- list(	solution=intersect(which(abs(getArmReactionFluxes(model, opt)) >= 1e-6), addid),
							relGap=pa(opt$postP)[[c(1, 1)]],
							absGap=pa(opt$postP)[[c(2, 1)]]-pa(opt$postP)[[c(3, 1)]],
							solStat=opt$stat)
		
		if(is.null(sba)){
			sba <- opt$sba
		}
		
		count <- count + 1
		if(count %% 20 == 0){
			cat("saving split results\n")
			save(res, file=splitfile)
		}
	}
	
	return(res)
}


















