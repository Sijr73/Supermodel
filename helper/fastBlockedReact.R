fastBlockedReact <- function(model, tol=SYBIL_SETTINGS("TOLERANCE"), exex=TRUE){
	stopifnot(is(model, "modelorg"))
	
	obj_coef(model) <- rep(0, react_num(model))
	sba <- sysBiolAlg(model, algorithm="fba")
	
	toCheck <- rep(T, react_num(model))
	blocked <- rep(F, react_num(model))
	
	if(exex){
		toCheck[react_pos(findExchReact(model))] <- F
	}
	pb <- txtProgressBar(min=0, max=react_num(model), style=3)
	while(any(toCheck)){
		setTxtProgressBar(pb, react_num(model)-sum(toCheck))
#		rid <- sample(which(toCheck), 1)
		rid <- which(toCheck)[1]
		max <- optimizeProb(sba, react=rid, obj_coef=1, lpdir="max")
		if(length(checkSolStat(max$stat))!=0){
			if(SYBIL_SETTINGS("SOLVER") == "cplexAPI"){
				if(max$stat == CPX_STAT_OPTIMAL_INFEAS){
					warning(paste0("Status was:", getMeanStatus(max$stat)))
				}
			}else{
				stop(paste0("Solution is not optimal", getMeanStatus(max$stat)))
			}
		}
		toCheck[toCheck][abs(max$fluxes[toCheck]) > tol] <- F
		min <- optimizeProb(sba, react=rid, obj_coef=1, lpdir="min")
		if(length(checkSolStat(min$stat))!=0){
			if(SYBIL_SETTINGS("SOLVER") == "cplexAPI"){
				if(min$stat == CPX_STAT_OPTIMAL_INFEAS){
					warning(paste0("Status was:", getMeanStatus(min$stat)))
				}
			}else{
				stop(paste0("Solution is not optimal", getMeanStatus(min$stat)))
			}
		}
		toCheck[toCheck][abs(min$fluxes[toCheck]) > tol] <- F
		
		
		if(abs(max$obj) < tol && abs(min$obj) < tol ){
			blocked[rid] <- T
		}
		
		toCheck[rid] <- F
	}
	close(pb)
	return(blocked)
}

#irrevBlockedReact <- function(model, exex=TRUE, tol=SYBIL_SETTINGS("TOLERANCE")){
#	stopifnot(is(model, "modelorg_irrev"))
#	
#	sba <- sysBiolAlg(model, algorithm="fba")
#	toCheck <- rep(T, react_num(model))
##	blocked <- rep(F, react_num(model))
#	
#	if(exex){
#		toCheck[react_pos(findExchReact(model))] <- F
#	}
#	while(any(toCheck)){
#		print(sum(toCheck))
#		checkVector <- rep(0, react_num(model))
#		checkVector[toCheck] <- 1
#		fba <- optimizeProb(sba,
#							react=1:(1*react_num(model)),
#							obj_coef=c(checkVector))
#		stopifnot(length(checkSolStat(fba$stat))==0)
#		
#		f <- fba$fluxes
#		toCheck[toCheck][abs(f[toCheck]) > tol] <- F
##		blocked[toCheck][abs(f[toCheck]) > tol] <- F
#		
#		print(fba[c("obj", "stat")])
#		if(abs(fba$obj) < tol){
#			break
#		}
#	}
#	return(toCheck)
#}



#superFastBlockedReact <- function(model, tol=SYBIL_SETTINGS("TOLERANCE"), exex=TRUE){
#	stopifnot(is(model, "modelorg"))
#	
#	sba <- sysBiolAlg(model, algorithm="mtf", wtobj=0, lpdir="max")
#	
#	toCheck <- rep(T, react_num(model))
##	blocked <- rep(F, react_num(model))
#	
#	if(exex){
#		toCheck[react_pos(findExchReact(model))] <- F
#	}
#	while(any(toCheck)){
#		print(sum(toCheck))
#		checkVector <- rep(0, react_num(model))
#		checkVector[toCheck] <- 1
#		mtf <- optimizeProb(sba, fldind=FALSE,
#							react=1:(3*react_num(model)),
#							obj_coef=c(rep(0, react_num(model)), checkVector, checkVector))
#		stopifnot(length(checkSolStat(mtf$stat))==0)
#		
#		f <- mtf$fluxes[1:react_num(model)]
#		toCheck[toCheck][abs(f[toCheck]) > tol] <- F
##		blocked[toCheck][abs(f[toCheck]) > tol] <- F
#		
#		print(mtf[c("obj", "stat")])
#		browser()
#		if(abs(mtf$obj) < tol){
#			toCheck <- F
#		}
#		
#	}
#	return(toCheck)
#}

#library(sybil)
#data(Ec_core)
#source("helperScripts/fastBlockedReact.R")
#superFastBlockedReact(Ec_core)

#all(fastBlockedReact(Ec_core) == blockedReact(Ec_core))
#load("~/workspace/models/iAF1260.Rdata")
#data.frame(round(f), checkVector, round(mtf$fluxes[react_num(model)+(1:react_num(model))]))













