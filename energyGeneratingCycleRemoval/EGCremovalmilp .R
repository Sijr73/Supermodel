#!/usr/bin/Rscript

library(methods, verbose=T)
library(sybil, verbose=T)
library(cplexAPI, verbose=T)
#library(sybilSWITCH, verbose=T)
#source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")
source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")

SYBIL_SETTINGS("SOLVER", "cplexAPI")

if(!interactive()){
	args <- commandArgs(trailingOnly = TRUE)
	calcOnly <- as.integer(args[1])
}else{
	calcOnly <- 1
}
print(load("~/DATA/universalBiGG.ver1.2.Rdata"))
#print(load("~/EGC/universalBiGG.ver1.2.Rdata"))
uni1.1 <- uni1.2
uni <- uni1.1
ex <- findExchReact(uni)
modelReactMapId <- lapply(modelReactMap, match, table=react_id(uni))

stopifnot(calcOnly <= length(modelReactMap))
outputFile <- paste0("~/EGC/resultsmilp/", names(modelReactMap)[calcOnly],".Rdata")
if(file.exists(outputFile)){
	cat(">>> file already existing, nothing to calculat\n")
	quit(save="no")
}

egcOff <- scan(text="CAT,DHPTDNR,DHPTDNRN,FHL,SPODM,SPODMpp,SUCASPtpp,SUCFUMtpp,SUCMALtpp,SUCTARTtpp", sep=",", what="character")
egcOff <- unlist(lapply(egcOff, function(x) grep(paste0("^", x, "(_.+)*$"), react_id(uni1.1), value=T, ignore.case=T)))
egcOffId <- na.omit(match(egcOff, react_id(uni1.1)))

#print(load("globalFitSuggestions.Rdata"))
#gfsId <- unique(unlist(lapply(globalFitSuggestions, match, table=react_id(uni))))

fcc <- grep("^FCC_", react_id(uni))

lowbnd(uni)[react_pos(ex)] <- 0
uni <- changeObjFunc(uni, fcc)

#fccFba <- optimizeProb(uni, retOptSol=F)

c <- integer(0)
while(T){
	if(length(c)==0){
		erg <- optimizeProb(uni, solverParm=list(CPX_PARAM_EPRHS=1e-8), retOptSol=F)
	}
	else{
		erg <- optimizeProb(uni, react=c, lb=rep(0, length(c)), ub=rep(0, length(c)), solverParm=list(CPX_PARAM_EPRHS=1e-8), retOptSol=F)
	}
	b <- which(abs(erg$fluxes)> 1e-6)
	b <- setdiff(b, fcc)
	print(length(b))
	if(length(b)==0){
		break
	}
	c <- union(c, b)
}
egcParticipatingReactions <- c


fccResult <- lapply(modelReactMapId, function(reacts){
	onId <- union(reacts, fcc)
	offId <- setdiff(1:react_num(uni), onId)
	optimizeProb(uni, react=offId, lb=rep(0, length(offId)), ub=rep(0, length(offId)), retOptSol=F)
})
egcModelOrder <- order(sapply(fccResult, "[[", "obj"), decreasing=F)


sp <- suggestedArmSolverSettings(threads=8, timelimit=120, workMemLimit=10, treeMemLimit=10)
sp$CPX_PARAM_EPINT <- 1e-9
sp$CPX_PARAM_EPRHS <- 1e-9
sp$CPX_PARAM_MIPEMPHASIS <- CPX_MIPEMPHASIS_FEASIBILITY
poCmd = list(c("getMIPrelGapCPLEX", "LP_PROB@oobj@env", "LP_PROB@oobj@lp"),
			 c("getBestObjValCPLEX", "LP_PROB@oobj@env", "LP_PROB@oobj@lp"),
			 c("getObjValCPLEX", "LP_PROB@oobj@env", "LP_PROB@oobj@lp")
		)




activeReactions <- integer(0)
excludeReactions <- fcc

cycles <- list()
removedFwd <- list()
removedBwd <- list()

for(basemodel in names(modelReactMap)[calcOnly]){
	cat(">>> basemodel is:\n")
	print(basemodel)
	cycles[[basemodel]] <- list()
	removedFwd[[basemodel]] <- list()
	removedBwd[[basemodel]] <- list()
	order <- c(basemodel, setdiff(names(modelReactMap)[egcModelOrder], basemodel))
	
	for(i in seq(along=order)){
		m <- order[i]
		cat(">>> now adding model\n")
		print(m)
		activeReactions <- unique(c(unlist(modelReactMapId[order[1:i]]), fcc))
		
		cat(">>> building problem Object\n")
		uniArm <- sysBiolAlg(uni, biomassThreshold=1, additionalReact=1:react_num(uni), algorithm="arm", solverParm=sp)
		cat(">>> problem Object built\n")
		
		cycles[[c(basemodel, m)]] <- list()
		removedBwd[[c(basemodel, m)]] <- integer(0)
		removedFwd[[c(basemodel, m)]] <- integer(0)
		egcPresent <- T
		modelReactions <- modelReactMap[[m]]
		
		while(egcPresent){
			off <- setdiff(1:react_num(uni), activeReactions)
			offFwd <- union(unlist(removedFwd[[c(basemodel)]]), off)
			offBwd <- union(unlist(removedBwd[[c(basemodel)]]), off)
			off <- union(offFwd, offBwd)
			
			lb <- ifelse(off %in% offBwd, rep(0, length(off)), lowbnd(uni)[off])
			ub <- ifelse(off %in% offFwd, rep(0, length(off)), uppbnd(uni)[off])
			
			print(system.time(arm <- optimizeProb(uniArm, react=off, lb=lb, ub=ub, poCmd=poCmd)))
			
			if(arm$stat == CPXMIP_INFEASIBLE){
				break
			}
			if(!arm$stat %in% c(CPXMIP_OPTIMAL, CPXMIP_OPTIMAL_TOL, CPXMIP_TIME_LIM_FEAS, CPXMIP_OPTIMAL_INFEAS)){
				print(arm[setdiff(names(arm), "fluxes")])
				stop("problem solving")
			}
			cat(paste0("REL_GAP: ", arm$postP@pa[[c(1, 1)]], "\n"))
			cat(paste0("BEST_SOL: ", arm$postP@pa[[c(2, 1)]], "\n"))
			cat(paste0("OBJ_VAL: ", arm$postP@pa[[c(3, 1)]], "\n"))
			cat(paste0("ABS_GAP: ", arm$postP@pa[[c(3, 1)]] - arm$postP@pa[[c(2, 1)]] , "\n"))
			
			cycleId <- which(abs(getArmReactionFluxes(uni, arm)[1:react_num(uni)]) > 1e-6)
			cycle <- getArmReactionFluxes(uni, arm)[cycleId]
			names(cycle) <- react_id(uni)[cycleId]
			cat(">>> cycle:\n")
			print(cycle)
			
			rmR <- setdiff(cycleId, fcc)
			if(m != basemodel){
				rmR <- setdiff(rmR, modelReactMapId[[basemodel]])
			}
			if(any(rmR %in% egcOffId)){
				cat(">>> filtering for iJO paper suggestions\n")
				rmR <- intersect(rmR, egcOffId)
			}
#			if(any(rmR %in% gfsId)){
#				cat(">>> filtering for egcPaper suggestions\n")
#				rmR <- intersect(rmR, gfsId)
#			}
#			if(any(rmR %in% egcParticipatingReactions)){
#				rmR <- intersect(rmR, egcParticipatingReactions)
#			}
			
			stopifnot(length(rmR)>=1)
			rmR <- rmR[1] # because the solving process is non-deterministic, the first index is some kind of random, too.
			cycles[[c(basemodel, m)]] <- append(cycles[[c(basemodel, m)]], list(cycle))
			if(cycle[react_id(uni)[rmR]] > 0){
				cat(">>> removing forward:\n")
				print(react_id(uni)[rmR])
				removedFwd[[c(basemodel, m)]] <- append(removedFwd[[c(basemodel, m)]], rmR)
			}else{
				cat(">>> removing backward:\n")
				print(react_id(uni)[rmR])
				removedBwd[[c(basemodel, m)]] <- append(removedBwd[[c(basemodel, m)]], rmR)
			}
		}
	}
}

save(cycles, removedFwd, removedBwd, file=outputFile)






























