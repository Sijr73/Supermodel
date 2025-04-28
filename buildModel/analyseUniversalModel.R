#!/usr/bin/Rscript
library(methods)
library(sybil)

print(load("~/DATA/universalBiGG.ver1.Rdata"))
#for e.coli supermodel
#print(load("~/DATA/E.coli/universalBiGG.ver1.Rdata"))

source("~/helper/fastBlockedReact.R")

SYBIL_SETTINGS("SOLVER", "cplexAPI")
sp <- list(	CPX_PARAM_EPINT=1e-10, 		# tolerance to integer values
			CPX_PARAM_EPRHS=1e-9, 		# tolerance to bounds
			CPX_PARAM_TILIM=60*10, 		# timelimit in seconds
			CPX_PARAM_WORKMEM=1024*5, 	# workspace memory limit
			CPX_PARAM_TRELIM=1024*5, 	# tree memory limit
			CPX_PARAM_THREADS=7L, 		# only use 1 thread
			CPX_PARAM_PARALLELMODE=CPX_PARALLEL_OPPORTUNISTIC # no deterministic search
)

uni=uni1.2
cat(">>> resetting directions exchanges, demand reactions, sinks (EX, DM, SK)\n")
ex <- findExchReact(uni)
#ex=ex[-(1443:1444)]
stopifnot(all(grepl("^DM|^SK|^EX", react_id(ex))))

lowbnd(uni)[react_pos(ex)[grep("^DM", react_id(ex))]] <- 0
uppbnd(uni)[react_pos(ex)[grep("^DM", react_id(ex))]] <- 1000

lowbnd(uni)[react_pos(ex)[grep("^SK", react_id(ex))]] <- 0
uppbnd(uni)[react_pos(ex)[grep("^SK", react_id(ex))]] <- 1000

lowbnd(uni)[react_pos(ex)[grep("^EX", react_id(ex))]] <- -1000
uppbnd(uni)[react_pos(ex)[grep("^EX", react_id(ex))]] <- 1000

cat(">>> updating reaction reversibility slot\n")

react_rev(uni) <- sapply(1:react_num(uni), function(i){
	if(lowbnd(uni)[i] == -1000 && uppbnd(uni)[i] == 1000 ){
		return(T)
	}
	if(lowbnd(uni)[i] == 0 && uppbnd(uni)[i] == 1000 ){
		return(F)
	}
	if(lowbnd(uni)[i] == -1000 && uppbnd(uni)[i] == 0 ){
		return(F)
	}
	print(i)
	print(c(lowbnd(uni)[i], uppbnd(uni)[i]))
	stop("unknown direction")
})



ex <- findExchReact(uni)
lowbnd(uni)[react_pos(ex)] <- -1000
uppbnd(uni)[react_pos(ex)] <- 1000

if(!exists("bl")){
	cat(">>> calculating blocked reactions\n")
	bl <- fastBlockedReact(uni)
	save(bl, file="~/DATA/blockedReact.Rdata")
	cat(">>> blocked reactions calculated.\n")
}

rmReactId <- react_id(uni)[bl]


lowbnd(uni)[react_pos(ex)] <- lowbnd(ex)
uppbnd(uni)[react_pos(ex)] <- uppbnd(ex)
uniBlkRm <- rmReact(uni, which(bl))
ex <- findExchReact(uniBlkRm)
obj_coef(uniBlkRm) <- rep(0, react_num(uniBlkRm))

modelBiomassMap <- lapply(modelBiomassMap, function(x){
	x[x %in% react_id(uniBlkRm)]
})
modelReactMap <- lapply(modelReactMap, function(x){
	x[x %in% react_id(uniBlkRm)]
})

uni1.1 <- uniBlkRm
cat(">>> removing duplicated reactions\n")

if(!exists("dr")){
	dr <- doubleReact(uni1.1, checkRev=T, linInd=T)
}
skipped <- list()
toRemove <- integer(0)
for(p in dr){
	if(length(unique(lowbnd(uni1.1)[p]))>1){
		print(p)
		skipped <- append(skipped, list(p))
		next
	}
	if(length(unique(uppbnd(uni1.1)[p]))>1){
		print(p)
		skipped <- append(skipped, list(p))
		next
	}
	
	p <- react_id(uni1.1)[p]
	
	toRemove <- append(toRemove, p[-1])
	modelReactMap <- lapply(modelReactMap, function(x){
		if(any(p[-1] %in% x)){
			return(union(p[1], setdiff(x, p[-1])))
		}
		return(x)
	})
}
cat(">>> skipping these pairs due different reaction direction.\n")
print(skipped)
source("~/buildModel/addDefaultBOF.R")
toRemove <- setdiff(toRemove, c(grep("^FCC_", react_id(uni)), "ENERGY_BOF", "GENERAL_BOF"))
uni1.1 <- rmReact(uni1.1, toRemove)
modelReactMap <- lapply(modelReactMap, function(x) union(x, c("ENERGY_BOF", "GENERAL_BOF")))



ex <- findExchReact(uni1.1)
exPos <- react_pos(ex)[grep("_(LR|RL|B)$", react_id(ex), value=F)]
exId <- grep("_(LR|RL|B)$", react_id(ex), value=T)
react_id(uni1.1)[exPos] <- gsub("_(LR|RL|B)$", "", react_id(uni1.1)[exPos])


modelReactMap <- lapply(modelReactMap, function(x){
	x[x %in% exId] <- gsub("_(LR|RL|B)$", "", x[x %in% exId])
	return(x)
})

stopifnot(anyDuplicated(react_id(uni1.1)) == FALSE)




uniFCC <- uni1.1
uniFCC <- changeObjFunc(uniFCC, grep("^FCC_", react_id(uniFCC)))
lowbnd(uniFCC)[react_pos(ex)] <- 0

uniFCCSba <- sysBiolAlg(uniFCC, algorithm="fba", solverParm=sp)


print(optimizeProb(uniFCCSba)["obj"])

fccReact <- grep("^FCC_", react_id(uniFCC), value=T)
fccResult <- lapply(modelReactMap, function(reacts){
	on <- union(reacts, fccReact)
	onId <- match(on, react_id(uniFCC))
	onId <- na.omit(onId)
	offId <- setdiff(1:react_num(uniFCC), onId)
	optimizeProb(uniFCCSba, react=offId, lb=rep(0, length(offId)), ub=rep(0, length(offId)))
})

# if all these conditions are true, single models dont have cycles


stopifnot(all(sapply(fccResult, "[[", "ok")==0))
stopifnot(all(sapply(fccResult, "[[", "stat")==1))
#stopifnot(all(sapply(fccResult, "[[", "obj")==0))

fcc <- match(fccReact, react_id(uniFCC))
a <- sapply(fccResult, function(x) x$fluxes[fcc])
b <- apply(a, 2, function(x) fccReact[round(x, digits=6) != 0])

cat(">>> these models have energy cycles\n")
print(b[sapply(b, length) > 0 ])


biomassReacts <- unique(unlist(modelBiomassMap))
uniSba <- sysBiolAlg(uni1.1, algorithm="fba", solverParm=sp)

biomassResults <- lapply(names(modelBiomassMap), function(n){
	r <- modelBiomassMap[[n]]
	mid <- match(modelReactMap[[n]], react_id(uni1.1))
	mid <- setdiff(1:react_num(uni1.1), mid)
	
	id <- match(r, react_id(uni1.1))
	id=na.omit(id)
	res <- lapply(id, function(singleid) optimizeProb(uniSba, 
							react=c(singleid, mid),
							lb=rep(0, length(mid)+1), 
							ub=c(1000, rep(0,length(mid))), 
							obj_coef=c(1, rep(0, length(mid))))[c("stat", "obj")])
	for(ri in res){
		if(length(checkSolStat(ri$stat))!=0){
			print(n)
			print(id)
			print(ri$stat)
#			browser()
			stop()
		}
	}
	sapply(res, "[[", "obj")
})
names(biomassResults) <- names(modelBiomassMap)


biomassReactionTable <- lapply(names(modelBiomassMap[sapply(modelBiomassMap, length)!=0]), function(x) data.frame(modelBiomassMap[[x]], x, length(modelBiomassMap[[x]])))
biomassReactionTable <- do.call(rbind, biomassReactionTable)

colnames(biomassReactionTable) <- c("react_id", "modelName", "numberBiomassReactions")
biomassReactionTable$use <- "X"
biomassReactionTable$reactionName <- react_name(uni)[match(biomassReactionTable$react_id, react_id(uni))]
write.csv(biomassReactionTable, file="~/DATA/selectBiomassTable.csv", row.names=F, quote = F)
#write.csv(biomassReactionTable, file="~/DATA/E.coli/selectBiomassTable.csv", row.names=F, quote = F)

################################################################################

# Manual Step

if(!file.exists("selectBiomassTable_mod.csv")){
	stop("select biomass reactions for each model by placing a X and save selectBiomassTable.csv as selectBiomassTable_mod.csv")
}

################################################################################


selectedBiomasses <- read.table("~/DATA/selectBiomassTable.csv", sep=",", header=T, stringsAsFactors=F)
#selectedBiomasses$use[188]="X"
stopifnot(all(selectedBiomasses$use %in% c("X", "")))
l <- split(selectedBiomasses[c("use", "react_id")], selectedBiomasses$modelName)
biomassReactions <- sapply(l, function(x) with(x, react_id[use=="X"]))

modelBiomassMapNew <- rep(NA, length(modelBiomassMap))
names(modelBiomassMapNew) <- names(modelBiomassMap)

modelBiomassMapNew[names(biomassReactions)] <- biomassReactions

r <- unlist(modelBiomassMap)[!unlist(modelBiomassMap) %in% modelBiomassMapNew]
r <- match(r, react_id(uni1.1))
lowbnd(uni1.1)[r] <- uppbnd(uni1.1)[r] <- 0
lowbnd(uniFCC)[r] <- uppbnd(uniFCC)[r] <- 0







modelBiomassMapSelection <- modelBiomassMapNew[names(modelBiomassMap)]
save(file="~/DATA/universalBiGG.ver1.1.Rdata", uni1.1, modelReactMap, modelBiomassMap, modelBiomassMapSelection)
#save(file="~/DATA/E.coli/universalBiGG.ver1.1.Rdata", uni1.1, modelReactMap, modelBiomassMap, modelBiomassMapSelection)


# removing of compartments not necessary:

##reactions in compartments:

metCompMatrix <- sapply(seq(along=mod_compart(uni1.1)), function(x){
	return(met_comp(uni1.1)==x)
})

sbin <- S(uni1.1) != 0


reactionCompMatrix <- t(metCompMatrix) %*% sbin
reactionsPerComp <- rowSums(reactionCompMatrix!=0)
names(reactionsPerComp) <- mod_compart(uni1.1)

keepComp <- c("c", "e", "p")
removeComp <- setdiff(mod_compart(uni1.1), keepComp)
rownames(reactionCompMatrix) <- mod_compart(uni1.1)
removeCompReact <- which(colSums(reactionCompMatrix[removeComp, ])>0)

uniFCC <- rmReact(uniFCC, react=react_id(uni1.1)[removeCompReact])
uni1.1 <- rmReact(uni1.1, react=react_id(uni1.1)[removeCompReact])
save(file="New analysis bacillus/DATA/universalBiGG.ver1.1_compRemoved.Rdata", uniFCC, uni1.1, modelReactMap, modelBiomassMap, modelBiomassMapSelection)



uni1.1=uni1.2
r <- lapply(modelReactMap, function(x) 1:react_num(uni1.1))
off <- lapply(modelReactMap, function(x) !react_id(uni1.1) %in% x)
u <- lapply(off, function(x) ifelse(x, rep(0, react_num(uni1.1)), uppbnd(uni1.1)))
l <- lapply(off, function(x) ifelse(x, rep(0, react_num(uni1.1)), lowbnd(uni1.1)))
o <- lapply(modelBiomassMapSelection, function(x){
	v <- rep(0, react_num(uni1.1))
	if(is.na(x)){
		return(v)
	}
	v[grep(x, react_id(uni1.1))] <- 1
	return(v)
})

fbaResult <- optimizer(uni1.1, react=r, ub=u, lb=l, obj_coef=o,algorithm = "fba",fld="all")
# / remove of compartments





