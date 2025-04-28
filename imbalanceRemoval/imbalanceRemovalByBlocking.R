#!/usr/bin/Rscript
library(methods)
library(parallel)
library(sybil)
library(sybilSWITCH)
SYBIL_SETTINGS("SOLVER", "cplexAPI")
source("~/energyGeneratingCycleRemoval/sysBiolAlg_armClass.R")


print(load("~/DATA/universalBiGG.ver1.1.Rdata"))
#Ecoli supermodel
#print(load("~/DATA/E.coli/universalBiGG.ver1.1.Rdata"))

#print(load("~/DATA/E.coli/atomsMatrix.Rdata"))
print(load("~/DATA/atomsMatrix.Rdata"))

#print(load("~/DATA/E.coli/wronglyBalancedReact.Rdata"))
print(load("~/DATA/wronglyBalancedReact.Rdata"))

atomsMatrix <- atomsMatrix[!is.na(atomsMatrix[,1]),]
uni <- uni1.1

pos <- match(wronglyBalancedReact, react_id(uni))
lowbnd(uni)[pos] <- uppbnd(uni)[pos] <- 0
uniAtom <- uni

for(i in colnames(atomsMatrix)){
	#cavetas: lower bound is set to 0 because we only want to consider atom creation. if a atom vanishes we dont care.
	uniAtom <- addReact(uniAtom, id=paste0("CHECK_ATOM_", i, "_REACTION"), Scoef=-1, met=paste0("ATOM_", i), reversible=T, lb=0, ub=1000, obj=0)
}
################################################################################
	threshold <- 0.01
########################################################################
	testReact <- grep("CHECK_ATOM_.+_REACTION", react_id(uniAtom))
	biomassReact <- na.omit(unlist(modelBiomassMap))
	biomassReact <- unique(c(biomassReact, "GENERAL_BOF", "ENERGY_BOF"))
################################################################################
#add Atom metaboltes to atom matrix
for(i in colnames(atomsMatrix)){
	atomsMatrix <- rbind2(atomsMatrix,rep(0, ncol(atomsMatrix)))
	rownames(atomsMatrix)[nrow(atomsMatrix)] <- paste0("ATOM_", i)
	atomsMatrix[paste0("ATOM_", i), i] <- 1
}
#solver config
node <- Sys.info()["nodename"]
p_list <- NULL
if(grepl("^jarvis$", node)){
	cat(">>> config for jarvis\n")
	sp <- suggestedArmSolverSettings(threads=7, timelimit=10*60, workMemLimit=20, treeMemLimit=20)
}
if(grepl("^jedi|sith$", node)){
	cat(">>> config for jedi|sith\n")
	sp <- suggestedArmSolverSettings(threads=32, timelimit=72*3600, workMemLimit=200, treeMemLimit=200)
}
if(grepl("^hilbert", node)){
	cat(">>> config for hilbert\n")
	sp <- suggestedArmSolverSettings(threads=12, timelimit=7*24*3600, workMemLimit=100, treeMemLimit=100)
}
if(is.null(sp)){
	print(node)
	stop("no config for this node!!")
}
################################################################################

if(any(!biomassReact %in% react_id(uni))){
	stop("biomassreactions given, that are not present in the model!")
	biomassReact <- biomassReact[biomassReact %in% react_id(uni)]
}

biomassReact <- match(biomassReact, react_id(uniAtom))
knownMets <- (rownames(atomsMatrix))
ex <- findExchReact(uniAtom)
exmet <- met_id(ex)
exReactOff <- react_pos(ex)[!exmet %in% knownMets]
cat(">>> insufficient data for exchange reactions with these mets\n")
print(met_id(ex)[!exmet %in% knownMets])
lowbnd(uniAtom)[exReactOff] <- uppbnd(uniAtom)[exReactOff] <- 0

bioReactOff <- biomassReact[sapply(biomassReact, function(r){
  s <- shrinkMatrix(uniAtom, j=r)
	metIds <- rownames(s)
	notKnown <- metIds[!metIds %in% knownMets]
	if(length(notKnown)>0){
		cat(">>> no data for these Metabolites:\n")
		print(notKnown)
	}
	return(length(notKnown)>0)
})]

cat(">>> insufficient data for these biomass reactions\n")
print(react_id(uniAtom)[bioReactOff])
lowbnd(uniAtom)[bioReactOff] <- uppbnd(uniAtom)[bioReactOff] <- 0

# these reactions are also in/outflux of the model.
inOutReact <- setdiff(c(react_pos(ex), biomassReact), c(exReactOff, bioReactOff))
# met X react
s <- shrinkMatrix(uniAtom, j=inOutReact)
# met X atoms
ac <- atomsMatrix[rownames(s), ]
# react X atoms
reactBalance <- t(s) %*% ac


uniAtom <- changeObjFunc(uniAtom, testReact)

for(cA in colnames(atomsMatrix)){
	S(uniAtom)[which(met_id(uniAtom)==paste0("ATOM_", cA)), inOutReact] <- reactBalance[, cA]
}

source("~/helper/fastBlockedReact.R")
uniAtomBlocked <- uniAtom
lowbnd(uniAtomBlocked)[testReact] <- uppbnd(uniAtomBlocked)[testReact] <- 0

ab <- fastBlockedReact(uniAtomBlocked)
atomBlocked <- react_id(uniAtomBlocked)[ab]
print(length(atomBlocked))

uniAtom2 <- rmReact(uniAtom, atomBlocked)
print(opt <- optimizeProb(uniAtom2))
imbalances <- round(fluxes(opt)[obj_coef(uniAtom2)==1, ], digits=6)
names(imbalances) <- react_id(uniAtom2)[obj_coef(uniAtom2)==1]
cat(">>>following imbalances still exist\n")
print(imbalances)

imbalances["CHECK_ATOM_H_REACTION"] <- 0
if(!all(imbalances==0)){
	stop("there are imbalances other than hydrogen")
}
bioReactOff <- react_id(uniAtom)[bioReactOff]
exReactOff <- react_id(uniAtom)[exReactOff]
save(atomBlocked, exReactOff, bioReactOff, uniAtom, file="~/DATA/atomBlocked.Rdata")
#save(atomBlocked, exReactOff, bioReactOff, uniAtom, file="~/DATA/E.coli/atomBlocked.Rdata")







toremove <- unique(c(wronglyBalancedReact, atomBlocked, exReactOff, bioReactOff))
cat(">>> removing reactions because auf mass imbalance\n")
print(length(toremove))
cat(paste(toremove, collapse="\n"), file="~/DATA/removed_ver1.1_to_ver1.2.txt")
cat(paste(toremove, collapse="\n"), file="~/DATA/E.coli/removed_ver1.1_to_ver1.2.txt")

uni1.2 <- rmReact(uni1.1, toremove)

modelReactMap <- lapply(modelReactMap, function(x) x[x %in% react_id(uni1.2)])
modelBiomassMap <- lapply(modelBiomassMap, function(x) x[x %in% react_id(uni1.2)])
modelBiomassMapSelection[!modelBiomassMapSelection %in% react_id(uni1.2)] <- NA
save( list=c("uni1.2", "modelReactMap", "modelBiomassMap", "modelBiomassMapSelection"), file="~/DATA/universalBiGG.ver1.2.Rdata")
#save( list=c("uni1.2", "modelReactMap", "modelBiomassMap", "modelBiomassMapSelection"), file="~/DATA/E.coli/universalBiGG.ver1.2.Rdata")

cat(">>> EOF\n")








































