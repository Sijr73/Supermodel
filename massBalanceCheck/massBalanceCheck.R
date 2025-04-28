#!/usr/bin/Rscript
library(methods)
library(parallel)
library(sybil)
library(CHNOSZ)
library(rjson)
SYBIL_SETTINGS("SOLVER", "cplexAPI")


print(load("~/DATA/universalBiGG.ver1.1.Rdata"))
#for ecoli supermodel
#print(load("~/DATA/E.coli/universalBiGG.ver1.1.Rdata"))
metsInModel <- met_id(uni1.1)
print(load("~/DATA/formulaDFByMet.Rdata"))
#print(load("~/DATA/E.coli/formulaDFByMet.Rdata"))
stopifnot(all(metsInModel %in% names(metaboliteFormula)))

sumFormulaUnparsed <- metaboliteFormula[metsInModel]
#sumFormulaUnparsed[7434]='C46H66O2'


cat(">>> parse chemical formula\n")
data(thermo)
thermo()
options(warn=2)
sumFormula <- mclapply(sumFormulaUnparsed, function(x){
		tryCatch(makeup(x), error=function(e) NULL)
})
options(warn=1)


elements <- unique(unlist(lapply(sumFormula, names)))

atomsTable <- lapply(sumFormula, function(x){
	if(is.null(x)){
		return(rep(NA, length(elements)))
	}
	line <- rep(0, length(elements))
	names(line) <- elements
	
	line[names(x)] <- x
	return(line)
})

atomsTable <- do.call(rbind, atomsTable)
atomsMatrix <- Matrix(atomsTable)

rownames(atomsMatrix) <- metsInModel



save(atomsMatrix, file="~/DATA/E.coli/atomsMatrix.Rdata")
save(atomsMatrix, file="~/DATA/atomsMatrix.Rdata")


cat(">>>  validate mass balances\n")

uni <- uni1.1
source("~/buildModel/addDefaultBOF.R")
ex <- findExchReact(uni)
modelBiomassMap=modelBiomassMap[names(modelBiomassMap) %in%  "iAF987" == FALSE]   
exMets <- union(met_id(ex), rownames(shrinkMatrix(uni, j=c("GENERAL_BOF", unlist(modelBiomassMap)))))
exReact <- union(react_id(ex), c("GENERAL_BOF", unlist(modelBiomassMap)))
invalidExMets <- exMets[is.na(atomsMatrix[exMets,1])]

# deactivate exchanges with unknown metabolite composition
pos <- react_pos(ex)[met_id(ex) %in% invalidExMets]
invalidExReact <- react_id(ex)[met_id(ex) %in% invalidExMets]
lowbnd(uni)[pos] <- uppbnd(uni)[pos] <- 0

for(r in c("GENERAL_BOF", unlist(modelBiomassMap))){
	s <- shrinkMatrix(uni, j=r)
	
	if(any(rownames(s) %in% invalidExMets)){
		pos <- match(r, react_id(uni))
		lowbnd(uni)[pos] <- uppbnd(uni)[pos] <- 0
		invalidExReact <- union(invalidExReact, r)
	}
}
# end model preparing





# overall mass balance
reactBalance <- round(t(atomsMatrix[met_id(uni),]) %*% S(uni), digits=8)
colnames(reactBalance) <- react_id(uni)

reactIsBalancedKnown <- !is.na(reactBalance[1,])
names(reactIsBalancedKnown) <- react_id(uni)

unbalancedInModels <- do.call(rbind, lapply(modelReactMap, function(x){
	data.frame(sum(reactIsBalancedKnown[x]), length(x), sum(reactIsBalancedKnown[x])/length(x))
}))

Hcolumn <- which(rownames(reactBalance) == "H")
reactIsBalanced <- apply(reactBalance, 2, function(x){
	if(any(is.na(x))){
		return(NA)
	}
	return(all(x[-Hcolumn]==0))
})

print(table(reactIsBalanced[-match(exReact, react_id(uni))], useNA="always"))


balancedReactDF <- data.frame(id= react_id(uni), isEX=F, balanced=reactIsBalanced)
balancedReactDF$isEX[match(exReact, react_id(uni))] <- T

#balancedReactDF %>% filter(!isEX, balanced==F)

wronglyBalancedReact <- which(with(balancedReactDF, !isEX & balanced==F ))
cat(">>> number of reactions in models wrongly balanced:\n")
print(sapply(modelReactMap, function(x) sum(react_id(uni)[wronglyBalancedReact] %in% x)))
lowbnd(uni)[wronglyBalancedReact] <- uppbnd(uni)[wronglyBalancedReact] <- 0
wronglyBalancedReact <- react_id(uni)[wronglyBalancedReact]
print(wronglyBalancedReact)
save(wronglyBalancedReact, file="~/DATA/wronglyBalancedReact.Rdata")
#save(wronglyBalancedReact, file="~/DATA/E.coli/wronglyBalancedReact.Rdata")

if(F){
	cat("checking balance with FVA\n")
	fv <- fluxVar(uni, fld="fluxes")

	exReactS <- shrinkMatrix(uni, j=exReact)
	exReactBalance <- t(atomsMatrix[rownames(exReactS),]) %*% exReactS
	fluxes <- (fluxes(fv)[match(exReact, react_id(uni)), ])
	stopifnot(all(fluxes(fv)[match(invalidExReact, react_id(uni)), ]==0))

	fluxes[is.na(fluxes)] <- 0
	checkBalance <- exReactBalance %*% fluxes

	print(apply(checkBalance, 1, function(x) any(abs(x) > 1e-6)))

	if(any(abs(checkBalance)>1e-6)){
		stop(">>> mass blanace in the model not correct!\n")
	}
}





