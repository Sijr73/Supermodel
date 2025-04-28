#!/usr/bin/Rscript
library(methods)
library(sybil)
library(parallel)


print(load("~/DATA/checkedModels.Rdata"))
##For E. coli supermodel
#print(load("~/DATA/E.coli/checkedModels.Rdata"))

source("~/helper/loadPropTables.R")

uni <- modelorg(id = "uni1", name = "BiGG Universal Model Ver 1")

subSys(uni) <- Matrix(F, nrow=0, ncol=length(allModels), dimnames=list(NULL, names(allModels)))

mod_compart(uni) <- unique(unlist(lapply(allModels, mod_compart)))


for(i in names(allModels)){
	print(i)
	print(react_num(uni))
	m <- allModels[[i]]
	for(j in 1:react_num(m)){
		uid <- match(react_id(m)[j], react_id(uni))
		if(is.na(uid)){
			# add this reaction
			s <- S(m)[,j]
			sid <- which(s!=0)
			s <- s[sid]
			
			
			uni <- addReact(uni,
				id=react_id(m)[j],
				met=met_id(m)[sid],
				Scoef=s,
				reversible=ifelse(lowbnd(m)[j]<0, T, F),
				lb=lowbnd(m)[j],
				ub=uppbnd(m)[j],
				obj=0,
				subSystem=NA,
				gprAssoc="",
				reactName=react_name(m)[j],
				metName="",
				metComp=mod_compart(m)[met_comp(m)[sid]]
			)
			
		}
		else{
			# check if already added reaction is equal
			
#			stopifnot(react_rev(uni)[uid]==react_rev(m)[j])
			
			stopifnot(lowbnd(uni)[uid]==lowbnd(m)[j])
			stopifnot(uppbnd(uni)[uid]==uppbnd(m)[j])
			
			s1 <- shrinkMatrix(uni, j=uid)
			s2 <- shrinkMatrix(m, j=j)
			
			stopifnot(all(dim(s1)==dim(s2)))
			
			s1 <- s1[order(rownames(s1)),]
			s2 <- s2[order(rownames(s2)),]
			
			stopifnot(all(s1==s2))
			
		}
	}
}

# remove empty compartments:

compartmentCount <- table(met_comp(uni))
names(compartmentCount) <- mod_compart(uni)[as.integer(names(compartmentCount))]

metComp <- mod_compart(uni)[met_comp(uni)]
mod_compart(uni) <- mod_compart(uni)[mod_compart(uni) %in% names(compartmentCount)]

met_comp(uni) <- match(metComp, mod_compart(uni))


# ensure exchange reactions in the correct direction (they should always be -1)

ex <- findExchReact(uni)
exToFix <- sapply(react_pos(ex), function(j){
	i <- which(S(uni)[,j]!=0)
	stopifnot(length(i) == 1)
	return(S(uni)[i, j] == 1)
})

S(uni)[, react_pos(ex)[exToFix]] <- S(uni)[, react_pos(ex)[exToFix]] * -1

tmp <- lowbnd(uni)[react_pos(ex)[exToFix]]
lowbnd(uni)[react_pos(ex)[exToFix]] <- uppbnd(uni)[react_pos(ex)[exToFix]] * -1
uppbnd(uni)[react_pos(ex)[exToFix]] <- tmp * -1
# end





met_name(uni) <- metProp[match(gsub("\\[.+\\]", "", met_id(uni)), metProp$universal_bigg_id), "name"]


print(load("Downloads/freeCofactorsBiGG.Rdata"))

for(n in names(freeCofactors)){
	f <- freeCofactors[[n]]
	f$metIds <- gsub("glu-L", "glu__L", f$metIds)
	
	uni <- addReact(uni, id=paste0("FCC_", n), met=f$metIds, Scoef=f$Svalues, lb=0, ub=1000, reversible=F, reactName=paste0("FCC_", n))
}

modelReactMap <- lapply(allModels, react_id)

modelBiomassMap <- lapply(allModels, function(m) react_id(m)[(obj_coef(m)!=0)])



source("~/buildModel/addDefaultBOF.R")


ex <- findExchReact(uni)

#exSplit <- split(react_id(ex), react_id(uni)[react_pos(ex)])
#modelReactMapNew <- lapply(modelReactMap, function(x){
#	for(y in names(exSplit)){
#		if(any(exSplit[[y]] %in% x)){
#			x <- setdiff(x, exSplit[[y]])
#			x <- union(x, y)
#		}
#	}
#	return(x)
#})

#exNew <- findExchReact(uni)

#exSplitNew <- split(react_pos(exNew), react_id(exNew))
#exSplitNew[sapply(exSplitNew, length)>1]

#rmEx <- unique(unlist(lapply(exSplitNew, "[", -1)))

#uni <- rmReact(uni, rmEx)

stopifnot( all(sapply(modelReactMap, function(x) all(x %in% react_id(uni)))))

save(uni, modelReactMap, modelBiomassMap, file="~/DATA/universalBiGG.ver1.Rdata")

#save(uni, modelReactMap, modelBiomassMap, file="~/DATA/E.coli/universalBiGG.ver1.Rdata")






























