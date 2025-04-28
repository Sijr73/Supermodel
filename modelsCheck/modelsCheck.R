#!/usr/bin/Rscript
library(methods)
library(sybil)
library(parallel)
library(CHNOSZ)
options(warn=1) # report warnings as they occur
source("~/modelsCheck/checkFunctions.R")

source("~/helper/loadPropTables.R")



SYBIL_SETTINGS("SOLVER", "cplexAPI")
modelsTable <- read.table("~/sourceData/model/modelList.tsv", sep="\t", header=T, stringsAsFactors=F, quote="", comment.char="")
#excludeModels <- scan("../sourceData/excludeModels.txt", what="character", comment.char = "#")

cat("::: excluding the following models from the project:\n")
#print(excludeModels)

#modelsTable <- modelsTable[!(modelsTable$bigg_id %in% excludeModels),]


stripCompartment <- function(x){
	gsub("\\[(\\w)\\]$", "_\\1", x)
}



print(load("~/DATA/models.Rdata"))
##for ecoli supermodel
#print(load("New analysis bacillus/DATA/E.coli/models.Rdata"))

allModels <- m
rm(m)
load("~/sourceData/buchnera.Rdata")
allModels$iSM199 <- buchnera
gc()
names(allModels) <- modelsTable$bigg_id
m<- allModels[modelsTable$bigg_id]
allModels <- m
#rm(m)
modelsTable


cat("::: fixing exchange reactions with coefficients != -1\n")

exchToFix <- mclapply(allModels, function(m){
	oneEntry <- which(apply(S(m)!=0, 2, sum) == 1)
	toFix <- oneEntry[colSums(abs(S(m)[,oneEntry, drop=F])) != 1]
	toFix
})

for(i in names(exchToFix[sapply(exchToFix, length)>0])){
	for(j in exchToFix[[i]]){
		S(allModels[[i]])[, j] <- S(allModels[[i]])[, j] / abs(sum(S(allModels[[i]])[, j]))
	}
}


cat("::: removing reactions blocked by bounds (lb==0 && ub==0)\n")
blockedByBound <- mclapply(allModels,
	function(m){
		react_id(m)[(uppbnd(m)==0 & lowbnd(m)==0)]
	}
)

allModels <- lapply(names(allModels),
	function(m){
		if(length(blockedByBound[[m]]) > 0){
			return(rmReact(allModels[[m]], blockedByBound[[m]], rm_met=TRUE))
		}
		else{
			return(allModels[[m]])
		}
	}
)
names(allModels) <- modelsTable$bigg_id

cat("::: finding exchange reactions\n")

ex <- mclapply(allModels, findExchReact)
options(warn=2) # report all warnings as error

exErr <- which(!sapply(ex, is, "reactId_Exch"))

if(length(exErr)!=0){
	print(exErr)
	stop("there are models without exchange reactions")
}



modelMedia <- lapply(ex, function(x) met_id(x)[lowbnd(x) < 0])
save(modelMedia, file="New analysis bacillus/DATA1/modelMedia.Rdata")


cat("::: settings bounds to either 1000 or -1000 \n")
for(m in names(allModels)){
	
	uppbnd(allModels[[m]])[uppbnd(allModels[[m]]) >  0] <- 1000
	lowbnd(allModels[[m]])[lowbnd(allModels[[m]]) >= 0] <- 0
	
	uppbnd(allModels[[m]])[uppbnd(allModels[[m]]) <= 0] <- 0
	lowbnd(allModels[[m]])[lowbnd(allModels[[m]]) <  0] <- -1000
}


cat("::: check if reactions are same in the models\n")

reactFromModel <- lapply(names(allModels), 
	function(x){
		m <- allModels[[x]]
		data.frame(model=x, id=react_id(m), lb=lowbnd(m), ub=uppbnd(m), rev=react_rev(m),
			 stringsAsFactors=F)
	}
)
reactFromModel <- do.call("rbind", reactFromModel)




#reactFromModel$idDir <- sapply(regmatches(reactFromModel$id, regexec("_(LR|RL|B)", reactFromModel$id)), "[", 2)
#reactFromModel$onlyId <- sapply(regmatches(reactFromModel$id, regexec("^(R.+)_", reactFromModel$id)), "[", 2)


reactFromModelByReact <- split(reactFromModel, reactFromModel$id)

cat("::: check if reactions have the same coefficients\n")

ism <- allModels$iSM199

for(i in react_id(allModels$iSM199)){
	
	mns <- reactFromModelByReact[[i]]$model
	if(length(mns) == 1){
		next # reaction is only in iSM
	}
	
	s1 <- shrinkMatrix(allModels[[setdiff(mns, "iSM199")[1] ]], j=i)
	s2 <- shrinkMatrix(ism, j=i)
	
	
	
	s1 <- s1[order(rownames(s1)),]
	s2 <- s2[order(rownames(s2)),]
	
	if((!all(dim(s1)==dim(s2))) || (!all(s1==s2))){
		print(s1)
		print(s2)
		
		if(all(s1==s2*-1)){
			S(ism)[,react_id(ism)==i] <- S(ism)[,react_id(ism)==i] * -1
			tmp <- lowbnd(ism)[react_id(ism)==i] * -1
			lowbnd(ism)[react_id(ism)==i] <- uppbnd(ism)[react_id(ism)==i] *-1
			uppbnd(ism)[react_id(ism)==i] <- tmp
		}else{
			react_id(ism)[react_id(ism)==i] <- paste0(i, "_iSM")
		}
	}
}

allModels$iSM199 <- ism




reactFromModelByReact <- split(reactFromModel, reactFromModel$id)


checkedReact <- sapply(reactFromModelByReact, function(x){
	if(length(unique(x$lb)) > 1){
		return(F)
	}
	if(length(unique(x$ub)) > 1){
		return(F)
	}
	if(length(unique(x$rev)) > 1){
		return(F)
	}
	return(T)
})



dirSuffix <- function(l, u){
	if(l == -1000 & u == 1000){
		return("_B")
	}
	if(l == 0 & u == 1000){
		return("_LR")
	}
	if(l == -1000 & u == 0){
		return("_RL")
	}
	stop("not supported direction")
}


cat("::: reactions with different bounds have to get different ids.\n")
for(t in reactFromModelByReact[!checkedReact]){
	for(i in 1:nrow(t)){
		m <- t[i, "model"]
		id <- t[i, "id"]
		lb <- t[i, "lb"]
		ub <- t[i, "ub"]
		
		mod <- allModels[[m]]
		
		react_id(mod)[react_id(mod) == id] <- paste0(id, dirSuffix(lb, ub))
		
		allModels[[m]] <- mod
	}
}









options(warn=1)


cat("::: identifying objective reactions\n")

obj <- sapply(allModels, function(m){
	byName <- grep("biomass", react_id(m), ignore.case=T)
	
	sbin <- S(m) != 0
	byConnection <- which.max(colSums(sbin))
	
	unique(c(byName, byConnection))
})


for(i in names(obj)){
	obj_coef(allModels[[i]])[obj[[i]]] <- 1
}
#objCoef=objCoef[!objCoef[["iAF987"]][1]]

# collect information about the objective functions
cat("::: collect information about objective reactions\n")
objCoef <- lapply(allModels, function(m){
	which(obj_coef(m)!=0)
})

obj <- lapply(allModels, function(m){
	o <- which(obj_coef(m)!=0)
	names(o) <- react_id(m)[o]
	obj_coef(m) <- rep(0, react_num(m))
	
	sba <- sysBiolAlg(m, algorithm="fba")
	lapply(o ,function(x) optimizeProb(sba, react=x, obj_coef=1)) # optimize Each biomass reaction once
})

g <- sapply(lapply(lapply(obj, function(x) sapply(x, "[[", "obj")), unlist), function(y) if(length(y)==0) return(-1) else return(max(y)))
modelSummary <- data.frame(modelsTable, objLength=sapply(objCoef, length), objValue=round(g, digits=6), exchangeReactions=sapply(ex, length))
rownames(modelSummary) <- modelsTable$bigg_id


# recalc exchange reactions
cat("::: recalc exchange reactions\n")
ex <- mclapply(allModels, findExchReact)


cat("::: check metabolite sum formulae\n")



formulaDF <- do.call(rbind, lapply(names(allModels), function(n){
	f <- met_attr(allModels[[n]])$chemicalFormula
	if(is.null(f)){
		f <- ""
	}
	data.frame(met_id=met_id(allModels[[n]]), formula=f, model=n)
}))

formulaDFByMet <- split(formulaDF, formulaDF$met_id)

formulaAlternatives <- lapply(formulaDFByMet, function(x){
	return(setdiff(unique(x$formula), ""))
})
formulaDFUniques <- (formulaDFByMet[sapply(formulaAlternatives, length)<=1])
formulaDFByMet <- formulaDFByMet[sapply(formulaAlternatives, length)>1]

thermo <- get("thermo", CHNOSZ)
data(thermo)
thermo()
options(warn=2)

formulaDFByMet <- lapply(formulaDFByMet, function(df){
	df$formula <- factor(df$formula)
#	browser(expr="ppi[n]"==df$met_id[1])
	
	parsed <- lapply(levels(df$formula), function(x) tryCatch(makeup(x), error=function(e) NULL))
	elements <- unique(unlist(lapply(parsed, names)))
	parsed <- lapply(parsed, function(x){
		v <- rep(0, length(elements))
		if(is.null(x)){
			return(v)
		}
		names(v) <- elements
		v[names(x)] <- x
		return(v)
	})
	parsed <- do.call(rbind, parsed)
	df$version <- as.integer(df$formula)
	for(i in 1:nrow(parsed)){
		for(j in 1:nrow(parsed)){
			if(i < j){
				if(all(parsed[i,] == parsed[j,])){
					#version i and j are the same
					if(!all(parsed[i,]==0)){
						#change to version i#
						print(levels(df$formula)[c(i, j)])
						df$version[df$version == j] <- i
					}
				}
			}
		}
	}
	df
})
formulaDFUniques <- append(formulaDFUniques, formulaDFByMet[sapply(formulaDFByMet, function(x) length(unique(x$version))==1)])
formulaDFByMet <- formulaDFByMet[sapply(formulaDFByMet, function(x) length(unique(x$version))>1)]

options(warn=1)


# make empty strings each as a version itself.
formulaDFByMet <- lapply(formulaDFByMet, function(df){
	if(length(unique(df$version))==2 && "" %in% df$formula){
		print(as.character(df$met_id[1]))
		return(df)
	}
	versionOfEmptyFormula <- df[df$formula=="", "version"]
	versionOfEmptyFormula <- seq(from=max(df$version+1), length.out=length(versionOfEmptyFormula))
	df[df$formula=="", "version"] <- versionOfEmptyFormula
	return(df)
})


cat(":::rename metabolites with double meaning.\n")
metaboliteFormula <- character(0)

for(df in formulaDFByMet){
	if(length(unique(df$version))==2 && "" %in% df$formula){
		cat("skipping, only alternative is empty string.\n")
		print(as.character(df$met_id[1]))
		formula <- setdiff(unique(df$formula), "")[1] # maybe there is an equal alternative.
		stopifnot(length(formula)==1)
		metaboliteFormula[as.character(df$met_id[1])] <- formula
		next
	}
	t <- table(df$version)
	met <- as.character(df$met_id[1])
	cat("renaming metabolite:\n")
	print(met)
	allV <- as.integer(names(t[order(t, decreasing=T)]))
	mainV <- allV[1]
	
	for(v in allV){
		# make new met name
		if(v == mainV){
			newMetName <- met
		}else{
			newMetName <- gsub("^(.+)(\\[.\\])$", paste0("\\1_V", v, "\\2"), met)
		}
		formula <- unique(as.character(df[df$version == v,]$formula))
		#stopifnot(length(formula)==1)
		metaboliteFormula[newMetName] <- formula
		if(v == mainV){
			next
		}
		for(mn in df$model[df$version==v]){
			m <- allModels[[mn]]
			#rename metabolite with version
			metPos <- which(met_id(m)==met)
			met_id(m)[metPos] <-  newMetName
			#find affected reactions
			reactPos <- which(S(m)[metPos,] != 0)
			#as long as reaction id doesnt has the model name in, append the name.
			react_id(m)[reactPos] <- ifelse(grepl(mn, react_id(m)[reactPos], fixed=T), react_id(m)[reactPos], paste0(react_id(m)[reactPos], "_", mn))
			allModels[[mn]] <- m
		}
	}
}

formulae <- sapply(formulaDFUniques, function(x) (as.character(x$formula[1])))


stopifnot(sum(names(metaboliteFormula) %in% names(formulae))==0)
metaboliteFormula[names(formulae)] <- formulae

save(formulaDFByMet, metaboliteFormula, file="~/DATA/formulaDFByMet.Rdata")
##For ecoli supermodel
#save(formulaDFByMet, metaboliteFormula, file="~/DATA/E.coli/formulaDFByMet.Rdata")

#
# Final steps
#
#


cat("::: summarize and save everything\n")
ms <- modelSummary[modelSummary$objValue> 1e-6,]



biomassTable <- sapply(modelSummary$bigg_id[modelSummary$objValue> 1e-6], function(x){
	lapply(objCoef[[x]], function(y)shrinkMatrix(allModels[[x]], j=y))
})
biomassTable <- unlist(biomassTable)

mets <- unique(unlist(sapply(biomassTable, rownames)))
bio <- (unlist(sapply(biomassTable, colnames)))
bio <- paste0(names(bio), "_", bio)
m <- Matrix(0, nrow=length(mets), ncol=length(bio),
dimnames= list(rownames=mets, colnames=bio))

for(i in seq(along=biomassTable)){
	m[rownames(biomassTable[[i]]), i] <- biomassTable[[i]][,1]
}

m <- m[order(rowSums(m!=0), decreasing=T),]
biomassTable <- as.data.frame(as.matrix(m))
mets <- stripCompartment(rownames(m))


write.csv(data.frame(mets, metProp[mets, "name"], biomassTable), file="~/DATA/biomassTable.csv")
#write.csv(data.frame(mets, metProp[mets, "name"], biomassTable), file="~/DATA/E.coli/biomassTable.csv")
save(allModels, file="~/DATA/checkedModels.Rdata")
#save(allModels, file="~/DATA/E.coli/checkedModels.Rdata")




























