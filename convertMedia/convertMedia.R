#!/usr/bin/Rscript
library(methods)
library(sybil)
library(rjson)
library(parallel)
library(dplyr)
library(CHNOSZ)

options(warn=2)

print(load("../sourceData/seed.mobd.media.Rdata"))
print(load("../sourceData/balazs.media.Rdata"))

print(load("~/DATA/formulaDFByMet.Rdata"))

mets <- unique(c(unlist(balazs.media), unlist(seed.mobd.media)))
mets=met_id(uni1.1)
names(mets) <- mets
mets=na.omit(mets)
media <- c(balazs.media, seed.mobd.media)
########################################################

print(load("~/models/universalBiGG.ver1.1.Rdata"))
ex <- findExchReact(uni1.1)
exMetsInModel <- gsub("\\[.+$", "", met_id(ex))
isMetWithVersion <- grep("_V\\d$", exMetsInModel)
exMetsInModel <- gsub("_V\\d$", "", exMetsInModel)
exMetsComparts <- gsub("^.+(\\[.\\])$", "\\1", met_id(ex))

##http://bigg.ucsd.edu/api/v2/universal/metabolites/g3p
metInfo <- read.table("~/sourceData/bigg_models_metabolites.txt", sep="\t", comment.char="", quote="", header=T, stringsAsFactors=F)
metInfo$metIds <- gsub("_.$", "", metInfo$bigg_id)
metInfo <- metInfo[!duplicated(metInfo$metIds),]
rownames(metInfo) <- metInfo$metIds

##metInfo <- metInfo[metsInModel,]
##stopifnot(all(!is.na(metInfo$metIds)))

xref <- lapply(metInfo$database_links, function(x) fromJSON(json_str=x))
#names(xref) <- metInfo$metIds

#seedIds <- lapply(xref, function(x){
#	if(!is.null(x[["SEED Compound"]])){
#		return(sapply(x[["SEED Compound"]], "[[", c("id")))
#	}else{
#		return(NULL)
#	}
#})

#seedMap <- do.call("rbind", lapply(names(seedIds[!sapply(seedIds, is.null)]), function(x) data.frame(bigg=x, seed=seedIds[[x]])))
#seedMapDict <- as.character(seedMap$seed)
#names(seedMapDict) <- seedMap$bigg

#seedMapDict["cpd00027"] <- "glc__D"


##XREF   MNX_ID  Evidence        Description
#t <- read.table("chem_xref.tsv", sep="\t", comment.char="", quote="", skip=126, stringsAsFactors=F)
#colnames(t) <- c("XREF", "MNX_ID", "Evidence", "Description")
#t <- t[grep("^seed:|^bigg:", t$XREF),]
#t$DB <- gsub("^(.+):.+$", "\\1", t$XREF)
#t$XREF <- gsub("^.+:", "", t$XREF)
#t$Description <- ""

#ts <- split(t, (t$MNX_ID))

#seedBiggMapping <- lapply(ts, function(x) expand.grid(seed=x$XREF[x$DB == "seed"], bigg=x$XREF[x$DB == "bigg"]))
#seedBiggMapping <- do.call(rbind, seedBiggMapping)
#seedBiggMapping <- split(seedBiggMapping, seedBiggMapping$seed)

#seedTable <- read.table("ModelSEED-compounds-db.tsv", sep="\t", quote="", comment.char="", header=T, stringsAsFactors=F)
#seedMapDict <- seedTable$ABBREVIATION
#names(seedMapDict) <- seedTable$DATABASE
#seedMapDict <- na.omit(seedMapDict[mets])

seedTable <- read.table("New analysis bacillus/networkComplexityBigg-master/convertMedia/SEEDCompounds.tsv", sep="\t", quote="", comment.char="", header=T, stringsAsFactors=F)
rownames(seedTable) <- seedTable$id
seedTable <- seedTable[mets, ]
seedTable <- seedTable[!is.na(seedTable$id),]

seedAliases <- lapply(seedTable[,"aliases"], function(x){
	text <- scan(text=x, sep=";", quote="\"", what="character", quiet=T)
	m <- regexec("^([^:]+):(.+)$", text)
	table <- do.call(rbind, regmatches(text, m))[,2:3]
	res <- unique(table[table[,1] %in% names(modelReactMap), 2])
	res <- gsub("(\\w)-([DLRMS])$", "\\1__\\2", res)
	#browser(condition=print(res), expr=any(grepl("^ala", res)))
	gsub("-", "_", res)
})
names(seedAliases) <- seedTable$id

seedAliases <-seedTable$id

#seedBiggMapping <- seedBiggMapping[mets]
#seedBiggMapping <- seedBiggMapping[!sapply(seedBiggMapping, is.null)]
#seedBiggMapping <- sapply(seedBiggMapping, function(x){
#	if(nrow(x)> 1){
#		x <- x[x$bigg %in% exMetsInModel, ]
#		print(x)
#	}
#	return(x)
#})

metsTrans <- lapply(mets, function(x){
	#new <- as.character(na.omit(union((seedBiggMapping[[x]]$bigg), seedMapDict[x])))
	all <- na.omit(as.character(seedAliases[[x]]))
	
	new <- all[all %in% exMetsInModel]
	
	if(length(new)==0){
		if(length(all) > 0){
			print(all)
		}
		return(x)
	}
	unique(new)
})


versionMetabolites <- exMetsInModel[isMetWithVersion]
mappingToCheck <- metsTrans[sapply(metsTrans, function(x) any(x %in% versionMetabolites))]
stopifnot(all(sapply(mappingToCheck, length)==1))
data(thermo)
validVersionMappings <- lapply(names(mappingToCheck), function(n){
	seedid <- n
	biggid <- mappingToCheck[[seedid]]
	
	df <- formulaDFByMet[[paste0(biggid, "[e]")]]
	
	stopifnot(is.data.frame(df))
	
	parseSumFormula <- function(x) tryCatch(makeup(x), error=function(e) NULL)
	parsed <- lapply(levels(df$formula), parseSumFormula)
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
	
	seedFormula <- parseSumFormula(seedTable[seedid,"formula"])
	equalFormula <- sapply(1:nrow(parsed), function(i){
		comp <- parsed[i, names(seedFormula)] == seedFormula
		all(comp[names(comp)!="H"])
	})
	# return versions that are compatible with seed formula
	numberV <- unique(df$version[equalFormula[as.integer(df$formula)]])
	
	mainV <- as.integer(names(sort(decreasing=T, table(df$version)))[1])
	numberV[numberV==mainV] <- ""
	if(length(numberV)> 0){
		ids <- paste0(biggid, ifelse(numberV=="", "", "_V"), numberV, "")
		return(ids)
	}else{
		return(character(0))
	}
})
names(validVersionMappings) <- names(mappingToCheck)

metsTrans[names(validVersionMappings)] <- validVersionMappings

#metsTrans <- as.character(seedMap$bigg[match(mets, seedMap$seed)])
#names(metsTrans) <- mets

save(metsTrans, file="metsTrans.Rdata")

#stopifnot(setequal(mets[(is.na(metsTrans))], c("cpd81911", "cpd81814")))
# these two dont exist any more in the seed
#> mets[(is.na(metsTrans))]
#[1] "cpd81911" "cpd81814"


mediaBiGG <- lapply(media, function(x) unlist(metsTrans[x]))
#mediaBiGG <- lapply(mediaBiGG, function(x) x[!is.na(x)])
stopifnot(all(!is.na(unlist(mediaBiGG))))

# add essential metabolites to the media
mediaBiGG <- lapply(mediaBiGG, function(x) union(x, "ni2"))

#[1] "MNXM128"   "MNXM43"    "MNXM95"    "MNXM2255"  "MNXM1026"  "MNXM89361"
#> chemProb[a, ]
#                 V1        V2   V3 V4      V5                      V6                   V7          V8
#MNXM128     MNXM128    Ca(2+)   Ca  2  40.078         InChI=1S/Ca/q+2               [Ca++] chebi:29108
#MNXM43       MNXM43     Cl(-)   Cl -1  35.453    InChI=1S/ClH/h1H/p-1                [Cl-] chebi:17996
#MNXM95       MNXM95      K(+)    K  1  39.098          InChI=1S/K/q+1                 [K+] chebi:29103
#MNXM2255   MNXM2255    Mn(2+)   Mn  2  54.938         InChI=1S/Mn/q+2               [Mn++] chebi:29035
#MNXM1026   MNXM1026 molybdate MoO4 -2 159.938 InChI=1S/Mo.4O/q;;;2*-1 [O-][Mo]([O-])(=O)=O chebi:36264
#MNXM89361 MNXM89361 zinc atom   Zn  0  65.409             InChI=1S/Zn                 [Zn] chebi:27363
#essential <- c("MNXM128", "MNXM43", "MNXM95", "MNXM2255", "MNXM1026", "MNXM89361")

# TODO change this later?!?!
#essential <- union(essential, c("MNXM51417", "MNXM78334", "MNXM81014"))

#mediaMNX <- lapply(mediaMNX, function(x) return(unique(c(x, essential))))


mediaBiGG$FULLmedium <- exMetsInModel
mediaBiGG$NULLmedium <- character(0)

#m <- unique(unlist(mediaMNX))
#mediaMNX <- lapply(m, function(x) return(setdiff(m, x)))
#names(mediaMNX) <- m

biggMets <- unique(unlist(mediaBiGG))
print(biggMets[!biggMets %in% exMetsInModel])


save(mediaBiGG, file="mediaBiGG.Rdata")





