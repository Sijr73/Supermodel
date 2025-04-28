#!/usr/bin/Rscript
library(methods)

# set seed to some random number to be consisten over multiple script runs.
set.seed(-1208076)

lengthRandomMedia <- 5000

print(load("~/convertMedia/mediaBiGG.Rdata"))

bMedia <- mediaBiGG[grep("^[CNPS]_", names(mediaBiGG))]


print(load("~/metsTrans.Rdata"))
#print(load("~/DATA/metsTrans.Rdata"))
#print(load("New analysis bacillus/DATA/Bacillus/byExchange.Rdata"))
metsTrans=unlist(metsTrans)
byExchange <- split(names(bMedia), gsub("_cpd.+$", "", names(bMedia)))
byExchange <- lapply(byExchange, function(x){
	metsTrans[gsub("^[CNPS]_", "", x)]
})

byExchange <- lapply(byExchange, function(x){
	x[!sapply(x, function(y) any(grepl(y, pattern="^cpd")))]
})
byExchange=lapply(byExchange, function(x) x[!is.na(x)])

save(byExchange, file="~/DATA/byExchange.Rdata")
#print(load("~/DATA/byExchange.Rdata"))

essentials <- setdiff(unique(unlist(bMedia)), unique(unlist(byExchange)))
essentials <- essentials[!grepl("^cpd", essentials)]
#essentials= essentials[1:13]
#essential1 <- scan(text="fe3 ca2 co2 h2o h k mg2 na1 o2",what="character", sep=" ")

allComb <- do.call(expand.grid, lapply(byExchange, function(x) seq(along=x)))
print(dim(allComb))

set.seed(length(bMedia))
randomSelection <- sample.int(nrow(allComb), lengthRandomMedia)

selection <- allComb[randomSelection, ]
selection <- apply(selection, 1, function(x){
	v <- seq(along=byExchange)
	names(v) <- names(byExchange)
	
	lapply(v, function(y) byExchange[[c(y, x[y])]])
})

names(selection) <- sapply(selection, function(x){
	paste(names(byExchange), sapply(x, paste, collapse="/"), sep="|", collapse="|")
})

selectionWithEssentials <- lapply(selection, function(x){
	c(unlist(x, use.names=F), essentials)
})

randomMedia <- selectionWithEssentials
names(randomMedia) <- paste0("R|", names(selectionWithEssentials))
mediaBiGG <- c(mediaBiGG, randomMedia)

save(mediaBiGG, file="~/DATA/mediaBiGGwithRandom.Rdata")


mediaComposition <- data.frame(mediumName=names(randomMedia), mediumType="random")
cnps <- c("C", "N", "P", "S")
names(cnps) <- c("C", "N", "P", "S")

for(comp in cnps){
	mediaComposition[paste0("source", comp)] <- sapply(selection, function(x) paste(x[[comp]], collapse="/"))
}



occurences <- lapply(cnps, function(comp) lapply(bMedia, function(x) intersect(byExchange[[comp]], x)))
mains <- lapply(occurences, function(x){
	t <- table(unlist(x))
	names(t)[which.max(t)]
})

minMediaComposition <- data.frame(mediumName=names(bMedia), mediumType="minimal")
for(comp in cnps){
	minMediaComposition[paste0("source", comp)] <- sapply(occurences[[comp]], function(x){
		s <- setdiff(x, unlist(mains))
		return(ifelse(length(s)==0, NA, s))
	})
}

mediaComposition <- rbind(mediaComposition, minMediaComposition)






for(comp in cnps){
	mediaComposition[[paste0("sourceCount", comp)]] <- sapply(as.character(mediaComposition$mediumName), function(x) sum(unlist(byExchange[[comp]]) %in% mediaBiGG[[x]]))
}



save(mediaComposition, file="~/DATA/mediaComposition.Rdata")












media1 <- sapply(mediaBiGG, function(x){
  mediaBiGG=gsub("...$", "", x)
})








