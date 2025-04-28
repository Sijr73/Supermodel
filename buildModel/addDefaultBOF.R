generalBOF <- read.table("generalBOF.tsv", header=T, sep="\t")
energyBOF <- read.table("energyBOF.tsv", header=T, sep="\t")

bof <- energyBOF
bof <- bof[bof$bigg_id!="", ]
bof <- sapply(split(bof$Scoef, bof$bigg_id), sum)
mod <- addReact(uni, id="ENERGY_BOF", met=names(bof), Scoef=bof)

bof <- generalBOF
bof <- bof[bof$bigg_id!="", ]
bof <- sapply(split(bof$Scoef, bof$bigg_id), sum)
uni <- addReact(mod, id="GENERAL_BOF", met=names(bof), Scoef=bof)

newBOF <- c("ENERGY_BOF", "GENERAL_BOF")
modelReactMap <- lapply(modelReactMap, function(x) union(x, newBOF))
