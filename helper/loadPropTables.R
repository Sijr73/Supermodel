
metProp <- read.table("~/sourceData/bigg_models_metabolites.txt", sep="\t", quote="", header=T, stringsAsFactors=F,fill = TRUE)
rownames(metProp) <- metProp$bigg_id
reacProp <- read.table("~/sourceData/bigg_models_reactions.txt", sep="\t", quote="", header=T, stringsAsFactors=F,fill=TRUE)







