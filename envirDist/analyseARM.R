#!/usr/bin/Rscript

library(methods)
library(sybil)
library(parallel)
library(dplyr)
library(ggplot2)
library(cplexAPI)
print(load("~/DATA/compareDf.Rdata"))
#print(load("~/DATA/E.coli/compareDf.Rdata"))
print(load("~/DATA/mediaBiGGwithRandom.Rdata"))
print(load("~/DATA/universalBiGG.ver2.Rdata"))
#print(load("~/DATA/E.coli/universalBiGG.ver2.Rdata"))
media <- names(mediaBiGG)


files <- dir("~/DATA/envirdata//", pattern="*.Rdata", full.names=T)
#files <- dir("~/DATA/E.coli/eni/Normal//", pattern="*.Rdata", full.names=T)
data <- mclapply(files, function(f){
	get(load(f))
})
names(data) <- gsub("~/DATA/envirdata///", "", gsub(".Rdata", "", files, fixed=T), fixed=T)
#names(data) <- gsub("~/DATA/E.coli/eni/Normal///", "", gsub(".Rdata", "", files, fixed=T), fixed=T)
writePDF <- T

saveAsPDF <- function(filename=null, plot=g){
	stopifnot(!is.null(filename))
	if(!interactive() || writePDF){
		ggsave(plot, filename=filename, path="~/DATA", width=297, height=210, units="mm")
	}else{
		#print(plot)
	}
}

addedReactions <- lapply(names(data), function(mn){
	e <- data[[mn]]
	erg <- lapply(names(e), function(etype){
		lapply(e[[etype]], function(x){
			react_id(uni2)[x$solution]
		})
	})
	names(erg) <- names(e)
	return(erg)
})
names(addedReactions) <- names(data)

df <- lapply(names(data), function(mn){
	print(mn)
	e <- data[[mn]]
	do.call(rbind, lapply(names(e), function(etype){
		d <- e[[etype]]
		if(is.null(d) || length(d) == 0){
			return(data.frame())
		}
		data.frame(
			model=mn,
			medium=names(d),
			bof=etype,
			add = sapply(d, function(x) length(x$solution)),
			solStat = sapply(d, function(x) paste(unique(x$solStat), collapse="-"))
		)
	}))
})
df <- do.call(rbind, df)
rownames(df) <- NULL


growthSummarized <- compareDf %>% 
			select(model, bof, mediumType, NN, NG, GN, GG) %>% 
			group_by(model, bof, mediumType) %>% 
			summarize(NN=sum(NN), NG=sum(NG), GN=sum(GN), GG=sum(GG), count=n())



df %>% filter(solStat != 1) %>% group_by(model) %>% summarize(n()) %>% print()
result <- merge.data.frame(compareDf, df, all.x=T, by=c("model", "medium", "bof"))
result <- result %>% group_by(model, bof, mediumType) %>% mutate(addCor= add - min(add, na.rm=T))
save(result, file="~/DATA/result.Rdata")
#save(result, file="~/DATA/E.coli/result.Rdata")
stopifnot((result %>% filter(NG, is.na(add)) %>% nrow) == 0)

modelTable <- read.csv("~/DATA/modelTable.csv")

#result %>% filter(model=="e_coli_core", NG, !is.na(add)) %>% head() %>% print()



g <- (ggplot(result, aes(add, fill=model))+ 
		geom_histogram(alpha=1, binwidth=5) +
		facet_grid(mediumType ~ bof) +
		guides(fill=guide_legend(title.position = "top", ncol=9)) +
		theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1))
)
saveAsPDF("overlappingHistograms.pdf", g)


# plot histogram for each model
pdf(file="~/DATA/plots/histPerModel.pdf", width=10, height=7, onefile=T)
maxAdd <- max(result$add, na.rm=T)
for (mn in unique(result$model)){
	gm <- (ggplot(result %>% filter(model == mn), aes(add))+ 
		geom_histogram(alpha=0.7, binwidth=5, fill="blue") +
		facet_grid(mediumType ~ bof) +
		xlim(0, maxAdd) +
		ggtitle(paste0(mn, " (", modelTable[mn, "organism"], ")"))
	)
	print(gm)
}
dev.off()






# calculate statistic values for number of added reactions and merge with model table.
cat(">>> building addSummary table\n")
addSummary <- result %>% group_by(model, bof, mediumType) %>%
	summarize(max=max(add, na.rm=T), 
		min=min(add, na.rm=T),
		avg=mean(add, na.rm=T),
		avgCor=mean(addCor, na.rm=T),
		var=var(add, na.rm=T))
addSummary <- merge.data.frame(addSummary, modelTable, all.x=T, by.x=c("model"), by.y=c("bigg_id"))
addSummary <- merge.data.frame(addSummary, growthSummarized, all.x=T, by=c("model", "bof", "mediumType"))


# average number of added reactions vs variance of #added reactions
g <- ggplot(addSummary, aes(x=avg, y=var, colour=model, shape=isEco55)) +
	geom_point() +
	guides(colour=guide_legend(title.position = "top", ncol=7), 
			shape=guide_legend(title.position = "top", ncol=1)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("avgAddLengthVsVariance")

saveAsPDF("avgAddLengthVsVariance.pdf", g)





# average number of added reactions vs reaction count
g <- ggplot(addSummary, aes(x=reaction_count, y=avg, colour=model, shape=isEco55)) +
	geom_point() +
	guides(colour=guide_legend(title.position = "top", ncol=7), 
			shape=guide_legend(title.position = "top", nrow=2)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("avgAddLengthVsReactionCount")

saveAsPDF("avgAddLengthVsReactionCount.pdf", g)


# average number of added reactions vs non viable environments
g <- ggplot(addSummary, aes(x=NG/count, y=avg, colour=model, shape=isEco55)) +
	geom_point() +
	guides(colour=guide_legend(title.position = "top", ncol=7), 
			shape=guide_legend(title.position = "top", nrow=2)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("avgAddLengthVsSubNonViableEnvir")

saveAsPDF("avgAddLengthVsSubNonViableEnvir.pdf", g)


# average number of added reactions vs viable environments
g <- ggplot(addSummary, aes(x=GG/count, y=avg, colour=model, shape=isEco55)) +
	geom_point() +
	guides(colour=guide_legend(title.position = "top", ncol=7), 
			shape=guide_legend(title.position = "top", nrow=2)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("avgAddLengthVsSubViableEnvir")

saveAsPDF("avgAddLengthVsSubViableEnvir.pdf", g)



# average number of added reactions vs reaction count for eco 55 only
g <- ggplot(addSummary %>% filter(isEco55), aes(x=reaction_count, y=avg, colour=model)) +
	geom_point() +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("avgAddLengthVsReactionCountEco55")

saveAsPDF("avgAddLengthVsReactionCountEco55.pdf", g)


cat(">>> merging full results with model table \n")
resultWithModel <- merge.data.frame(result, modelTable, all.x=T, by.x=c("model"), by.y=c("bigg_id"))


# single reaction additions vs reaction numbers
g <- ggplot(resultWithModel, aes(x=reaction_count, y=add, colour=model, shape=isEco55)) +
	geom_point() +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("singleAddsVsReactNum")

saveAsPDF("singleAddsVsReactNum.pdf", g)





# violin of added reactions vs number reactions
g <- ggplot(resultWithModel, aes(x=(reaction_count), y=add, colour=model, group=factor(model))) +
	geom_violin(position = "identity", scale="width") +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 90, hjust = 1)) +
	ggtitle("singleAddsVsReactNumViolin")

saveAsPDF("singleAddsVsReactNumViolin.pdf", g)


# violins of added reactions
g <- ggplot(resultWithModel, aes(x=factor(model), y=add, colour=model)) +
	geom_violin() +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	scale_y_log10() +
	facet_grid(mediumType ~ bof) +
	theme(legend.position = "none",
		axis.text.x = element_text(angle = 90, hjust = 1, size=7)) +
	ggtitle("singleAddsViolin")

saveAsPDF("singleAddsViolin.pdf", g)









save(addSummary, resultWithModel, result, modelTable, addedReactions, file="~/DATA/mergedResults.Rdata")
#save(addSummary, resultWithModel, result, modelTable, addedReactions, file="~/DATA/E.coli/mergedResults.Rdata")
#save(addSummary, resultWithModel, result, modelTable, addedReactions, file="~/DATA/Rebuilding_model/mergedResults.Rdata")
















