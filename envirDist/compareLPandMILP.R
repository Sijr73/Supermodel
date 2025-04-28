#!/usr/bin/Rscript

library(methods)
library(ggplot2)
library(dplyr)

print(load("mergedResults.Rdata"))
print(load("envirDistMILP.iAF1260.Rdata"))

writePDF <- T
saveAsPDF <- function(filename=null, plot=g){
	stopifnot(!is.null(filename))
	if(!interactive() || writePDF){
		ggsave(plot, filename=filename, path="plots", width=297, height=210, units="mm")
	}else{
		#print(plot)
	}
}

resultsMILP <- lapply(names(df), function(e){
	lapply(names(df[[e]]), function(m){
		x <- df[[c(e, m)]]
		return(data.frame(model="iAF1260", bof=e, medium=m, addMILP=length(x$solution), relGap=x$relGap, 
			absGap=-x$absGap, # sign error.
			x$solStat))
	})
})
resultsMILP <- do.call(rbind, unlist(resultsMILP, recursive=F))


compareResults <- merge.data.frame(resultsMILP, resultWithModel, by=c("model", "bof", "medium"), all.x=T) %>%
					mutate(exact=absGap == 0, addDiff = add - addMILP, addDiffRel = addDiff / ((add + addMILP) / 2))

save(compareResults, file="compareResults.Rdata")


# compare of LP and MILP solutions
g <- ggplot(compareResults, aes(x=add, y=addMILP, colour=relGap)) +
#	geom_density2d() +
#	stat_density2d(geom="tile", aes(fill=..level..), contour=F) +
	geom_point() + 
	guides(colour=guide_legend(title.position = "top")) + 
#			shape=guide_legend(title.position = "top", nrow=1)) +
#	scale_y_log10() +
	facet_grid(. ~ exact) +
	geom_abline(colour="blue") +
#	theme(legend.position = "none",
#		axis.text.x = element_text(angle = 90, hjust = 1, size=7)) +
	ggtitle("Comparison of LP and MILP solutions for normal BOF and seed medium type.")


saveAsPDF("compareLPandMILP.pdf", g)










