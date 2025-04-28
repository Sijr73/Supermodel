#!/usr/bin/Rscript

library(methods)
library(parallel)
library(dplyr)
library(ggplot2)
library(cplexAPI)

#print(load("~/DATA/E.coli/universalBiGG.ver2.Rdata"))
print(load("~/DATA/universalBiGG.ver2.Rdata"))
bioThreshold <- 1e-2


#files <- dir("~/DATA/E.coli/FBAdata/", pattern="*.Rdata", full.names=T)
files <- dir("~/DATA/FBAdata/", pattern="*.Rdata", full.names=T)
df <- do.call(rbind, mclapply(files, function(f){
	get(load(f))
}))

df <- df %>% mutate(growth = ifelse(status==CPX_STAT_INFEASIBLE, 0, growth))
stopifnot(all(df$status %in% c(CPX_STAT_INFEASIBLE, CPX_STAT_OPTIMAL_INFEAS, CPX_STAT_OPTIMAL)))


compareDf <- df %>% select(-status) %>% tidyr::spread(run, growth)
compareDf <- compareDf %>% mutate(fullG = full > bioThreshold, subG = sub > bioThreshold)
compareDf$mediumType <- "seed"
compareDf$mediumType[grep("^R\\|", compareDf$medium)] <- "random"
compareDf$mediumType[grep("^[CNPS]_", compareDf$medium)] <- "minimal"

growthCounts <- compareDf %>% group_by(model, bof, mediumType) %>%
					summarize(
						GG=sum(subG & fullG)/n(),
						NG=sum(!subG & fullG)/n(),
						GN=sum(subG & !fullG)/n(),
						NN=sum(!subG & !fullG)/n()
						)






g <- (ggplot(reshape2::melt(growthCounts, id.vars=c("model", "bof", "mediumType"), variable.name="type"), aes(x=model, y=value, fill=type))+ 
		geom_bar(stat="identity", alpha=1) +
		coord_flip() +
		scale_fill_manual(name = "subset growth / full growth", 
			labels = c(GG="growth/growth", NG="no growth/growth", GN="growth/ no growth", NN="no growth/ no growth"),
			values = c("steelblue", "springgreen", "tomato", "slategrey")
			) +
		facet_grid(bof ~ mediumType) +
		xlab("") +
		ylab("percent") +
		theme(#text = element_text(size = 20),
			legend.position=c(1, 0),
			legend.justification = c(1, 0),
			legend.background = element_rect(fill="#ffffff77"),
			axis.text.y=element_text(angle=45, hjust=1, size=5))
)

if(interactive()){
	print(g)
}else{
	ggsave(plot=g, filename="fbaMediumGrowth.pdf", width = 297, height = 210, units = "mm")
}



if(!all(compareDf %>% filter(medium=="NULLmedium", bof=="energy") %>% select(full, sub) ==0)){
	stop("some model can produce energy without nutrition")
}

if(!all(compareDf %>% filter(medium=="FULLmedium", bof=="energy") %>% select(full, sub) >=10)){
	warning("some models cannot produce energy in FULL medium.")
	
	cat(">>> Following models cannot produce energy in FULL medium:\n")
	print(compareDf %>% filter(medium=="FULLmedium", bof=="energy") %>% select(model, full, sub) %>% filter(sub<=10 | full <= 10))
}


compareDf <- compareDf %>% mutate(NN=!fullG & !subG, NG=!subG & fullG, GN=subG & !fullG, GG=subG & fullG)

nnMedia <- as.character((compareDf %>% filter(mediumType=="minimal", NN, bof=="general") %>% select(medium) %>% unique())[,1])

print(load("~/DATA/metsTrans.Rdata"))
gsub("[CNPS]_", "", nnMedia)



nnMets <- metsTrans[gsub("[CNPS]_", "", nnMedia)]
cat("Percent environments with NN growth and important met is cpd\n")
print(sum(sapply(nnMets, function(x) all(grepl("^cpd\\d+$", x)))) / length(nnMets)*100)


stopifnot(all(names(modelReactMap) %in% df$model))

save(compareDf, file="~/DATA/compareDf.Rdata")
#save(compareDf, file="~/DATA/E.coli/compareDf.Rdata")
#save(compareDf, file="~/DATA/Rebuilding_model/compareDf1.Rdata")








