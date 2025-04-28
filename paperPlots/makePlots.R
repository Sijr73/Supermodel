#!/usr/bin/Rscript

library(methods)
library(parallel)
library(sybil)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(gtable)
library(ggrepel)
library(nls2)
library(xlsx)
library(xtable)




modelSelector <- scan(text="iAM_Pb448 iCN718 iCN900 iNF517 iLB1027_lipid iJB785 iIS312 iEK1008 iAF692 iIT341 iLJ478 iSM199 iSB619 iHN637 iJN678 iJN746 iYO844 iAF987 iMM904 iPC815 iRC1080 iYL1228 STM_v1_0 iML1515",
					what="character", sep=" ")


myBlue <- "#006ab3"
myBlue2 <- "#008cc2"
myGreen <- "#97bf0d"
myGrey <- "#D9DADB"
myRed <- "#BE0A26"
myYellow <- "#FFCC00"

file.remove(list.files("~/plots/", full.names=T), showWarnings=F)
file.remove("~/plot/plots.pdf", showWarnings=F)


if(!exists("dataLoaded")){
	print(load("~/DATA/Rebuilding_model/universalBiGG.ver21.Rdata"))
	print(load("~/DATA/Rebuilding_model/mergedResults1.Rdata"))
	print(load("~/DATA/Rebuilding_model/compareDf1.Rdata"))
	print(load("~/DATA/innovationIndex1.Rdata"))
	#print(load("../envirDist/compareResults.Rdata"))
	print(load("~/DATA/innovationIndexClusters1.Rdata"))
	print(load("~/DATA/exaptation1.Rdata"))
	print(load("~/DATA/metaboliteParticipation1.Rdata"))
	
	growthCounts <- compareDf %>% group_by(model, bof, mediumType) %>% summarize(NN=sum(NN)/n(), NG=sum(NG)/n(), GN=sum(GN)/n(), GG=sum(GG)/n())
	
	
	
	modelTable <- read.csv("~/DATA/modelTable1.csv") %>% filter(bigg_id %in% names(modelReactMap))
	modelTable <- modelTable %>% left_join(read.csv("~/paperPlots/shortOrganismNames1.csv") %>% select(bigg_id, shortO, sO))
	
	modelTable <- modelTable %>% mutate(label= paste0(organism, " (", bigg_id, ")"))
	modelTable <- modelTable %>% mutate(labelSO= paste0(sO, " (", bigg_id, ")"))
	modelTable <- modelTable %>% mutate(labelSOE= paste0("italic(", gsub(" ", "~", sO), ")", "~(", bigg_id, ")"))
	modelTable <- modelTable %>% mutate(ecoLabel= factor(ifelse(isEco55, "55~italic(E.Coli)", "other~species"), levels=c("other~species", "55~italic(E.Coli)")))
	modelTable <- modelTable %>% mutate(#shortO = sapply(strsplit(as.character(organism), "\\s"), function(x) paste(x[1:2], collapse=" ")),
										shortLabel = paste0(shortO, " (", bigg_id, ")"))
	
	auxotrophic <- c("ic_1306", "iECDH10B_1368", "iECIAI39_1322", "iECUMN_1333", 
		"iSbBS512_1146", "iSBO_1134", "iSFV_1184", "iSFxv_1172", "iS_1188", 
		"iSF_1195", "iSSON_1240", "iSDY_1059",
		"iSM199", "iPC815", "iIT341", "iAF692", "iNJ661", "iSB619","iEcDH1_1363") # manually added for non eco55
	modelTable <- modelTable %>% mutate(auxo = bigg_id %in% auxotrophic)
	modelTable <- modelTable %>% left_join(metaboliteParticipation %>% filter(y==2) %>% select(name, x) %>% rename(bigg_id=name, twoReactMetRatio = x))
	
	addSummary <- addSummary %>% mutate(auxo = model %in% auxotrophic)
	
	
	resultWithModel <- resultWithModel %>% left_join(modelTable %>% select(bigg_id, shortLabel, labelSO, labelSOE, ecoLabel), by=c(model="bigg_id"))
	growthCounts <- growthCounts %>% left_join(modelTable %>% select(bigg_id, shortLabel, labelSO, labelSOE), by=c(model="bigg_id"))
	
	
	
	df <- result %>% filter(bof=="general", mediumType=="random")
	df2 <- innovationIndex %>% filter(bof=="general", mediumType=="random") %>% select(-NN, -NG, -GG, -GN)
	
	innovationIndexSingle <- merge.data.frame(df, df2, by.x=c("model", "bof", "mediumType", "medium"), by.y=c("model", "bof", "mediumType", "addedForMedium"), all.x=T)
#	innovationIndexSingle <- innovationIndexSingle %>% select(-c(organism.y, X.y, label.y, ))
	innovationIndexSingle <- merge.data.frame(innovationIndexSingle, modelTable, by.x="model", by.y="bigg_id", all.x)
	
	innovationIndex <- merge.data.frame(innovationIndex, modelTable, by.x="model", by.y="bigg_id", all.x)
	df <- innovationIndex %>% group_by(model, bof, mediumType) %>% summarize(innoIndexRelMean=mean(innoIndexRel), innoIndexDistinctRelMean=mean(innoIndexDistinctRel))
	addSummary <- merge.data.frame(addSummary, df, by=c("model", "bof", "mediumType"), all.x=T)
	addSummary <- addSummary %>% left_join(modelTable %>% select(bigg_id, labelSO, labelSOE, ecoLabel) %>% rename(model=bigg_id))
	
	innovationIndexWithClusters <- innovationIndexWithClusters %>% left_join(modelTable %>% select(bigg_id, gene_count, labelSOE, ecoLabel), by=c(model="bigg_id"))
	
	
	
	innovationIndex <- innovationIndex %>% left_join(exaptation, by=c(model="model", addedForMedium="second"))
	innovationIndex <- innovationIndex %>% left_join(result %>% select(model, medium, bof, mediumType, add, addCor), c(model="model", mediumType="mediumType", addedForMedium="medium", bof="bof"))
	
	
	
	mycolor <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB")
	CCB_plot_style <- function(...)
	{
	  
	  theme_bw()+
	    
	    
	    theme(
	      text = element_text(size=16,face="bold",color="black",family="sans"),
	      axis.title.x = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 5, r = 0, b = 0, l = 0)),
	      axis.title.y.left  = element_text(color = "black", size = 20, face = "bold",margin = margin(t = 0, r = 5, b = 0, l = 0)),
	      axis.text.x=element_text(size=18,family="sans",face="bold",colour="#666666"),
	      axis.text.y=element_text(size=18,family="sans",face="bold",colour="#666666"),
	      plot.margin = margin(10, 10, 10, 10),
	      axis.title.y = element_text(color = "black", size = 16, face = "bold", margin = margin(t = 0, r = 10, b = 0, l = 0)),
	      legend.box.background = element_rect(colour = "black",linewidth = 1),
	      panel.grid.major = element_blank(),
	      panel.grid.minor = element_blank(),
	      panel.border = element_blank() ,
	      panel.background = element_blank(),
	      axis.line = element_line(color = '#666666',linewidth = 0.8),
	      axis.ticks =element_line(color="#666666",linewidth =1),
	      axis.ticks.length=unit(c(-0.3,0.3), "cm"))
	  
	  
	}
	

	
	dataLoaded <- T
}


table1 <- modelTable %>% arrange(gene_count) %>% select(bigg_id, organism, ends_with("_count"), isEco55, taxid, pubMedId, year)
#write.csv(table1, row.names = F, file="tableS1.csv")
rownames(table1) <- NULL
table1 <- table1 %>% select(-year) %>% mutate(isEco55=ifelse(isEco55, "X", "")) %>%
					 rename("Model ID"=bigg_id, 
							"Organism"=organism, 
							"Gene count"=gene_count, 
							"Metabolite count"=metabolite_count, 
							"Reaction count"=reaction_count, 
							"55\\textit{E.coli}"=isEco55,
							"Taxonomy ID"=taxid,
							"PubMed ID"=pubMedId)
write.xlsx(table1, showNA=F, file=paste0("~/Final plot Feb2024 supermodel/tableS1.xlsx"))


addtorow          <- list()
addtorow$pos      <- list()
addtorow$pos[[1]] <- c(0)
addtorow$command  <- c(paste("\\hline \n",
							"\\endhead \n",
							"\\hline \n",
							"{\\footnotesize Continued on next page} \n",
							"\\endfoot \n",
							"\\endlastfoot \n",sep=""))
print(xtable(table1, align=c("l", "p{4cm}", "p{6cm}", "r", "r", "r", "c", "r", "r"),
					 caption="Organism specific models and their properties"
					),
					file="~/plot/tableS1.tex", table.placement="",
											include.rownames=F,
											rotate.colnames=F,
											floating=F,
											sanitize.colnames.function = function(x) {x},
											tabular.environment = "longtable",
											floating.environment="",
											caption.placement = "top",
											add.to.row = addtorow,
											hline.after=c(-1))

myParse <- function(x) parse(text=x)

gS1 <- ggplot(modelTable %>% arrange(reaction_count), aes(x=reorder(labelSOE, order(reaction_count)), y=reaction_count, fill=isEco55)) +
	geom_bar(stat="identity") + 
	labs(x="", y="reaction count") + 
	scale_fill_manual(name="is E.Coli 62", values=c(mycolor[1], mycolor[5]), breaks=c(FALSE, TRUE)) +
	scale_x_discrete(labels=myParse) +
	coord_flip() +
  theme_bw()+
  
  
  theme(
    text = element_text(size=10,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x=element_text(size=10,family="sans",face="bold",colour="#666666"),
    axis.text.y=element_text(size=8.5,family="sans",face="bold",colour="#666666",angle = 0, hjust = 1, vjust=0.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 10, face = "bold", margin = margin(t = 5, r = 10, b = 5, l = 5)),
    legend.box.background = element_rect(colour = "black",linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank() ,
    
    panel.background = element_blank(),
    axis.line = element_line(color = '#666666',linewidth = 0.8),
    axis.ticks =element_line(color="#666666",linewidth =1),
    axis.ticks.length=unit(c(-0.3,0.3), "cm"))+
  theme(plot.title=element_text(hjust=0, size=15)) +
  
  
	facet_grid(ecoLabel ~ ., scales="free_y", space="free_y", labeller=label_parsed) +
	theme(#axis.text.x = element_text(size=10, angle = 0, hjust = 1, vjust=0.5),
			legend.position="none") +
	ggtitle(bold("Networks sizes")~"")

gS1
ggsave(gS1, filename="FigS1.pdf", path="~/plot/", width=210, height=297, units="mm")


myLabels <- modelTable %>% filter(bigg_id %in% modelSelector) %>% arrange(reaction_count) %>% select(labelSOE) %>% .[[1]]

g1 <- ggplot(modelTable %>% filter(bigg_id %in% modelSelector) %>% arrange(reaction_count)%>%
               mutate(color = ifelse(labelSOE == "italic(S.~cerevisiae)~(iMM904)" | labelSOE == "italic(Chlamydomonas)~(iRC1080)"
                                     | labelSOE == "italic(P.~tricornutum)~(iLB1027_lipid)",
                                     mycolor[2], mycolor[1])), aes(x=reorder(labelSOE, order(reaction_count)), y=reaction_count,fill=color)) +
	geom_bar(stat="identity") + 
  scale_fill_manual("legend",values=c("#4477AA"="#4477AA","#EE6677"="#EE6677"))+
	scale_x_discrete(labels=parse(text=myLabels)) +
	coord_flip() +
	labs(x="", y="reaction count") + 
 CCB_plot_style()+

	theme(	axis.text.y = element_text(size=14, angle = 0,family="sans",face="bold",colour="#666666", hjust = 1, vjust=0.5),  

			legend.position="none")#+
#	ggtitle("Number of reactions per model")
g1
ggsave(g1, filename="Fig1.pdf", path="~/plot/")



count <- table(met_comp(uni2))
compCount <- data.frame(id=mod_compart(uni2)[as.integer(names(count))], value=as.integer(count)) %>% 
			arrange(desc(value)) %>% mutate(id=factor(id, levels=id))

gS2a <- ggplot(compCount, aes(x=id, y=value)) +
	geom_bar(stat="identity", fill=mycolor[1]) + 
	labs(x="compartment id", y="metabolite count") + 
#	scale_y_log10() +
	ggtitle(bold("Number of metabolites in compartments")~"") +
  CCB_plot_style()+   

	theme(plot.title=element_text(color = "black", face = "bold",hjust=0, size=15))
gS2a
ggsave(gS2a, filename="FigS2a.pdf", path="~/plot/", width=210, height=297/2, units="mm")

gS2b <- ggplot(modelTable %>% filter(bigg_id %in% modelSelector) %>% group_by(year) %>% summarize(yearCount=n()), aes(x=year, y=yearCount)) +
	geom_bar(stat="identity", fill=mycolor[1]) +
	xlab("year") +
	ylab("count") +
	ggtitle(bold("Published models per year")~"") +
	scale_x_continuous(breaks=seq(from=2005, to=2022, by=2), minor_breaks=seq(from=2005+1, to=2022, by=2), labels=seq(from=2005, to=2022, by=2)) +
  CCB_plot_style()+
  
  theme(plot.title=element_text(hjust=0, size=15))
gS2b
ggsave(gS2b, filename="FigS2d.pdf", path="~/plot/", width=210, height=297/2, units="mm")

gS2 <- grid.arrange(gS2a, gS2b, ncol=1, nrow =2)
ggsave(gS2, filename="FigS2.pdf", path="~/plot/", width=210, height=297, units="mm")







count <- table(unlist(lapply(modelReactMap, function(x) unique(gsub("_((LR)|(RL)|B)", "", x)))))
reactionCount <- data.frame(id=names(count), value=as.integer(count))
a <- reactionCount %>% filter(id %in% gsub("_((LR)|(RL)|B)", "", modelReactMap$e_coli_core))

b <- a %>% filter(value < 20) %>% select(id) %>%  .[,1] %>% as.character()
c <- table(unlist(lapply(b, function(y) names(which(sapply(modelReactMap, function(x) y %in% x))))))

# alternatives for singleton E.coli reactions.
lapply(b, grep ,react_id(uni2), value=T)



g2 <- ggplot(reactionCount, aes(x=value)) +
	geom_histogram(binwidth=2, fill=myBlue) + 
	geom_histogram(data=a, binwidth=2, aes(x=value), fill=myGreen) +
	scale_y_log10() +
	labs(y="frequency (log)", x="# occurrences of reaction in all models") +
	#ggtitle("Reaction frequencies (log scale) over all models") +
	theme(plot.title=element_text(hjust=0, size=9),
			plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
g2
ggsave(g2, filename="Fig2.pdf", path="~/plot/", width=297/2, height=210/2, units="mm")





a %>% arrange(value) %>% 
		mutate(count=value >= 60) %>% 
		group_by(count) %>% 
		summarize(n=n()/nrow(.)) %>% 
		print()


nbof <- c(
	normal="model specific biomass reaction",
	energy="energy production",
	general="general biomass reaction"
)
nmt <- c(
	minimal="minimal environments",
	random="random minimal environments",
	seed="wet lab (the seed) environments"
)

gS3 <- (ggplot(reshape2::melt(growthCounts %>% select(-GN), id.vars=c("model", "bof", "mediumType", "shortLabel", "labelSO", "labelSOE"), variable.name="type"), aes(x=model, y=value, fill=type))+ 
		geom_bar(stat="identity", alpha=1) +
#		coord_flip() +
		scale_fill_manual(name = "subset growth / full growth:", 
			labels = c(GG="growth/growth", NG="no growth/growth", NN="no growth/ no growth"),
			values = c(mycolor[2], mycolor[1], mycolor[3])
			) +
		#facet_grid(bof ~ mediumType) +
		facet_wrap(bof ~ mediumType, nrow=9, ncol=1, labeller = labeller(bof=nbof, mediumType=nmt, .multi_line=FALSE)) +
		xlab("") +
		ylab("percent") +
		ggtitle(bold("Growth for submodels and supermodel")~"") +

		theme(#text = element_text(size = 20),
			#legend.position=c(1, 1),
			#legend.justification = c(1, 1),
			legend.position="bottom",
			legend.background = element_rect(fill="#ffffff77"),
			axis.text.x=element_text(angle=90, hjust=1, size=7.5),
			axis.text.y=element_text(size=8),
			axis.title.y.left =element_text(size=14), 
			plot.title=element_text(hjust=0, size=15),
			plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
)
gS3
ggsave(gS3, filename="FigS3.pdf", path="~/plot/", width=210, height=297, units="mm")



cat("not growing in FULLmedium:\n")
resultWithModel %>% filter(medium == "FULLmedium", bof == "general", !subG) %>% select(model, organism) %>% print()


initialGrowth <- addSummary %>% filter(mediumType=="seed", bof=="general") %>% mutate(initialGrowthRatio = GG/(NN+NG+GG))
#midPoint <- with(initialGrowth, min(initialGrowthRatio) + (max(initialGrowthRatio) - min(initialGrowthRatio)) / 2)
midPoint=0.2
#energy
#midPoint=0.75
specialists <- initialGrowth %>% filter(initialGrowthRatio <= midPoint) %>% select(model) %>% .[[1]] %>% as.character()
addSummary <- addSummary %>% mutate(sgCat=as.factor(ifelse(model %in% specialists, "specialist", "generalist")))
resultWithModel <- resultWithModel %>% mutate(sgCat=as.factor(ifelse(model %in% specialists, "specialist", "generalist")))
innovationIndex <- innovationIndex %>% mutate(sgCat=as.factor(ifelse(model %in% specialists, "specialist", "generalist")))
modelTable <- modelTable %>% mutate(sgCat=as.factor(ifelse(bigg_id %in% specialists, "specialist", "generalist")))




# order growth percentage
#order <- growthCounts %>% filter(model %in% modelSelector, bof=="general", mediumType=="seed") %>% arrange(GG, NN) %>% ungroup() %>% select(labelSO) %>% .[[1]] %>% as.character()

order <- modelTable %>% filter(bigg_id %in% modelSelector) %>% arrange(gene_count) %>% select(labelSOE) %>% .[[1]] %>% as.character()

g3 <- (ggplot(reshape2::melt(
			growthCounts %>% filter(model %in% modelSelector, bof=="general", mediumType %in% c("random", "seed")) %>% 
								select(-GN, -NN, -NG) %>% mutate(GG=GG*ifelse(mediumType == "random", -1, 1)),
								id.vars=c("model", "bof", "mediumType", "shortLabel", "labelSO", "labelSOE"), variable.name="type"),
 								aes(x=factor(labelSOE, order), y=value, group=mediumType, fill=paste0(mediumType, model %in% specialists) ))+ 
		geom_bar(stat="identity", alpha=1) +
		geom_hline(yintercept=midPoint) +
		coord_flip() +
		guides(fill=F) +
		scale_x_discrete(labels=parse(text=order)) +
#		scale_fill_manual(name = "submodel/fullmodel", 
#			labels = c(GG="growth/growth", NG="no growth/growth", NN="no growth/ no growth"),
#			values = c(myBlue, myGreen, myGrey)
#			) +
		scale_fill_manual(name = "medium type", 
			values = c(randomFALSE=mycolor[3], seedFALSE=mycolor[1], seedTRUE=mycolor[2], randomTRUE=mycolor[3])) +
		scale_y_continuous(breaks=c(-0.8,-0.5,-0.2, 0,0.2,0.5, midPoint ), labels=c("0.8","0.5","0.2", "0","0.2","0.5", round(midPoint, digits=2) )) +
		xlab("") +
		ylab("fraction of environments") +
		ggtitle(~bold("Random environments | Wet lab environments")) +
		theme(#text = element_text(size = 20),
			legend.position=c(1, 0),
			legend.justification = c(1, 0),
			legend.background = element_rect(fill="#ffffff77"),
			#plot.title=element_text(hjust=0.445, size=10),
			axis.text.y=element_text(angle=0, hjust=1,lineheight = 2, size=10),
		)+  CCB_plot_style()+
  theme( plot.title=element_text(hjust=0.22, size=15,face="bold"),
         axis.text.y = element_text(size=14, angle = 0,family="sans",face="bold",colour="#666666", hjust = 1, vjust=0.5))

)
g3

ggsave(g3, filename="Fig3.pdf", path="~/plot/", width=220, height=297/2, units="mm")
#ggsave(g3, filename="Fig3energy2.pdf", path="~/plot/", width=220, height=297/2, units="mm")


gS4 <- ggplot(initialGrowth %>% filter(model %in% modelSelector | isEco55), aes(x=initialGrowthRatio, fill=ifelse(initialGrowthRatio <=midPoint, "S", "G"))) +
	geom_histogram(binwidth=0.02) +
	geom_vline(xintercept=midPoint) +
	scale_fill_manual(values=c(G=mycolor[1], S=mycolor[2])) +
	guides(fill=FALSE) +
	xlab("fraction of viable seed environments") +
	ylab("count") +
	#coord_flip() +
  scale_y_continuous(expand=c(0,0))+
	ggtitle(bold("Growth in wet lab media")~"") +
	theme(plot.title=element_text(hjust=0, size=15))+CCB_plot_style()
	
	gS4
ggsave(gS4, filename="FigS4.pdf", path="~/plot/", width=220, height=297/2, units="mm")



#       (b-a)(x - min)
#f(x) = --------------  + a
#          max - min
#adddf <- adddf %>% mutate(x=(7768^2 - (reaction_count^2)    ),
#				a= min(avg/max(avg)),
#				b= max(avg/max(avg)),
#				toolboxScaled= ((b-a) * (x - min(x)) / (max(x) - min(x))) + a
#) %>% select(-x, -a, -b)





# average number of added reactions vs reaction count
adddf <- addSummary %>% filter(bof=="general", mediumType=="seed")
#adddf=adddf %>% filter(model !="iIS312")
#adddf=adddf %>% filter(model !="iIS312_Epimastigote")
#adddf=adddf %>% filter(model !="iIS312_Trypomastigote")
#adddf=adddf %>% filter(model !="iIS312_Amastigote")
#adddf=adddf %>% filter(model !="iSM199")
#adddf=adddf %>% filter(model !="iLJ478")

adddf=adddf[complete.cases(adddf$gene_count), ]
adddf=adddf[complete.cases(adddf$avg), ]
order <- resultWithModel %>% filter(model %in% modelSelector, bof=="general", mediumType=="seed") %>% arrange(gene_count) %>% select(labelSOE) %>% .[[1]]
g4a <- ggplot(resultWithModel %>% filter(model %in% modelSelector, bof=="general", mediumType=="seed"), 
				aes(x=reorder(factor(labelSOE), gene_count), y=add, fill=sgCat)) +
	geom_violin(colour=NA, scale="width") +
	stat_summary(fun.y=mean, geom="point", shape="|", size=5, colour="black") +
	scale_fill_manual(values=c(specialist=mycolor[2], generalist=mycolor[1])) +
	scale_x_discrete(labels=parse(text=modelTable %>% filter(bigg_id %in% modelSelector) %>% arrange(gene_count) %>% select(labelSOE) %>% .[[1]])) +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	ylab("# added reactions") +
	xlab("model - ordered by network size") + 
  CCB_plot_style()+
  
	theme(legend.position = "none",
		axis.text.y = element_text(angle = 0, hjust = 1, size=14)) +
	#ggtitle(bold(a)~"distribution of # added reactions") +
  ggtitle(bold(a)~"") +
	theme(plot.title=element_text(hjust=0, size=15),
			plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm")) +
	coord_flip()
g4a
ggsave(g4a, filename="Fig4aenergyseed.pdf", path="~/plot/", width=220, height=297/2, units="mm")

adddf=adddf[complete.cases(adddf$gene_count), ]
adddf=adddf[complete.cases(adddf$avg), ]
modelNLS <- nls(adddf, formula="avg ~ 1/(2*c*gene_count)", start=list(c=1), control=nls.control(minFactor=1e-10))
adddf$fit <- predict(modelNLS)


#modelLog <- nls(adddf, formula="log(avg) ~ lc - log(gene_count)", start=list(lc=1), control=nls.control(minFactor=1e-10))
#adddf$fitLog <- exp(predict(modelLog))
#modelLm <- lm(adddf, formula="log(avg) ~  offset(log(gene_count))")
#adddf$fitLm <- exp(predict(modelLm))
#cat("fittet value c is:\n")
#print(fitValueC <- exp(-modelLog$m$getPars()["lc"]) / 2)

#adddf <- adddf %>% filter(!auxo, model!="e_coli_core")


modelLm <- lm(adddf %>% filter(model !="iAF1260", model %in% modelSelector ), formula="log(avg) ~  offset(-log(gene_count))")
adddf=adddf[complete.cases(adddf$gene_count), ]
#modelFit <- adddf %>% filter(model !="iAF1260", model %in% modelSelector ) %>% select(model, avg, gene_count)
modelFit <- data.frame(gene_count=min(adddf$gene_count):max(adddf$gene_count),rm.na=T)
#modelLmPred <- as.data.frame(exp(predict(modelLm, interval="confidence")))
modelLmPred <- data.frame(fit=exp(predict(modelLm, newdata=modelFit)))
modelFit$fitLm <- modelLmPred$fit
#modelFit$fitLmLwr <- modelLmPred$lwr
#modelFit$fitLmUpr <- modelLmPred$upr


modelLmF <- lm(adddf %>% filter(model !="iAF1260", model %in% modelSelector ), formula="log(avg) ~  log(gene_count)")
#modelFit$fitLmF <- exp(predict(modelLmF))
modelFit$fitLmF <- exp(predict(modelLmF, newdata=modelFit))


cat("fittet value c is:\n")
print(fitValueC <- exp(-modelLm$coefficients["(Intercept)"]) / 2)
cat("fittet value c for free model is:\n")
print(fitValueCF <- exp(-modelLmF$coefficients["(Intercept)"]) / 2)
cat("factor for free model is:\n")
print(modelLmF$coefficients["log(gene_count)"])
cat("alpha for free model:\n")
calcAlpha <- function(m) -(m -1)
print( calcAlpha(modelLmF$coefficients["log(gene_count)"]))
cat("confidence interval for slope and alpha:\n")
print(confint(modelLmF, 'log(gene_count)', level=0.95))
print(calcAlpha(confint(modelLmF, 'log(gene_count)', level=0.95)))
cat("MSE for normal fit:\n")
print(mean(modelLm$residuals^2))
cat("MSE for free fit:\n")
print(mean(modelLmF$residuals^2))

addReactSpearman <- with(adddf %>% filter(model !="iAF1260", model %in% modelSelector ), cor(x=gene_count, y=avg, method="spearman"))
print(with(adddf %>% filter(model !="iAF1260", model %in% modelSelector ), cor.test(x=gene_count, y=avg, method="spearman")))
#data:  gene_count and avg
#S = 4030, p-value = 3.579e-05
#alternative hypothesis: true rho is not equal to 0
#sample estimates:
#  rho 
#-0.7521739 
g4b <- ggplot(adddf %>% filter(model %in% modelSelector), aes(x=gene_count, y=avg, colour=sgCat,
				shape=I(ifelse(isEco55, 17, 16) - ifelse(auxo, 15, 0)),
				alpha=I(ifelse(model %in% modelSelector, 1, 0.5)))) +
	geom_line(data=modelFit, aes(x=gene_count, y=fitLmF, colour=NULL, shape=NULL, alpha=NULL), size=.7, colour=mycolor[3]) +
	#geom_line(data=modelFit, aes(x=gene_count, y=fitLm, colour=NULL, shape=NULL, alpha=NULL), size=.7, colour=myYellow) +
	geom_point(size=3) +
	scale_colour_manual(name="is E.Coli 55", values=c(mycolor[2],mycolor[1]), breaks=c("specialist", "generalist")) +
	guides(shape=FALSE, alpha=FALSE,
		colour=FALSE
		#colour = guide_legend(order = 1, override.aes = list(linetype = c(1, 1), color = c(myGreen, myYellow))) 
		) +
	xlab("gene count") + 
	ylab("average # added reactions") +
	scale_y_log10(limits=c(1, 162), breaks=c(1, 5, 10,50, 100, 160), minor_breaks=c(1:10, seq(from=10, to=162, by=10))) +
  #for_energy
  #scale_y_log10(limits=c(1, 8), breaks=c(1, 2, 3, 5, 8), minor_breaks=c(1:10, seq(from=10, to=162, by=10))) +
  
  CCB_plot_style()+
  
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 0, hjust = 0.5)) +
	geom_text_repel(aes(label=ifelse(!isEco55 | model %in% c(modelSelector, "e_coli_core"), as.character(model), "")), 
							segment.color = "#cccccc",
							colour="black",
							max.iter=5000,
							size=5,
							max.overlaps = Inf,
							segment.curvature = 0,
							segment.ncp = 3,
							segment.angle = 20,
							point.padding= unit(0.5, 'lines'),
							box.padding = unit(0.5, 'lines')) +
	#ggtitle(bold(b)~"average # added reactions and model size") +
  ggtitle(bold(b)~"") +
  annotate("text", x = 1200, y = 120, vjust=0, label = paste0("Spearman's ", "\u03C1", "= ", format(addReactSpearman, digits=2)), parse=F, size=5,fontface = "bold") +
 #for energy
  #annotate("text", x = 1400, y = 8, vjust=0, label = paste0("Spearman's ", "\u03C1", "= ", format(addReactSpearman, digits=2)), parse=F, size=5,fontface = "bold") +
  
  #annotate("text", x = 1200, y = 100, vjust=0, label = paste0("spearmans correlation = ", format(addReactSpearman, digits=2)), parse=F, size=4) +
	#annotate("text", x = 100, y = 3, hjust=0, label = "f(g) == frac(g^{1-alpha}, alpha*c)", parse=T, size=4) +
	theme(plot.title=element_text(hjust=0, size=15),
			plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
g4b


#

ggsave(g4b, filename="Fig4b.pdf", path="~/plot/",width=220, height=297/2, units="mm")





gS5 <- ggplot(adddf %>% filter(isEco55), aes(x=gene_count, y=avg, shape=auxo, colour=sgCat)) +
	geom_smooth(method="lm", aes(x=gene_count, y=avg, color=NULL, shape=NULL), colour = "#777777", alpha=0.7, method.args = list(), se=F) +
	geom_point() +
	guides(colour=FALSE, shape=FALSE, alpha=FALSE) +
	xlab("gene count") + 
	ylab("average # added reactions") +
	scale_shape_manual(name="is auxotrophic", values=c(17, 2), breaks=c(FALSE, TRUE)) +
	scale_y_continuous(breaks=c(1:10), labels=1:10) +
	scale_colour_manual(values=c(specialist=mycolor[2], generalist=mycolor[1])) +
#	coord_trans(y="log10", limy=c(2, 200)) +
  CCB_plot_style()+
  
	theme(legend.position = "bottom",
		axis.text.x = element_text(angle = 0, hjust = 0.5)) +
	geom_text_repel(aes(label=ifelse(model %in% c("iML1515", "iEcDH1_1363") | auxo, as.character(model), "")), 
							segment.color = "#cccccc",
							size=4.5,
							max.overlaps = Inf,
							colour="black",
							point.padding= unit(0.5, 'lines'),
							box.padding = unit(0.5, 'lines')) +
	#ggtitle(bold(a)~"") +
	annotate("text", x = 1100, y = 5, vjust=0, hjust=0, label = paste0("Spearman's ", "\u03C1", "= ", format(with(adddf %>% filter(isEco55,model!="iEcDH1_1363"), cor(x=gene_count, y=avg, method="spearman")), digits=2)), parse=F, size=5,fontface = "bold") +
	ggtitle(bold("Average # added reactions and model size for 55"~italic(E.coli)~"only")~"") +
	theme(plot.title=element_text(hjust=0, size=15),
			plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
gS5
ggsave(gS5, filename="FigS5.pdf", path="~/plot/", width=210, height=297/2, units="mm")

###
###without auxotrophic
gS51 <- ggplot(adddf %>% filter(isEco55,!(model %in% auxotrophic)), aes(x=gene_count, y=avg, shape=auxo, colour=sgCat)) +
  geom_smooth(method="lm", aes(x=gene_count, y=avg, color=NULL, shape=NULL), colour = "#777777", alpha=0.7, method.args = list(), se=F) +
  geom_point() +
  guides(colour=FALSE, shape=FALSE, alpha=FALSE) +
  xlab("gene count") + 
  ylab("average # added reactions") +
  scale_shape_manual(name="is auxotrophic", values=c(17, 2), breaks=c(FALSE, TRUE)) +
  scale_y_continuous(limits=c(1,4),breaks=c(1:3), labels=1:3) +
  scale_colour_manual(values=c(specialist=myRed, generalist=myBlue)) +
  #	coord_trans(y="log10", limy=c(2, 200)) +
  theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 0, hjust = 0.5)) +
  geom_text_repel(aes(label=ifelse(model %in% c("iML1515", "iEcDH1_1363") | auxo, as.character(model), "")), 
                  segment.color = "#cccccc",
                  max.overlaps = Inf,
                  colour="black",
                  point.padding= unit(0.5, 'lines'),
                  box.padding = unit(0.5, 'lines')) +
  #ggtitle(bold(a)~"") +
  annotate("text", x = 1450, y = 3.5, vjust=0, hjust=0, label = paste0("Spearman correlation = ", format(with(adddf %>% filter(isEco55,!(model %in% auxotrophic)), cor(x=gene_count, y=avg, method="spearman")), digits=2)), parse=F, size=4) +
  ggtitle(bold("Figure S5.")~"average # added reactions and model size for 62"~italic(E.coli)~"only") +
  theme(plot.title=element_text(hjust=0, size=10),
        plot.margin=unit(c(0.5, 0.5, 0.5, 0.5), "cm"))
gS51
ggsave(gS51, filename="FigS51.pdf", path="~/plot/", width=210, height=297/2, units="mm")

####

indexData <- innovationIndex %>% filter(bof=="general", mediumType=="seed") %>% group_by(model, bof, mediumType, isEco55, organism, sgCat, auxo, gene_count) %>% 
		summarize(innoIndexDistinctRelMean = mean(innoIndexDistinctRel),
					innoIndexDistinctMean = mean(innoIndexDistinct, na.rm=T),
					v1=mean(startEnv/distinctEnv), 
					v50rel=mean(startEnv50/distinctEnv),
					v50=mean(startEnv50),
					vMean=mean(startEnvMean),
					meanAdd=mean(addCor)) %>% ungroup()
indexData <- indexData %>% ungroup() %>% left_join(modelTable %>% select(bigg_id, twoReactMetRatio) %>% rename(model=bigg_id))


g4 <- grid.arrange(g4a, g4b, ncol = 1, nrow = 2)
ggsave(g4, filename="Fig4newfeb.pdf", path="~/plot/", width=297/2, height=260, units="mm")






networkLinaritySpearman <- with(modelTable %>% filter(bigg_id !="iAF1260", bigg_id %in% modelSelector ), cor.test(x=gene_count, y=twoReactMetRatio, method="spearman"))
print(networkLinaritySpearman)

g5 <- ggplot(modelTable %>% filter(bigg_id %in% c(modelSelector, "e_coli_core")), aes(x=gene_count, y=twoReactMetRatio, colour=sgCat, 
		shape=I(ifelse(bigg_id == "e_coli_core", 3, ifelse(isEco55, 17, 16) - ifelse(auxo, 15, 0)))
		#alpha=I(ifelse(bigg_id %in% modelSelector, 1, 0.5))
		)) +
	geom_hline(mapping=NULL, yintercept=0.4, colour="#cccccc") +
	geom_point() +
	ylab("network linearity") + 
	xlab("gene count") +
	xlim(c(100, 1710))+
#	annotate("text", x = 800, y = 0.6, vjust=0, label = paste0("spearman correlation = ", 
#						format(cor(x=modelTable$gene_count, y=modelTable$twoReactMetRatio, method="spearman"), digits=2)),#
#						parse=F, size=4) +
	guides(colour=FALSE, shape=FALSE) +
	scale_colour_manual(values=c(specialist=myRed, generalist=myBlue), breaks=c("specialist", "generalist")) +
	geom_text_repel(aes(label=as.character(bigg_id)), colour="black",max.overlaps = Inf, segment.color = "#cccccc",
		 point.padding= unit(0.5, 'lines'), box.padding = unit(0.5, 'lines')) +
  annotate("text", x = 1000, y = 0.6, vjust=0, hjust=0, label = paste0("spearman correlation = ", format(networkLinaritySpearman[["estimate"]], digits=2)), parse=F, size=4)+
	#ggtitle("network structure and gene count") +
	theme(plot.title=element_text(hjust=0, size=10))
g5
ggsave(g5, filename="Fig5.pdf", path="~/plot/", width=297/2, height=210/2, units="mm")


gS6 <- ggplot(resultWithModel %>% filter(bof=="general", mediumType=="seed"), aes(x=reorder(factor(labelSOE), gene_count), y=add, fill=sgCat)) +
	geom_violin(colour=NA, scale="width") +
	stat_summary(fun.y=mean, geom="point", shape="|", size=4, colour="black") +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	#scale_fill_manual(name="is E.Coli 55", values=c(myBlue, myGreen), breaks=c(FALSE, TRUE)) +
	scale_fill_manual(values=c(specialist=mycolor[2], generalist=mycolor[1]), breaks=c("specialist", "generalist")) +
	scale_y_continuous(name="# added reactions", breaks=c(seq(from=0, to=150, by=10))) + #, labels=c(1:5, rep("", 4), 10, rep("", 8), 100) ) +
	scale_x_discrete(label=myParse) +
	xlab("model - ordered by network size") + 
	facet_grid(ecoLabel ~ ., scales="free_y", space="free_y", labeller=label_parsed) +
	ggtitle(bold("Distributions of added reactions per submodel")~"") +
  theme_bw()+
  
  
  theme(
    text = element_text(size=10,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x=element_text(size=10,family="sans",face="bold",colour="#666666"),
    axis.text.y=element_text(size=8.5,family="sans",face="bold",colour="#666666",angle = 0, hjust = 1, vjust=0.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 10, face = "bold", margin = margin(t = 5, r = 10, b = 5, l = 5)),
    legend.box.background = element_rect(colour = "black",linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank() ,
    
    panel.background = element_blank(),
    axis.line = element_line(color = '#666666',linewidth = 0.8),
    axis.ticks =element_line(color="#666666",linewidth =1),
    axis.ticks.length=unit(c(-0.3,0.3), "cm"))+
  
  
  theme(legend.position = "none",
		axis.text.y = element_text(angle = 0, hjust = 1, size=8)) +
	theme(plot.title=element_text(hjust=0, size=15)) +
	coord_flip()
	gS6
ggsave(gS6, filename="FigS6.pdf", path="~/plot/", width=210, height=297, units="mm")



gS7 <- ggplot(innovationIndex %>% filter(bof=="general", mediumType=="random"), aes(x=reorder(factor(labelSOE), gene_count), y=innoIndexRel, fill=sgCat)) +
	geom_violin(colour=NA, scale="width") +
	stat_summary(fun=mean, geom="point", shape="|", size=4, colour="black") +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	scale_fill_manual(values=c(specialist=mycolor[2], generalist=mycolor[1]), breaks=c("specialist", "generalist")) +
	ylab("collateral adaptation index") +
	xlab("model - ordered by network size") + 
	scale_x_discrete(label=myParse) +
	ggtitle(bold("Distributions of the collateral adaptation index per submodel")~"") +
	facet_grid(ecoLabel ~ ., scales="free_y", space="free_y", labeller=label_parsed) +
  theme_bw()+
  
  
  theme(
    text = element_text(size=10,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x=element_text(size=10,family="sans",face="bold",colour="#666666"),
    axis.text.y=element_text(size=8.5,family="sans",face="bold",colour="#666666",angle = 0, hjust = 1, vjust=0.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 10, face = "bold", margin = margin(t = 5, r = 10, b = 5, l = 5)),
    legend.box.background = element_rect(colour = "black",linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank() ,
    
    panel.background = element_blank(),
    axis.line = element_line(color = '#666666',linewidth = 0.8),
    axis.ticks =element_line(color="#666666",linewidth =1),
    axis.ticks.length=unit(c(-0.3,0.3), "cm"))+
  
	theme(legend.position = "none",
		axis.text.y = element_text(angle = 0, hjust = 1, size=8)) +
	theme(plot.title=element_text(hjust=0, size=15)) +
	coord_flip()
gS7
ggsave(gS7, filename="FigS7.pdf", path="~/plot/", width=210, height=297, units="mm")




gS8 <- ggplot(innovationIndex %>% filter(bof=="general", mediumType=="random"), aes(x=reorder(factor(labelSOE), gene_count), y=startEnvMean, fill=sgCat)) +
	geom_violin(colour=NA, scale="width") +
	stat_summary(fun.y=mean, geom="point", shape="|", size=4, colour="black") +
	guides(colour=guide_legend(title.position = "top", ncol=9), 
			shape=guide_legend(title.position = "top", nrow=1)) +
	scale_fill_manual(values=c(specialist=mycolor[2], generalist=mycolor[1]), breaks=c("specialist", "generalist")) +
	ylab("exaptation index") +
	xlab("model - ordered by network size") + 
	scale_x_discrete(label=myParse) +
	ggtitle(bold("Distributions of the exaptation index per submodel")~"") +
	facet_grid(ecoLabel ~ ., scales="free_y", space="free_y", labeller=label_parsed) +
  theme_bw()+
  
  
  theme(
    text = element_text(size=10,face="bold",color="black",family="sans"),
    axis.title.x = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.title.y.left  = element_text(color = "black", size = 14, face = "bold",margin = margin(t = 0, r = 0, b = 0, l = 0)),
    axis.text.x=element_text(size=10,family="sans",face="bold",colour="#666666"),
    axis.text.y=element_text(size=8.5,family="sans",face="bold",colour="#666666",angle = 0, hjust = 1, vjust=0.5),
    plot.margin = margin(10, 10, 10, 10),
    axis.title.y = element_text(color = "black", size = 10, face = "bold", margin = margin(t = 5, r = 10, b = 5, l = 5)),
    legend.box.background = element_rect(colour = "black",linewidth = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank() ,
    
    panel.background = element_blank(),
    axis.line = element_line(color = '#666666',linewidth = 0.8),
    axis.ticks =element_line(color="#666666",linewidth =1),
    axis.ticks.length=unit(c(-0.3,0.3), "cm"))+
  
	theme(legend.position = "none",
		axis.text.y = element_text(angle = 0, hjust = 1, size=8)) +
	theme(plot.title=element_text(hjust=0, size=15)) +
	coord_flip()
gS8
ggsave(gS8, filename="FigS8.pdf", path="~/plot/", width=210, height=297, units="mm")








system("pdfunite plots/Fig{,S}?.pdf plots.pdf")
cat("converting to jpg\n")
file.remove(list.files("jpg/", full.names=T))
system("./makeJpg.sh")






















