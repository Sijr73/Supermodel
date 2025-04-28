#!/usr/bin/Rscript
library(methods)
library(dplyr)
library(XML)

modelTable <- read.csv("../envirDist/modelTable.csv")


xmlURL <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=taxonomy&retmode=xml&id=",
					paste(modelTable$taxid, collapse=","))
xmlText <- paste(scan(xmlURL, what="",sep="\n"),"\n",collapse="\n")
doc <- xmlParse(xmlText)


getNodeSet(doc, "/TaxaSet/Taxon/Lineage")


modelTableTaxonInformation <- data.frame(bigg_id=modelTable$bigg_id,
taxid=xpathSApply(doc, "/TaxaSet/Taxon/TaxId", xmlValue),
name=xpathSApply(doc, "/TaxaSet/Taxon/ScientificName", xmlValue),
division=xpathSApply(doc, "/TaxaSet/Taxon/Division", xmlValue),
lineage=xpathSApply(doc, "/TaxaSet/Taxon/Lineage", xmlValue)
)












