#!/usr/bin/Rscript
library(sybil)

load("~/workspace/models/iSM199/iSM199_e_compart.Rdata")
iSM199 <- upgradeModelorg(iSM199)

source("../helper/loadPropTables.R")

munique <- unique(metProp$universal_bigg_id)


m <- gsub("\\[[cep]\\]", "", met_id(iSM199))

m1 <- gsub("_D_D$", "_DD", m)
m2 <- gsub("_([LDMRS]$)", "__\\1", m1)

m3 <- gsub("^_", "", m2)

mod_compart(iSM199) <- gsub("\\W", "", mod_compart(iSM199))

met_id(iSM199) <- paste0(m3, "[", mod_compart(iSM199)[met_comp(iSM199)], "]" )


S(iSM199)[met_id(iSM199)=="Biomass[c]", react_id(iSM199)=="VGRO"] <- 0

iSM199 <- rmReact(iSM199, "EX_Biomass[c]", rm_met=T)

react_id(iSM199)[react_id(iSM199)=="VGRO"] <- "Biomass"


buchnera <- iSM199

save(buchnera, file="buchnera.Rdata")














