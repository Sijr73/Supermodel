.PHONY : all
.SUFFIXES:
.SECONDARY: 
MAKEFLAGS += --no-builtin-rules
CC=
CXX=
SHELL:=/bin/bash -o pipefail

CORES=7

# task dependency files (*.d) are orderd sequential.
# targets below are ordered in reverse
all: 	download.d\
		readModels.d\
		modelsCheck.d\
		buildModel.d\
		analyseUniversalModel.d\
		convertMedia.d\
		mediaMapping/mediaBiGG.Rdata\
		models/universalBiGG.ver1.2.Rdata\
		models/universalBiGG.ver1.1.Rdata\
		models/universalBiGG.ver1.Rdata\
		compGrowth.d\
		models/universalMNX.ver2.Rdata\
		imbalanceRemovalByBlocking.d\
		massBalanceCheck.d\
		modelTable.d\
		analyseEGC.d\
		analyseFBA.d\
		analyseARM.d\
		analyseCross.d\
		innovationIndex.Rdata

clean: rm -f *.d


#.d: .d
#	@echo ">>> "
#	@rm -f .out
#	exit 1 |& tee -a ../.out
#	touch .d

################################################################################
innovationIndex.Rdata: analyseCross.d
analyseCross.d: analyseARM.d envirDist/crossdata/*.Rdata
	@echo ">>> analysing cross fba data"
	@rm -f analyseCross.out
	cd envirDist && ./analyseCross.R |& tee -a ../analyseCross.out
	touch analyseCross.d
################################################################################
mergedResults.Rdata result.Rdata: analyseARM.d
analyseARM.d: envirDist/compareDf.Rdata envirDist/envirdata/*.Rdata modelTable.d
	@echo ">>> analysing envirDist data and creating plots"
	@rm -f analyseARM.out
	cd envirDist && ./analyseARM.R |& tee -a ../analyseARM.out
	touch analyseARM.d
################################################################################
modelTable.d: dataAnalysis/tax_report.txt sourceData/modelList.tsv dataAnalysis/updateModelTable.R
	@echo ">>> updating model table"
	@rm -f modelTable.out
	cd dataAnalysis && ./updateModelTable.R |& tee -a ../modelTable.out
	touch modelTable.d
################################################################################
envirDist/compareDf.Rdata: analyseFBA.d
analyseFBA.d: analyseEGC.d models/universalMNX.ver2.Rdata envirDist/analyseFBA.R
	@echo ">>> analysing fba results"
	@rm -f analyseFBA.out
	cd envirDist && ./analyseFBA.R |& tee -a ../analyseFBA.out
	touch analyseFBA.d
################################################################################
# egcRemoval has to be done on a compute node, here only results will be summarized
models/universalMNX.ver2.Rdata: analyseEGC.d qc/testGrowth.R
analyseEGC.d: imbalanceRemovalByBlocking.d models/universalBiGG.ver1.2.Rdata energyGeneratingCycleRemoval/results/*.Rdata energyGeneratingCycleRemoval/analyseResults.R compGrowth.d
	@echo ">>> analysing EGC removal data"
	@rm -f analyseEGC.out testGrowth.out
	@rm -f energyGeneratingCycleRemoval/*.pdf
	cd energyGeneratingCycleRemoval && ./analyseResults.R |& tee -a ../.out
	cd qc && ./testGrowth.R |& tee -a ../testGrowth.out
	touch analyseEGC.d
################################################################################
# mass balances removal by finding blocked reactions.
imbalanceRemovalByBlocking.d: massBalanceCheck.d models/universalBiGG.ver1.1.Rdata imbalanceRemoval/imbalanceRemovalByBlocking.R
	@echo ">>> mass imbalances removal by finding blocked reactions"
	@rm -f imbalanceRemovalByBlocking.out
	cd imbalanceRemoval && ./imbalanceRemovalByBlocking.R |& tee -a ../imbalanceRemovalByBlocking.out
	touch imbalanceRemovalByBlocking.d
################################################################################
# mass balance check
massBalanceCheck.d: convertMedia.d models/universalBiGG.ver1.1.Rdata massBalanceCheck/massBalanceCheck.R
	@echo ">>> mass balance check"
	@rm -f massBalanceCheck.out
	cd massBalanceCheck && ./massBalanceCheck.R |& tee -a ../massBalanceCheck.out
	touch massBalanceCheck.d

################################################################################
# compare growth of previous seed adaption to latest bigg version
compGrowth.d: convertMedia.d analyseUniversalModel.d convertMedia/mediaBiGG.Rdata models/universalBiGG.ver1.1.Rdata
	@echo ">>> comparing bigg model to growth from previous model"
	@rm -f compGrowth.out
	cd qc && ./compGrowth.R |& tee -a ../.out
	touch compGrowth.d
################################################################################
#map media from seed namespace to BiGG namespace
mediaMapping/mediaBiGG.Rdata: convertMedia.d
convertMedia.d: analyseUniversalModel.d\
				models/universalBiGG.ver1.1.Rdata\
				convertMedia/convertMedia.R
	@echo ">>> converting media from seed to bigg"
	@rm -f convertMedia.out
	cd convertMedia && ./convertMedia.R |& tee -a ../convertMedia.out
	cd convertMedia && ./generateRandomMedia.R |& tee -a ../convertMedia.out
	cd convertMedia && ./mdsMedia.R |& tee -a ../convertMedia.out
	touch convertMedia.d
################################################################################
#analysemodel
models/universalBiGG.ver1.1.Rdata: analyseUniversalModel.d
analyseUniversalModel.d: 	buildModel/analyseUniversalModel.R\
							models/universalBiGG.ver1.Rdata\
							buildModel/selectBiomassTable_mod.csv\
							buildModel.d
	@echo ">>> running analyseUniversalModel"
	@rm -f analyseUniversalModel.out
	cd buildModel && ./analyseUniversalModel.R |& tee -a ../analyseUniversalModel.out
	touch analyseUniversalModel.d

################################################################################
# build model
models/universalBiGG.ver1.Rdata: buildModel.d
buildModel.d: 	buildModel/addDefaultBOF.R\
				buildModel/buildModel.R\
				modelsCheck.d
	@echo ">>> running buildModel"
	@rm -f buildModel.out
	cd buildModel && ./buildModel.R |& tee -a ../buildModel.out
	touch buildModel.d

################################################################################
modelsCheck.d: 	modelsCheck/modelsCheck.R\
				sourceData/models.Rdata\
				sourceData/excludeModels.txt\
				download.d\
				readModels.d
	@echo ">>> running modelsCheck"
	@rm -f modelsCheck.out
	cd modelsCheck && exec ./modelsCheck.R |& tee -a ../modelsCheck.out
	touch modelsCheck.d
	

################################################################################
sourceData/models.Rdata: readModels.d

readModels.d: download.d
	@echo ">>> reading models"
	rm -rf readModels.out
	cd modelsCheck && ./readModels.R |& tee ../readModels.out
	touch readModels.d

################################################################################

download.d:	
	@echo ">>> fetching models from Bigg Database"
	cd download && ./getModels.py |& tee ../sourceData/modelList_source.tsv
	cd sourceData && cat modelList_source.tsv buchnera.tsv > modelList.tsv
	cd sourceData && wget -N http://bigg.ucsd.edu/static/namespace/bigg_models_reactions.txt
	cd sourceData && wget -N http://bigg.ucsd.edu/static/namespace/bigg_models_metabolites.txt
	touch download.d
	



















