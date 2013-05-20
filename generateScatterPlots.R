require(DreamBox7)
require(synapseClient)
require(fields)
require(rGithubClient)

synapseLogin()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/getGeneSymbols.R")
getGeneSymbolLink <- getPermlink(analysisRepo, "sourceScripts/getGeneSymbols.R")

sourceRepoFile(analysisRepo, "sourceScripts/coltransform.R")
coltransformLink <- getPermlink(analysisRepo, "sourceScripts/coltransform.R")

sourceRepoFile(analysisRepo, "sourceScripts/pancan12Scatter.R")
scatterLink <- getPermlink(analysisRepo, "sourceScripts/pancan12Scatter.R")

sourceRepoFile(analysisRepo, "sourceScripts/loadTCGAPancan12RNASeq.R")
loadDataLink <- getPermlink(analysisRepo, "sourceScripts/loadTCGAPancan12RNASeq.R")

figDir <- file.path(tempdir(), "scatterPlots")
dir.create(figDir)

synScatterFolder <- "syn1759352"

#==================================================
#
# Figure 1: Top three genes of the CIN attractor metagene in pancan 12
#
#==================================================


gene1 = "CENPA"
gene2 = "DLGAP5"
gene3 = "MELK"

fig_cin <- pancan12Scatter(gene1, gene2, gene3)

#==================================================
#
# Figure 2: Top three genes of the MES attractor metagene in pancan 12
#
#==================================================


gene1 = "COL5A2"
gene2 = "SPARC"
gene3 = "VCAN"

fig_mes <- pancan12Scatter(gene1, gene2, gene3)

#==================================================
#
# Figure 3: Top three genes of the LYM attractor metagene in pancan 12
#
#==================================================

gene1 = "PTPRC"
gene2 = "CD53"
gene3 = "LCP2"

fig_lym <- pancan12Scatter(gene1, gene2, gene3)

#==================================================
#
# Figure 4: Top three genes of the chr8q24.3 attractor metagene in pancan 12
#
#==================================================

gene1 = "EXOSC4"
gene2 = "PUF60"
gene3 = "BOP1"

fig_8q24_3 <- pancan12Scatter(gene1, gene2, gene3)

used <- list(
	list(url=getGeneSymbolLink, name=basename(getGeneSymbolLink), wasExecuted=TRUE),
	list(url=coltransformLink, name=basename(coltransformLink), wasExecuted=TRUE),
	list(url=scatterLink, name=basename(scatterLink), wasExecuted=TRUE),
	list(url=loadDataLink, name=basename(loadDataLink), wasExecuted=TRUE)
)

for(s in syn){
	used[[length(used)+1]] <- list(entity=s, wasExecuted=F)
}


activity <- Activity(name = "Scatter plots for top genes in attractor metagenes",
			used=used)

activity <- storeEntity(activity)

figFile <- File(fig_cin, synapseStore=TRUE, parentId=synScatterFolder)
generatedBy(figFile) <- activity
figFile <- storeEntity(figFile)
activity <- generatedBy(figFile)

figFile <- File(fig_mes, synapseStore=TRUE, parentId=synScatterFolder)
generatedBy(figFile) <- activity
figFile <- storeEntity(figFile)
activity <- generatedBy(figFile)

figFile <- File(fig_lym, synapseStore=TRUE, parentId=synScatterFolder)
generatedBy(figFile) <- activity
figFile <- storeEntity(figFile)
activity <- generatedBy(figFile)

figFile <- File(fig_8q24_3, synapseStore=TRUE, parentId=synScatterFolder)
generatedBy(figFile) <- activity
figFile <- storeEntity(figFile)
activity <- generatedBy(figFile)

