require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/clusterAttractorsFromSynapse.R")
clusterLink <- getPermlink(analysisRepo, "sourceScripts/clusterAttractorsFromSynapse.R")

sourceRepoFile(analysisRepo, "sourceScripts/createConsensusFromSynapse.R")
consensusLink <- getPermlink(analysisRepo, "sourceScripts/createConsensusFromSynapse.R")


# miRNA attractor matrices Synapse ID
synIDs <- c(
	"syn1875924",	#BLCA
	"syn1875934",	#BRCA
	"syn1875932",	#COAD
	"syn1875920",	#GBM
	"syn1875922",	#HNSC
	"syn1875926",	#KIRC
	"syn1875930",	#LUAD
	"syn1875914",	#LUSC
	"syn1875936",	#OV
	"syn1875918",	#READ
	"syn1875928"	#UCEC
)

datasetTags <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "READ", "UCEC")

rppaParentID <- "syn1875889"

used <- list(
	list(url=clusterLink, name=basename(clusterLink), wasExecuted=TRUE),
	list(url=consensusLink, name=basename(consensusLink), wasExecuted=TRUE)
)

for(s in synIDs){
  used <- c(used, s)
}


resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

# clustering attractors
alist <- clusterAttractorsFromSynapse(synIDs, numGenes=5, strength.pos=3, datasetTags=datasetTags, tempDir=tmpDir)

# creating consensus miRNA ranking from attractor clusters
synIDList <- as.list(synIDs)
names(synIDList) <- datasetTags
attractome <- createConsensusFromSynapse(alist, synIDList, minGenes=3, tempDir=tmpDir)


# save files and upload to Synapse
resultFile.cluster <- file.path(resultDir, "attractorClusters.rppa.rda")
save(alist, file=resultFile.cluster)
resultFile.consensus <- file.path(resultDir, "attractome.rppa.rda")
save(attractome, file=resultFile.consensus)

activity <- Activity(name="Attractor clustering", used=used, description="numGenes=5, strength.pos=3")

resultFile.cluster <- File(resultFile.cluster, synapseStore=TRUE, parentId=rppaParentID)
generatedBy(resultFile.cluster) <- activity
resultFile.cluster <- storeEntity(resultFile.cluster)
activity <- generatedBy(resultFile.cluster)

resultFile.consensus <- File(resultFile.consensus, synapseStore=TRUE, parentId=rppaParentID)
generatedBy(resultFile.consensus) <- activity
resultFile.consensus <- storeEntity(resultFile.consensus)
activity <- generatedBy(resultFile.consensus)

