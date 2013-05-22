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
	"syn1876059",	#BRCA
	"syn1876058",	#COAD
	"syn1876048",	#GBM
	"syn1876052",	#KIRC
	"syn1876050", 	#LAML
	"syn1876056",	#LUAD
	"syn1876044",	#LUSC
	"syn1876061",	#OV
	"syn1876046",	#READ
	"syn1876054"	#UCEC
)

datasetTags <- c("BRCA", "COAD", "GBM","KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")

methParentID <- "syn1876042"

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
alist <- clusterAttractorsFromSynapse(synIDs, numGenes=100, strength.pos=10, datasetTags=datasetTags, tempDir=tmpDir)

# creating consensus miRNA ranking from attractor clusters
synIDList <- as.list(synIDs)
names(synIDList) <- datasetTags
attractome <- createConsensusFromSynapse(alist, synIDList, minGenes=10, tempDir=tmpDir)


# save files and upload to Synapse
resultFile.cluster <- file.path(resultDir, "attractorClusters.meth.rda")
save(alist, file=resultFile.cluster)
resultFile.consensus <- file.path(resultDir, "attractome.meth.rda")
save(attractome, file=resultFile.consensus)

activity <- Activity(name="Attractor clustering", used=used, description="numGenes=100, strength.pos=10")

resultFile.cluster <- File(resultFile.cluster, synapseStore=TRUE, parentId=methParentID)
generatedBy(resultFile.cluster) <- activity
resultFile.cluster <- storeEntity(resultFile.cluster)
activity <- generatedBy(resultFile.cluster)

resultFile.consensus <- File(resultFile.consensus, synapseStore=TRUE, parentId=methParentID)
generatedBy(resultFile.consensus) <- activity
resultFile.consensus <- storeEntity(resultFile.consensus)
activity <- generatedBy(resultFile.consensus)

