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
	"syn1876029",	#BLCA
	"syn1876039",	#BRCA
	"syn1876037",	#COAD
	"syn1876023",	#GBM
	"syn1876027",	#HNSC
	"syn1876031",	#KIRC
	"syn1876025", 	#LAML
	"syn1876035",	#LUAD
	"syn1876019",	#LUSC
	"syn1876040",	#OV
	"syn1876021",	#READ
	"syn1876033"	#UCEC
)

datasetTags <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC","KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")

mrnaParentID <- "syn1875937"

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
alist <- clusterAttractorsFromSynapse(synIDs, numGenes=20, strength.pos=10, datasetTags=datasetTags, tempDir=tmpDir)

# creating consensus miRNA ranking from attractor clusters
synIDList <- as.list(synIDs)
names(synIDList) <- datasetTags
attractome <- createConsensusFromSynapse(alist, synIDList, minGenes=10, tempDir=tmpDir)


# save files and upload to Synapse
resultFile.cluster <- file.path(resultDir, "attractorClusters.rnaseq.rda")
save(alist, file=resultFile.cluster)
resultFile.consensus <- file.path(resultDir, "attractome.rnaseq.rda")
save(attractome, file=resultFile.consensus)

activity <- Activity(name="Attractor clustering", used=used, description="numGenes=20, strength.pos=10")

resultFile.cluster <- File(resultFile.cluster, synapseStore=TRUE, parentId=mrnaParentID)
generatedBy(resultFile.cluster) <- activity
resultFile.cluster <- storeEntity(resultFile.cluster)
activity <- generatedBy(resultFile.cluster)

resultFile.consensus <- File(resultFile.consensus, synapseStore=TRUE, parentId=mrnaParentID)
generatedBy(resultFile.consensus) <- activity
resultFile.consensus <- storeEntity(resultFile.consensus)
activity <- generatedBy(resultFile.consensus)

