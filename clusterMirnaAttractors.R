require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempDir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/clusterAttractorsFromSynapse.R")
clusterLink <- getPermlink(analysisRepo, "sourceScripts/clusterAttractorsFromSynapse.R")

sourceRepoFile(analysisRepo, "sourceScripts/createConsensusFromSynapse.R")
consensusLink <- getPermlink(analysisRepo, "sourceScripts/createConsensusFromSynapse.R")


# miRNA attractor matrices Synapse ID
synIDs <- c(
	"syn1875859",	#BLCA
	"syn1875869",	#BRCA
	"syn1875867",	#COAD
	"syn1875857",	#HNSC
	"syn1875861",	#KIRC
	"syn1875854",	#LAML
	"syn1875865",	#LUAD
	"syn1875850",	#LUSC
	"syn1875871",	#OV
	"syn1875852",	#READ
	"syn1875863"	#UCEC
)

datasetTags <- c("BLCA", "BRCA", "COAD", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")

miRNAParentID <- "syn1875845"

used <- list(
	list(url=clusterLink, name=basename(clusterLink), wasExecuted=TRUE)
	list(url=consensusLink, name=basename(consensusLink), wasExecuted=TRUE)
)

for(s in synIDs){
  used <- c(used, s)
}


resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

alist <- clusterAttractorsFromSynapse(synIDs, numGenes=5, strength.pos=3, datasetTags=datasetTags, tempDir=tmpDir)

synIDList <- as.list(synIDs)
names(synIDList) <- datasetTags
attractome <- createConsensusFromSynapse(alist, synIDList, minGenes=3, tempDir=tmpDir)

resultFile.cluster <- file.path(resultDir, "attractorClusters.mirna.rda")
save(alist, file=resultFile.cluster)
resultFile.consensus <- file.path(resultDir, "attractome.mirna.rda")
save(attractome, file=resultFile.consensus)

activity <- Activity(name="Attractor clustering, NG=5, Str=3rd MI", used=used)

resultFile.cluster <- File(resultFile.cluster, synapseStore=TRUE, parentId=miRNAParentID)
generatedBy(resultFile.cluster) <- activity
resultFile.cluster <- storeEntity(resultFile.cluster)
activity <- generatedBy(resultFile.cluster)

resultFile.consensus <- File(resultFile.consensus, synapseStore=TRUE, parentId=miRNAParentID)
generatedBy(resultFile.consensus) <- activity
resultFile.consensus <- storeEntity(resultFile.consensus)
activity <- generatedBy(resultFile.consensus)

