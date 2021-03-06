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


# meth attractor matrices Synapse ID
synIDs <- c(
	"syn1895440",	#BLCA
	"syn1876059",	#BRCA
	"syn1876058",	#COAD
	"syn1876048",	#GBM
	"syn1895920",	#HNSC
	"syn1876052",	#KIRC
	"syn1876050", 	#LAML
	"syn1876056",	#LUAD
	"syn1876044",	#LUSC
	"syn1876061",	#OV
	"syn1876046",	#READ
	"syn1876054"	#UCEC
)

datasetTags <- c("BLCA","BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")

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
alist <- clusterAttractorsFromSynapse(synIDs, numGenes=50, strength.pos=10, datasetTags=datasetTags, tempDir=tmpDir)

# creating consensus miRNA ranking from attractor clusters
synIDList <- as.list(synIDs)
names(synIDList) <- datasetTags
attractome <- createConsensusFromSynapse(alist, synIDList, minGenes=10, tempDir=tmpDir)

# remove X chromosome attractor
alist[[1]] <- NULL
attractome[[1]] <- NULL

names(attractome)[2:3] <- c("M+", "M-")

# removing weak attractor from attractor clusters using attractome data
topGenes <- unlist(lapply(attractome, function(x){x[1,1]}))
idx <-  unlist(lapply( lapply(alist, function(x){ names(x$getConsensus(10))}) , function(x){sum(topGenes %in% x)}))
alist <- alist[idx==1]

# save files and upload to Synapse
resultFile.cluster <- file.path(resultDir, "attractorClusters.meth.rda")
save(alist, file=resultFile.cluster)
resultFile.consensus <- file.path(resultDir, "attractome.meth.rda")
save(attractome, file=resultFile.consensus)

activity <- Activity(name="Attractor clustering", used=used, description="numGenes=50, strength.pos=10, *remove sex chromosome attractors")

resultFile.cluster <- File(resultFile.cluster, synapseStore=TRUE, parentId=methParentID)
generatedBy(resultFile.cluster) <- activity
resultFile.cluster <- storeEntity(resultFile.cluster)
activity <- generatedBy(resultFile.cluster)

resultFile.consensus <- File(resultFile.consensus, synapseStore=TRUE, parentId=methParentID)
generatedBy(resultFile.consensus) <- activity
resultFile.consensus <- storeEntity(resultFile.consensus)
activity <- generatedBy(resultFile.consensus)

