require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempDir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/clusterAttractorsFromSynapse.R")
clusterLink <- getPermlink(analysisRepo, "sourceScripts/clusterAttractorsFromSynapse.R")

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
	list(url=clusterLink, name=basename(clusterLink), wasExecuted=TRUE),
	as.list(synIDs)
)

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

alist <- clusterAttractorsFromSynapse(synIDs, numGenes=5, strength.pos=3, datasetTags=datasetTags, tempDir=tmpDir)

resultFile <- file.path(resultDir, "attractorClusters.mirna.rda")
save(alist, file=resultFile)
activity <- Activity(name="Attractor clustering, NG=5, Str=3rd MI", used=used)
resultFile <- File(resultFile, synapseStore=TRUE, parentId=miRNAParentID)
generatedBy(resultFile) <- activity
resultFile <- storeEntity(resultFile)
activity <- generatedBy(resultFile)


