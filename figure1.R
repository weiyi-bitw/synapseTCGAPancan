require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/getGeneSymbols.R")
getGeneSymbolLink <- getPermlink(analysisRepo, "sourceScripts/getGeneSymbols.R")
sourceRepoFile(analysisRepo, "sourceScripts/coltransform2.R")
coltransformLink <- getPermlink(analysisRepo, "sourceScripts/coltransform2.R")
sourceRepoFile(analysisRepo, "sourceScripts/createFigure1.R")
scatterRNASeqLink <- getPermlink(analysisRepo, "sourceScripts/createFigure1.R")

# load pancan syn ID table
syn <- synGet("syn1875837", downloadFile=T, downloadLocation=tmpDir)
pancanTable <- loadClin(getFileLocation(syn))
nf <- nrow(pancanTable)

# load RNASeq attractome
syn <- synGet("syn1876552", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
attractome <- env[[nm]]

resultDir <- file.path(tmpDir, "scatter")
dir.create(resultDir)

scatterParentID <- "syn1899337"
used <- list(
	list(url=scatterRNASeqLink, name=basename(scatterRNASeqLink), wasExecuted=TRUE),
	list(url=getGeneSymbolLink, name=basename(getGeneSymbolLink), wasExecuted=TRUE),
	list(url=coltransformLink, name=basename(coltransformLink), wasExecuted=TRUE),
	"syn1875837",
	"syn1876552"
)

synIDs <- pancanTable[,"RNASeq"]
names(synIDs) <- rownames(pancanTable)
synIDs <- synIDs[order(names(synIDs))]
idx <- !is.na(synIDs)
synIDs <- synIDs[idx]
synIDList <- as.list(synIDs)

for(s in synIDs){
  used <- c(used, s)
}

geneList <- lapply(attractome, function(a){a[1:3,1]})
geneList <- geneList[1:4]
fList <- createFigure1(synIDList, geneList, filePath=resultDir)

nf <- length(fList)

activity <- Activity(name="Create scatter plots", used=used)
for(i in 1:nf){
	syn <- File(fList[[i]], synapseStore=TRUE, parentId = scatterParentID)
	generatedBy(syn) <- activity
	syn <- storeEntity(syn)
	activity <- generatedBy(syn)
}




