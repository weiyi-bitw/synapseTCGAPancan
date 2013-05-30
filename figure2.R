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
sourceRepoFile(analysisRepo, "sourceScripts/createFigure2.R")
scatterRNASeqLink <- getPermlink(analysisRepo, "sourceScripts/createFigure2.R")

# load attractome profile
syn <- synGet("syn1899334", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
meta.pancan <- env[[nm]]

x <- "M-"
y <- "M+"
z <- "LYM"

resultDir <- file.path(tmpDir, "scatter")
dir.create(resultDir)

scatterParentID <- "syn1899337"
used <- list(
	list(url=scatterRNASeqLink, name=basename(scatterRNASeqLink), wasExecuted=TRUE),
	list(url=getGeneSymbolLink, name=basename(getGeneSymbolLink), wasExecuted=TRUE),
	list(url=coltransformLink, name=basename(coltransformLink), wasExecuted=TRUE),
	"syn1899334"
)

f <- createFigure2(meta.pancan, x, y, z)

activity <- Activity(name="Create scatter plots", used=used)
syn <- File(f, synapseStore=TRUE, parentId = scatterParentID)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)



