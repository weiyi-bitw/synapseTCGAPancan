require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/createFigureS2.R")
scatterRNASeqLink <- getPermlink(analysisRepo, "sourceScripts/createFigureS2.R")

# load attractome profile
syn <- synGet("syn1899334", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
meta.pancan <- env[[nm]]

x <- "MES"
y <- "END"

resultDir <- file.path(tmpDir, "scatter")
dir.create(resultDir)

scatterParentID <- "syn1899337"
used <- list(
	list(url=scatterRNASeqLink, name=basename(scatterRNASeqLink), wasExecuted=TRUE),
	"syn1899334"
)

f <- createFigureS2(meta.pancan, x, y)

activity <- Activity(name="Create scatter plots", used=used)
syn <- File(f, synapseStore=TRUE, parentId = scatterParentID)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)



