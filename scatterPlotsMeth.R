require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/coltransform2.R")
coltransformLink <- getPermlink(analysisRepo, "sourceScripts/coltransform2.R")
sourceRepoFile(analysisRepo, "sourceScripts/createScatterMeth.R")
scatterMethLink <- getPermlink(analysisRepo, "sourceScripts/createScatterMeth.R")
sourceRepoFile(analysisRepo, "sourceScripts/colorMapping3.R")
colmapLink <- getPermlink(analysisRepo, "sourceScripts/colorMapping3.R")

# load pancan syn ID table
syn <- synGet("syn1875837", downloadFile=T, downloadLocation=tmpDir)
pancanTable <- loadClin(getFileLocation(syn))
nf <- nrow(pancanTable)

# load methylation probe map
syn <- synGet("syn1875841", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
map <- env[[nm]]

# load methylation attractome
syn <- synGet("syn1876325", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
attractome <- env[[nm]]

resultDir <- file.path(tmpDir, "scatter")
dir.create(resultDir)

scatterParentID <- "syn1759352"
used <- list(
	list(url=scatterMethLink, name=basename(scatterMethLink), wasExecuted=TRUE),
	list(url=coltransformLink, name=basename(coltransformLink), wasExecuted=TRUE),
	list(url=colmapLink, name=basename(colmapLink), wasExecuted=TRUE),
	"syn1875837",
	"syn1875841",
	"syn1876325"
)

synIDs <- pancanTable[,"methylation_27k"]
names(synIDs) <- rownames(pancanTable)
synIDs <- synIDs[order(names(synIDs))]
synIDs['BLCA'] <- 'syn1889358'
synIDs['HNSC'] <- 'syn1889356'
synIDList <- as.list(synIDs)

for(s in synIDs){
  used <- c(used, s)
}

geneList <- lapply(attractome, function(a){a[1:3,1]})
fList <- createScatterMeth(synIDList, geneList, map, filePath=resultDir)

nf <- length(fList)

activity <- Activity(name="Create scatter plots", used=used)
for(i in 1:nf){
	syn <- File(fList[[i]], synapseStore=TRUE, parentId = scatterParentID)
	generatedBy(syn) <- activity
	syn <- storeEntity(syn)
	activity <- generatedBy(syn)
}




