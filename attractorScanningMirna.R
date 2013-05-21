require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/mirnaAttractorScanning.R")
mirnaAttractorScanningLink <- getPermlink(analysisRepo, "sourceScripts/mirnaAttractorScanning.R")

# load pancan syn ID table
syn <- synGet("syn1875837", downloadFile=T, downloadLocation=tmpDir)
pancanTable <- loadClin(file.path(tmpDir, syn$properties$name))
nf <- nrow(pancanTable)

# load miRNA stem name map
syn <- synGet("syn1875840", downloadFile=T, downloadLocation=tmpDir)
nm <- load(file.path(tmpDir, syn$properties$name), env)
map <- env[[nm]]

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

miRNAParentID <- "syn1875845"
used <- list(
	list(url=mirnaAttractorScanningLink, name=basename(mirnaAttractorScanningLink), wasExecuted=TRUE),
	"syn1875840"
)


for(i in 1:nf){
	syn <- pancanTable[i, "miRNA"]
	if(is.na(syn) | rownames(pancanTable)[i]=="GBM") next
	x <- mirnaAttractorScanning(syn, map, tmpDir)
	resultFile <- file.path(resultDir, paste(rownames(pancanTable)[i], ".attractorMatrix.mirna.rda", sep=""))
	save(x, file=resultFile)
	activity <- Activity(name="miRNA attractor scanning", used=c(used, syn))
	resultFile <- File(resultFile, synapseStore=TRUE, parentId=miRNAParentID)
	generatedBy(resultFile) <- activity
	resultFile <- storeEntity(resultFile)
	activity <- generatedBy(resultFile)
}


