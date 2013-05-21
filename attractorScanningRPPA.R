require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/rppaAttractorScanning.R")
rppaAttractorScanningLink <- getPermlink(analysisRepo, "sourceScripts/rppaAttractorScanning.R")

# load pancan syn ID table
syn <- synGet("syn1875837", downloadFile=T, downloadLocation=tmpDir)
pancanTable <- loadClin(file.path(tmpDir, syn$properties$name))
nf <- nrow(pancanTable)

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

rppaParentID <- "syn1875889"
used <- list(
	list(url=rppaAttractorScanningLink, name=basename(rppaAttractorScanningLink), wasExecuted=TRUE)
)

for(i in 4:nf){
	syn <- pancanTable[i, "RPPA"]
	if(is.na(syn)) next
	x <- rppaAttractorScanning(syn, map, tmpDir)
	resultFile <- file.path(resultDir, paste(rownames(pancanTable)[i], ".attractorMatrix.rppa.rda", sep=""))
	save(x, file=resultFile)
	activity <- Activity(name="Protein attractor scanning", used=c(used, syn))
	resultFile <- File(resultFile, synapseStore=TRUE, parentId=rppaParentID)
	generatedBy(resultFile) <- activity
	resultFile <- storeEntity(resultFile)
	activity <- generatedBy(resultFile)
}


