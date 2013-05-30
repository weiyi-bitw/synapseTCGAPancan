require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

# RPPA attractor matrices Synapse ID
synIDs <- c(
	"syn1875924",	#BLCA
	"syn1875934",	#BRCA
	"syn1875932",	#COAD
	"syn1875920",	#GBM
	"syn1875922",	#HNSC
	"syn1875926",	#KIRC
	"syn1875930",	#LUAD
	"syn1875914",	#LUSC
	"syn1875936",	#OV
	"syn1875918",	#READ
	"syn1875928"	#UCEC
)

datasetTags <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC", "KIRC", "LUAD", "LUSC", "OV", "READ", "UCEC")

rppaParentID <- "syn1875889"
nc <- length(datasetTags)

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

for(i in 1:nc){
	synid <- synIDs[i]
	syn <- synGet(synid, downloadFile=TRUE, downloadLocation=tmpDir)
	nm <- load(getFileLocation(syn), env)
	x <- env[[nm]]
	fileName <- file.path(resultDir, paste(datasetTags[i], ".attractors.rppa.txt", sep=""))
	out <- outputAttractors(x, strength.pos=3, outputGeneNumber=5, write2File=T, fileName=fileName)
	
	outFile <- File(fileName, synapseStore=TRUE, parentId=rppaParentID)
	used <- list(synIDs[i])
	activity <- Activity(name="Output attractors", used=used, description="strength.pos=3, outputGeneNumber=5")

	outFile <- synStore(outFile, activity=activity)


}
