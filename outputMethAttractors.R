require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

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

nc <- length(datasetTags)
methParentID <- "syn1876042"

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

for(i in 1:nc){
	synid <- synIDs[i]
	syn <- synGet(synid, downloadFile=TRUE, downloadLocation=tmpDir)
	nm <- load(getFileLocation(syn), env)
	x <- env[[nm]]
	fileName <- file.path(resultDir, paste(datasetTags[i], ".attractors.meth.txt", sep=""))
	out <- outputAttractors(x, strength.pos=10, outputGeneNumber=50, write2File=T, fileName=fileName)
	
	outFile <- File(fileName, synapseStore=TRUE, parentId=methParentID)
	used <- list(synIDs[i])
	activity <- Activity(name="Output attractors", used=used, description="strength.pos=10, outputGeneNumber=50")

	outFile <- synStore(outFile, activity=activity)


}
