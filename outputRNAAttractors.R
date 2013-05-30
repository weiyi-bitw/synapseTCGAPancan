require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

# mRNA attractor matrices Synapse ID
synIDs <- c(
	"syn1876029",	#BLCA
	"syn1876039",	#BRCA
	"syn1876037",	#COAD
	"syn1876023",	#GBM
	"syn1876027",	#HNSC
	"syn1876031",	#KIRC
	"syn1876025", 	#LAML
	"syn1876035",	#LUAD
	"syn1876019",	#LUSC
	"syn1876040",	#OV
	"syn1876021",	#READ
	"syn1876033"	#UCEC
)

datasetTags <- c("BLCA", "BRCA", "COAD", "GBM", "HNSC","KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")
nc <- length(datasetTags)
mrnaParentID <- "syn1875937"

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

for(i in 1:nc){
	synid <- synIDs[i]
	syn <- synGet(synid, downloadFile=TRUE, downloadLocation=tmpDir)
	nm <- load(getFileLocation(syn), env)
	x <- env[[nm]]
	fileName <- file.path(resultDir, paste(datasetTags[i], ".attractors.rnaseq.txt", sep=""))
	out <- outputAttractors(x, strength.pos=10, outputGeneNumber=20, write2File=T, fileName=fileName)
	
	outFile <- File(fileName, synapseStore=TRUE, parentId=mrnaParentID)
	used <- list(synIDs[i])
	activity <- Activity(name="Output attractors", used=used, description="strength.pos=10, outputGeneNumber=20")

	outFile <- synStore(outFile, activity=activity)


}


