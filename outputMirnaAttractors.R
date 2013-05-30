require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

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
nc <- length(datasetTags)

resultDir <- file.path(tmpDir, "attractors")
dir.create(resultDir)

for(i in 1:nc){
	synid <- synIDs[i]
	syn <- synGet(synid, downloadFile=TRUE, downloadLocation=tmpDir)
	nm <- load(getFileLocation(syn), env)
	x <- env[[nm]]
	fileName <- file.path(resultDir, paste(datasetTags[i], ".attractors.mirna.txt", sep=""))
	out <- outputAttractors(x, strength.pos=3, outputGeneNumber=5, write2File=T, fileName=fileName)
	
	outFile <- File(fileName, synapseStore=TRUE, parentId=miRNAParentID)
	used <- list(synIDs[i])
	activity <- Activity(name="Output attractors", used=used, description="strength.pos=3, outputGeneNumber=5")

	outFile <- synStore(outFile, activity=activity)


}
