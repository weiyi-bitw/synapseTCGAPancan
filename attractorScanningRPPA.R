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
	list(url=rppaAttractorScanningLink, name=basename(mirnaAttractorScanningLink), wasExecuted=TRUE),
)

for(i in 1:nf){
	syn <- pancanTable[i, "RPPA"]
	x <- mirnaAttractorScanning(syn, map, tmpDir)
	resultFile <- file.path(resultDir, paste(rownames(pancanTable)[i], ".attractorMatrix.rda", sep=""))
	save(x, file=resultFile)
	activity <- Activity(name="Protein attractor scanning", used=c(used, syn))
	resultFile <- File(resultFile, synapseStore=TRUE, parentId=miRNAParentID)
	generatedBy(resultFile) <- activity
	resultFile <- storeEntity(resultFile)
	activity <- generatedBy(resultFile)
}

fileNames = list.files(path="attractors/", pattern="*.rda")

alist = clusterAttractors(filePath="attractors",fileNames=fileNames, numGenes=5, strength.pos=3, datasetTags=tags)

k = length(alist)
outFile = "attractors/matchTable.txt"
attractome = list()
for(i in 1:k){
	a = alist[[i]]
	if(class(a) != "AttractorSet") break
	write.table(a$getGeneMatrix(10), file=outFile, sep='\t', quote=F,  col.names=F, append=T)
	gt = a$getGeneTable(10)
	gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common genes"
	write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
	write("\n", file=outFile, append=T)
	temp = a$getConsensus(5)
	temp = cbind(names(temp), temp)
	rownames(temp) = NULL
	colnames(temp) = c("Gene.Symbol", "Score")
	attractome[[i]] = temp
	names(attractome)[i] = strsplit(temp[1,1],split="\\|")[[1]][1]
	
}

for(f in fileNames){
	load(paste("attractors/", f, sep=""))
	token = strsplit(f, "\\.")[[1]]
	mat <- outputAttractors(x, strength.pos=5, outputGeneNumber=10, write2File=TRUE, fileName=paste(token[1], ".", token[2], ".rppa.txt", sep=""))
}


