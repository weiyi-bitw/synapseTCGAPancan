require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

analysisRepo <- getRepo("weiyi-bitw/synapseTCGAPancan")

sourceRepoFile(analysisRepo, "sourceScripts/findCytoband.R")
cytobandLink <- getPermlink(analysisRepo, "sourceScripts/findCytoband.R")
data(genome)
rownames(genome) <- paste(rownames(genome), "|", gsub(" ", "", genome[,1]), sep="")

# load GL mRNA attractor clusters
syn <- synGet("syn1897839", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/glAttractorClusters.rnaseq.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	if(class(a) != "AttractorSet") break
	mat <- a$getGeneMatrix(10)
	mat <- findCytoband(mat, genome)
	write.table(mat, file=outFile, sep='\t', quote=F,  col.names=F, append=T)
	gt = a$getGeneTable(5)
	gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common genes"
	write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
	write("\n", file=outFile, append=T)
}

used <- list(
	"syn1897839",
	list(url=cytobandLink, name=basename(cytobandLink), wasExecuted=TRUE)
)
syn <- File(outFile, synapseStore=TRUE, parentId="syn1897732")
activity <- Activity(name="Output clusters", used=used, description="getGeneMatrix(10), geneTable(5)")
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


