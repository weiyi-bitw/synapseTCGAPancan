require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

# load meth attractor clusters
syn <- synGet("syn1876324", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractorClusters.meth.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	if(class(a) != "AttractorSet") break
	write.table(a$getGeneMatrix(50), file=outFile, sep='\t', quote=F,  col.names=F, append=T)
	gt = a$getGeneTable(10)
	gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common methylation sites"
	write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876324")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1876042")
activity <- Activity(name="Output clusters", used=used, description="getGeneMatrix(50), geneTable(10)")
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load mRNA attractor clusters
syn <- synGet("syn1876551", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractorClusters.rnaseq.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	if(class(a) != "AttractorSet") break
	write.table(a$getGeneMatrix(20), file=outFile, sep='\t', quote=F,  col.names=F, append=T)
	gt = a$getGeneTable(10)
	gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common genes"
	write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876551")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1875937")
activity <- Activity(name="Output clusters", used=used, description="getGeneMatrix(20), geneTable(10)")
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load miRNA attractor clusters
syn <- synGet("syn1876072", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractorClusters.mirna.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	if(class(a) != "AttractorSet") break
	write.table(a$getGeneMatrix(5), file=outFile, sep='\t', quote=F,  col.names=F, append=T)
	gt = a$getGeneTable(3)
	gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common miRNAs"
	write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876072")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1875845")
activity <- Activity(name="Output clusters", used=used, description="getGeneMatrix(5), geneTable(3)")
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load rppa attractor clusters
syn <- synGet("syn1876316", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractorClusters.rppa.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	if(class(a) != "AttractorSet") break
	write.table(a$getGeneMatrix(5), file=outFile, sep='\t', quote=F,  col.names=F, append=T)
	gt = a$getGeneTable(3)
	gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common proteins"
	write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876316")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1875889")
activity <- Activity(name="Output clusters", used=used, description="getGeneMatrix(5), geneTable(3)")
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


