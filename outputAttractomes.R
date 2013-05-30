require(synapseClient)
require(cafr)
require(rGithubClient)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

# load meth attractor clusters
syn <- synGet("syn1876325", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

# load methylation probe map
syn <- synGet("syn1875841", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
map.meth <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractome.meth.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	probe = rownames(map.meth)[sapply(a[,1], function(x){which(map.meth[,1]==x)})]
	a = cbind(probe, a)
	write(paste(names(alist)[i]), file=outFile, append=T)
	write.table(a, file=outFile, sep='\t', quote=F,  append=T)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876325", "syn1875841")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1876042")
activity <- Activity(name="Output attractome", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load mRNA attractor clusters
syn <- synGet("syn1876552", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractome.rnaseq.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	write(paste(names(alist)[i]), file=outFile, append=T)
	write.table(a, file=outFile, sep='\t', quote=F,  append=T)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876552")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1875937")
activity <- Activity(name="Output attractome", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load miRNA attractor clusters
syn <- synGet("syn1876073", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractome.mirna.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	write(paste(names(alist)[i]), file=outFile, append=T)
	write.table(a, file=outFile, sep='\t', quote=F,  append=T)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876073")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1875845")
activity <- Activity(name="Output attractome", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load rppa attractor clusters
syn <- synGet("syn1876317", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
alist <- env[[nm]]

k <- length(alist)
outFile <- paste(tmpDir, "/attractome.rppa.txt", sep="")
if(file.exists(outFile)) unlink(outFile)

for(i in 1:k){
	a = alist[[i]]
	write(paste(names(alist)[i]), file=outFile, append=T)
	write.table(a, file=outFile, sep='\t', quote=F,  append=T)
	write("\n", file=outFile, append=T)
}

used <- list("syn1876317")
syn <- File(outFile, synapseStore=TRUE, parentId="syn1875889")
activity <- Activity(name="Output attractome", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


