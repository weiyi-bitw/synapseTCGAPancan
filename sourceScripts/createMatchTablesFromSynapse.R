require(cafr)
require(synapseClient)

createMatchTablesFromSynapse <- function(synID, file="matchTable.txt", numGenes=20, geneTable=10){
	syn <- synGet(synID, downloadFile=TRUE)
	env <- new.env()
	nm <- load(getFileLocation(syn), env)
	alist <- env[[nm]]
	k <- length(alist)

	if(file.exists(file)) unlink(file)

	for(i in 1:k){
		a = alist[[i]]
		if(class(a) != "AttractorSet") break
		write.table(a$getGeneMatrix(numGenes), file=outFile, sep='\t', quote=F,  col.names=F, append=T)
		gt = a$getGeneTable(geneTable)
		gt = matrix(paste(names(gt), ":", gt), nrow=1); rownames(gt)="Common genes"
		write.table(gt, file=outFile, sep="\t", append=T, quote=F, col.names=F)
		write("\n", file=outFile, append=T)
	}
	
}
