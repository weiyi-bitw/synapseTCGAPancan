require(cafr)
require(impute)
require(limma)
require(synapseClient)

mirnaAttractorScanning <- function(synID, map, tmpDir=tempDir()){
	syn <- synGet(synID, downloadFile=TRUE, downloadLocation=tmpDir)
	mirna <- loadExpr(file.path(tmpDir, syn$files[[1]]))
	nz <- apply(mirna, 1, function(x){sum(x==0)})
	mirna <- mirna[nz < 0.5*(ncol(mirna)), ]

  # impute zero counts and missing values
	mirna[mirna==0] <- NA
	mirna <- log2(mirna)
	mirna <- impute.knn(mirna)$data
  
  # normalize expression values using quantile normalization
	mirna <- normalizeBetweenArrays(mirna)

  # summarize mirnas with the same family names into the same feature
	mirna <- probeSummarization(ge=mirna, map=map, threshold=0.7, gene.colname="miR_stem")
	
	x <- attractorScanning(mirna, a = 5, maxIter=500)
	return (x)
}



