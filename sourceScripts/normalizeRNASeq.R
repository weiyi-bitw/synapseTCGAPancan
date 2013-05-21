require(cafr)
require(impute)
require(limma)
require(synapseClient)

normalizeRNASeq <- function(synID, tmpDir=tempDir()){
	syn <- loadEntity(synID)
	ge <- loadExpr(file.path(syn$cacheDir, syn$files[[1]]))
	nz <- apply(ge, 1, function(x){sum(x==0)})
	ge <- ge[nz < 0.5*(ncol(ge)), ]

  # impute zero counts and missing values
	ge[ge==0] <- NA
	ge <- log2(ge)
	ge <- impute.knn(ge)$data
  
  # normalize expression values using quantile normalization
	ge <- normalizeBetweenArrays(ge)

	return (ge)
}
