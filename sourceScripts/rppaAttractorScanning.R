require(cafr)
require(impute)
require(limma)
require(synapseClient)

rppaAttractorScanning <- function(synID, map, tmpDir=tempDir()){
	syn <- loadEntity(synID)
	rppa <- loadExpr(file.path(syn$cacheDir, syn$files[[1]]))
  
	x <- attractorScanning(rppa, a = 2, maxIter=500)
	return (x)
}



