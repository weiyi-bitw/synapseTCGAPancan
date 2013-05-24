require(cafr)
require(synapseClient)

createConsensusFromSynapse <- function(alist, synIDList, minScore=0.5, minGenes=10, th=6, tempDir = tempdir()){
  attractome <- list()
  cMat <- list()
  env <- new.env()
  nas <- length(alist)
  cnt <- 1
  synList <- list()
  nf <- length(synIDList)
  message("Downloading files from Synapse...")
  for(i in 1:nf){
    syn <- synGet(synIDList[[i]], downloadFile=TRUE, downloadLocation=tempDir)
    synList[[ names(synIDList)[i] ]] <- syn
  }


  for(as in alist){
    message("Processing ", cnt , " / ", nas, " attractor cluster...")
    if(class(as) != "AttractorSet") {
    	cnt <- cnt + 1
	next
    }
    clust <- as$attractors
    na <- length(clust)
    if(na < th) next
    atlist <- list()
    for(c in names(clust)){
      tmp <- load(getFileLocation(synList[[c]]), env)[1]
      x <- env[[tmp]]
      topgene <- names(clust[[c]]$genes)[1]
      idx <- which.max(x[,topgene])
      atlist[[c]] <- x[idx,]
    }
    genes <- unique(unlist(lapply(atlist, function(x){names(x)})))
    ng <- length(genes)
    mat <- matrix(NA, nrow=ng, ncol=length(atlist))
    rownames(mat) <- genes
    colnames(mat) <- names(clust)

    for(c in names(atlist)){
      mat[names(atlist[[c]]), c] = atlist[[c]]
    }

    overlapSum <- unlist(lapply(as$attractors, function(a){sum(unlist(lapply(as$attractors, function(aa){aa$getOverlapNum(a)})))}))
    overlapSum = sort(overlapSum, decreasing=T)
    if(ncol(mat) > 2){
      chosen <- names(overlapSum)[1:th]
    }else{
      chosen <- 1:ncol(mat)
    }
    cssMI <- apply(mat[,chosen], 1, mean)
    cssMI <- sort(cssMI, decreasing=T)
    cssMI <- cbind(cssMI[cssMI >= minScore])
    if(length(cssMI) < minGenes) {
	cnt <- cnt + 1
	next
    }
    cssMI <- cbind(rownames(cssMI), cssMI)
    rownames(cssMI) <- NULL
    colnames(cssMI) <- c("Gene.Symbol", "Score")
    attractome[[length(attractome)+1]] <- cssMI
    names(attractome)[[length(attractome)]] <- cssMI[1,1]
    cnt <- cnt + 1
  }
  return (attractome)

}


