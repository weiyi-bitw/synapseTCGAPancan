require(cafr)

resortClusters <- function(alist, minSize, minOverlap=NULL){
	sizes <- unlist(lapply(alist, function(x){if(class(x)=="Attractor") return (1); return(length(x$attractors))}))
	alist[sizes < minSize] = NULL
	na <- length(alist)
	for(i in 1:na){
		aa = alist[[i]]
		if(class(aa)!="AttractorSet") next
		scores = as.numeric(as.vector(alist[[i]]$getGeneMatrix(5)[,6]))
		alist[[i]]$medStrength <- mean(scores[order(scores, decreasing=T)[1:minSize]])
	}
	scores <- unlist(lapply(alist, function(x){x$medStrength} ))
	alist <- alist[order(scores, decreasing=T)]

	if(!is.null(minOverlap)){
		topOvlp = unlist(lapply(alist, function(x){if(class(x)=="Attractor") return(1); return(x$getGeneTable(1)[1])}))
		alist[topOvlp < minOverlap] = NULL
	}
	
	return (alist)
}

