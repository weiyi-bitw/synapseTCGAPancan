require(cafr)

createAttractomeGLRNASeq <- function(alist){
	k <- length(alist)
	attractome <- list()
	data(genome)
	rownames(genome) <- paste(rownames(genome), "|", gsub(" ", "", genome[,1]), sep="")
	for(i in 1:k){
		a <- alist[[i]]
		if(class(a) != "AttractorSet") break
		temp <- a$getConsensus(10)
		temp2 <- rbind(c(names(temp), ""))
		ctbd <- findCytoband(temp2, genome)[1,1]
		temp <- cbind(names(temp), temp)
		rownames(temp) <- NULL
		colnames(temp) <- c("Gene.Symbol", "Score")
		attractome[[i]] <- temp
		names(attractome)[i] <- ctbd
	}
	return (attractome)
}
