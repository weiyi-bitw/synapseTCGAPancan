require(synapseClient)
require(cafr)
require(fields)


createScatterMeth <- function(synIDList, geneList, map, filePath="./"){
	nf <- length(synIDs)
	data <- list()
	syn <- list()
	figList <- list()
	pancan <- names(synIDList)
	geneList <- lapply(geneList, function(x){names(x) <- rownames(map)[sapply(x, function(xx){which(map[,1]==xx)} )];x})
	for(i in 1:nf){
		synid <- synIDList[[i]]
		tag <- names(synIDList)[i]
		syn[[tag]] <- synGet(synid, downloadFile=TRUE)
		message("Processing ", tag, " (", synid, ") ...\n")
		if(class(syn[[tag]])=="File"){
			ge <- loadExpr(getFileLocation(syn[[tag]]))
		}else{
			syn[[tag]] <- loadEntity(synid)
			ge <- loadExpr(file.path(syn[[tag]]$cacheDir, syn[[tag]]$files[[1]]))
		}

		data[[tag]] <- ge[intersect(unlist(lapply(geneList, names)), rownames(ge)),]
	}

	ng <- length(geneList)
	n <- floor(sqrt(nf))
	m <- ceiling(nf / n)

	for(j in 1:ng){
		gt <- geneList[[j]]
		gene1 <- names(gt)[1]
		gene2 <- names(gt)[2]
		gene3 <- names(gt)[3]
		fileName <- paste("scatter.", gene1, "x", gene2, "x", gene3, ".png", sep="")
		fig <- file.path(filePath, fileName)
		png(fig, width = 7.3, height = 8, units = "in", res = 300, pointsize=12)
	
		par(mar = c(4,4,2,5),       #plot margin
		mfrow = c(4, 3),
		oma=c(0, 0, 0, 0),
		mgp=c(2, 1, 0))


		for(i in 1:length(data)){
			ge = data[[i]]
			if(!prod(c(gene1, gene2, gene3) %in% rownames(ge))) next
			z <- ge[gene3,]
			z[is.infinite(z)] <- NA
			z <- z-median(z, na.rm=T)

			x <- ge[gene1,]
			x[is.infinite(x)] <- NA
			x <- x-median(x, na.rm=T)

			y <- ge[gene2,]
			y[is.infinite(y)] <- NA
			y <- y-median(y, na.rm=T)

			idx <- !is.na(x) & !is.na(y) & !is.na(z)
			cc <- coltransform2(z[idx])
	
			plot(x[idx], y[idx], pch=16, col=rgb(cc), 
				main = pancan[i] , xlab=paste(gene1, " (", map[gene1,1], ")", sep=""), ylab=paste(gene2, " (", map[gene2, 1], ")", sep=""), cex.axis=0.9, cex.lab=0.7)
			mtext(paste(gene3, " (", map[gene3, 1], ")", sep=""), side=4, line=0.35, cex=0.5)
			midpoint <- (median(z[idx])-min(z[idx])) / diff(range(z[idx])) 
			image.plot( legend.only=TRUE, legend.width=0.8, legend.mar=2, zlim= range(z[idx]),  
				col = colorMapping3(xbreak=c(min(z[idx]), median(z[idx]), max(z[idx]))), axis.args=list(cex.axis=0.9)) 
		}	

		dev.off()       #Write
		figList[[j]] <- fig
	}
	return (figList)
}
