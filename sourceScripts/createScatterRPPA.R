require(synapseClient)
require(cafr)
require(fields)

createScatterRPPA <- function(synIDList, geneList, filePath="./"){
	nf <- length(synIDs)
	data <- list()
	syn <- list()
	figList <- list()
	pancan <- names(synIDList)
	for(i in 1:nf){
		synid <- synIDList[[i]]
		tag <- names(synIDList)[i]
		syn[[tag]] <- loadEntity(synid)
		message("Processing ", tag, " (", synid, ") ...\n")
		ge <- loadExpr(file.path(syn[[tag]]$cacheDir, syn[[tag]]$files[[1]]))

		data[[tag]] <- ge[intersect(unlist(geneList), rownames(ge)),]
	}

	ng <- length(geneList)
	n <- floor(sqrt(nf))
	m <- ceiling(nf / n)

	for(j in 1:ng){
		gt <- geneList[[j]]
		gene1 <- gt[1]
		gene2 <- gt[2]
		gene3 <- gt[3]
		fileName <- paste("scatter.", gene1, "x", gene2, "x", gene3, ".png", sep="")
		fig <- file.path(filePath, fileName)
		png(fig, width = 7.3, height = 8, units = "in", res = 300, pointsize=12)
	
		par(mar = c(4,4,2,5),       #plot margin
		mfrow = c(4, 3),
		oma=c(0, 0, 0, 0),
		mgp=c(2, 1, 0))


		for(i in 1:length(data)){
			ge = data[[i]]
			map <- cbind(rownames(ge))
			rownames(map) <- rownames(ge)
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
				xlim = c(min(-3, min(x[idx])), max(3, max(x[idx]))), ylim=c(min(-3, min(y[idx])), max(3, max(y[idx]))), 
				main = pancan[i] , xlab=gene1, ylab=gene2, cex.axis=0.9, cex.label=1)
			mtext(gene3, side=4, line=0.35, cex=0.7)
			image.plot( legend.only=TRUE, legend.width=0.8, legend.mar=2, zlim= range(z[idx]),  
				col = colorMapping3(xbreak=c(min(z[idx]), median(z[idx]), max(z[idx]))), axis.args=list(cex.axis=0.9)) 
		}	

		dev.off()       #Write
		figList[[j]] <- fig
	}
	return (figList)
}
