generateScatterPlots <- function(gene1, gene2, gene3){
	fileName <- paste("scatter.", gene1, "x", gene2, "x", gene3,".png", sep="")
	fig <- path.file(figDir, fileName)
	
	png(fileName, width = 7.3, height = 8, units = "in", res = 300, pointsize=12)
	
	par(mar = c(4,4,2,5),       #plot margin
	mfrow = c(4, 3),
	oma=c(0, 0, 0, 0),
	mgp=c(2, 1, 0))


	for(i in 1:length(data)){
		ge = data[[i]]
		map <- cbind(rownames(ge))
		rownames(map) <- rownames(ge)
		z <- ge[gene3,]
		z[is.infinite(z)] <- NA
		z <- z-median(z, na.rm=T)
		iqr.z <- quantile(z, 0.75, na.rm=T) - quantile(z, 0.25, na.rm=T)
		z[z > (quantile(z, 0.75, na.rm=T) + 1.5*iqr.z)] <- NA
		z[z < (quantile(z, 0.25, na.rm=T) - 1.5*iqr.z)] <- NA

		x <- ge[gene1,]
		x[is.infinite(x)] <- NA
		x <- x-median(x, na.rm=T)
		iqr.x <- quantile(x, 0.75, na.rm=T) - quantile(x, 0.25, na.rm=T)
		x[x > (quantile(x, 0.75, na.rm=T) + 1.5*iqr.x)] <- NA
		x[x < (quantile(x, 0.25, na.rm=T) - 1.5*iqr.x)] <- NA

		y <- ge[gene2,]
		y[is.infinite(y)] <- NA
		y <- y-median(y, na.rm=T)
		iqr.y <- quantile(y, 0.75, na.rm=T) - quantile(y, 0.25, na.rm=T)
		y[y > (quantile(y, 0.75, na.rm=T) + 1.5*iqr.y)] <- NA
		y[y < (quantile(y, 0.25, na.rm=T) - 1.5*iqr.y)] <- NA

		idx <- !is.na(x) & !is.na(y) & !is.na(z)
		cc <- coltransform(z[idx])
	
		plot(x[idx], y[idx], pch=16, col=rgb(cc), 
			xlim = c(min(-3, min(x[idx])), max(3, max(x[idx]))), ylim=c(min(-3, min(y[idx])), max(3, max(y[idx]))), 
			main = pancan[i] , xlab=gene1, ylab=gene2, cex.axis=0.9, cex.label=1)
		mtext(gene3, side=4, line=0.35, cex=0.7)
		image.plot( legend.only=TRUE, legend.width=0.8, legend.mar=2, zlim= range(z[idx]), col = two.colors(start="blue", end="red", middle="gray"), axis.args=list(cex.axis=0.9)) 
	}	

	dev.off()       #Write
	return (fig)
}
