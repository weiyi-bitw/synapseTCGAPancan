createFigure2 <- function(meta.pancan, x, y, z){

	fileName <- paste("scatter.", x, "x", y, "x", z, ".png", sep="")

	png(fileName, width=7.3, height=8, units="in", res=300, pointsize=12)
	par(mar = c(4,4,2,5),       #plot margin
	  mfrow = c(4, 3),
	  oma=c(0, 0, 0, 0),
	  mgp=c(2, 1, 0)
	  )

	# find the features
	temp <- meta.pancan[[1]]
	idxx <- NULL
	idxy <- NULL
	idxz <- NULL
	for(d in names(temp)){
		if(x %in% rownames(temp[[d]])) idxx <- d
		if(y %in% rownames(temp[[d]])) idxy <- d
		if(z %in% rownames(temp[[d]])) idxz <- d
	}

	if(is.null(idxx) | is.null(idxy) | is.null(idxz)) stop("Cannot find all features!")

	for(c in names(meta.pancan)){
	  if(is.na(meta.pancan[[c]])){
	  }else{
	    message("Processing ", c, "...")
	    temp <- meta.pancan[[c]]

		dx <- temp[[idxx]]
		dy <- temp[[idxy]]
		dz <- temp[[idxz]]

	    commonid <- intersect(intersect(colnames(temp[[idxx]]), colnames(temp[[idxy]])), colnames(temp[[idxz]]) )
		if(length(commonid) < 10) next

		dx <- dx[,commonid]
		dy <- dy[,commonid]
		dz <- dz[,commonid]

	    cc <- coltransform2(dz[z,])
	    plot(dx[x,], dy[y,], pch=16, col=rgb(cc), 
		main = c , xlab=x, ylab=y, cex.axis=0.9, cex.label=1)
		mtext(z, side=4, line=0.35, cex=0.7)
	    midpoint <- (median(dz[z,])-min(dz[z,])) / diff(range(dz[z,])) 
			image.plot( legend.only=TRUE, legend.width=0.8, legend.mar=2, zlim= range(dz[z,]),  
			col = colorMapping3(xbreak=c(min(dz[z,]), median(dz[z,]), max(dz[z,]))), axis.args=list(cex.axis=0.9)) 
	  }
	}
	dev.off()

	return (fileName)
}
