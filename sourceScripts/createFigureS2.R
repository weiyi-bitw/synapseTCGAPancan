createFigureS2 <- function(meta.pancan, x, y){

	fileName <- paste("scatter.", x, "x", y, ".png", sep="")

	png(fileName, width=7.3, height=9, units="in", res=300, pointsize=12)
	par(mar = c(4,4,2,1),       #plot margin
	  mfrow = c(4, 3),
	  oma=c(0, 0, 0, 0),
	  mgp=c(2, 1, 0)
	  )

	# find the features
	temp <- meta.pancan[[1]]
	idxx <- NULL
	idxy <- NULL
	for(d in names(temp)){
		if(x %in% rownames(temp[[d]])) idxx <- d
		if(y %in% rownames(temp[[d]])) idxy <- d
	}

	if(is.null(idxx) | is.null(idxy)) stop("Cannot find all features!")

	for(c in names(meta.pancan)){
	  if(is.na(meta.pancan[[c]])){
	  }else{
	    message("Processing ", c, "...")
	    temp <- meta.pancan[[c]]

		dx <- temp[[idxx]]
		dy <- temp[[idxy]]

	    commonid <- intersect(colnames(temp[[idxx]]), colnames(temp[[idxy]])) 
		if(length(commonid) < 10) next

		dx <- dx[,commonid]
		dy <- dy[,commonid]

	    plot(dx[x,], dy[y,], pch=16, col="blue", 
		main = c , xlab=x, ylab=y, cex.axis=0.8)
	  }
	}
	dev.off()

	return (fileName)
}
