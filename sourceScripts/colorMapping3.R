colorMapping3 <- function(xbreak, cbreak=list( c(0, 0, 1), c(192/255, 192/255, 192/255), c(1, 0, 0)), n=256){
	if(length(xbreak) != length(cbreak)) stop("Color breaks have different length than value breaks!")
	cout <- matrix(NA, nrow=n, ncol=3)
	xbreak <- sort(xbreak)
	midp <- floor( (xbreak[2] - xbreak[1]) / (xbreak[3] - xbreak[1]) * n )
	cout[1,] <- cbreak[[1]]
	if(midp!=1){
		cout[midp,] <- cbreak[[2]]
	}
	cout[midp+1,] <- cbreak[[2]]
	cout[n, ] <- cbreak[[3]]

	if(midp-1 >= 2){
		for(i in 2:(midp-1)){
			cout[i, 1] <- cbreak[[1]][1] +  (i - 1)/(midp-1) * (cbreak[[2]][1] - cbreak[[1]][1])
			cout[i, 2] <- cbreak[[1]][2] +  (i - 1)/(midp-1) * (cbreak[[2]][2] - cbreak[[1]][2])
			cout[i, 3] <- cbreak[[1]][3] +  (i - 1)/(midp-1) * (cbreak[[2]][3] - cbreak[[1]][3])
		}
	}

	if(n >= midp+3 ){
		for(i in (midp+2) : (n-1)){
			cout[i, 1] <- cbreak[[2]][1] +  (i - midp - 1)/(n - midp) * (cbreak[[3]][1] - cbreak[[2]][1])
			cout[i, 2] <- cbreak[[2]][2] +  (i - midp - 1)/(n - midp) * (cbreak[[3]][2] - cbreak[[2]][2])
			cout[i, 3] <- cbreak[[2]][3] +  (i - midp - 1)/(n - midp) * (cbreak[[3]][3] - cbreak[[2]][3])

		}	
	}

	cc <- apply(cout, 1, function(x){rgb(x[1], x[2], x[3])})
	return (cc)
}
