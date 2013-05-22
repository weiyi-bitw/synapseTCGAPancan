coltransform2 <- function(data){
	dataL <- length(data)
	data <- data - median(data)

	rg <- diff(range(data))
	cs <- sapply(data, function(x){
		if(is.na(x)) return( NA )
		if(x >= 0){
			return (x / max(data))
		}else{
			return (x / (-min(data)) )
		}
	})
	cc <- array(NA, dim=c(dataL, 3))
	mid <- 192/255

	for (i in 1:dataL){
		if(is.na(cs[i])) next
		if(cs[i] > 0){
		cc[i, 1] <- mid + cs[i] * (1-mid)
		cc[i, 2] <- mid - cs[i]*mid
		cc[i, 3] <- cc[i, 2]
		} else{
		cc[i, 3] <- mid - cs[i] * (1-mid)
		cc[i, 2] <- mid + cs[i] * mid
		cc[i, 1] <- cc[i, 2]
		
		}
	
	}
	cc
}
