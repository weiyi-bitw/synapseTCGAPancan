
# transform expression value into color coding
coltransform <- function(data){
	dataL <- length(data)
	data <- data - median(data)
	rg = diff(range(data))
	cs = (data-min(data))/rg
	cc <- array(NA, dim=c(dataL, 3))
	mid <- 192/255

	for (i in 1:dataL){
		if(is.na(cs[i])) next
		if(cs[i] > 0.5){
		cc[i, 1] <- 2*cs[i] - 1 + 2*(1-cs[i])*mid
		cc[i, 2] <- 2*mid*(1-cs[i])
		cc[i, 3] <- cc[i, 2]
		
		} else{
		cc[i, 3] <- 1+2*cs[i]*(mid-1)
		cc[i, 2] <- 2*cs[i]*mid
		cc[i, 1] <- cc[i, 2]
		
		}
	
	}
	cc
}
