
# remove entrez ID from gene identifiers in RNASeq data
getGeneSymbols = function(innames){
	outnames <- sapply(innames, function(x){
	
		if(regexpr("\\?", x) > 0){
			o = strsplit(x, "\\|")[[1]][2]
		}else{
			o = strsplit(x, "\\|")[[1]][1]
		}
		return (o)

	}
	)
	return (outnames)
}
