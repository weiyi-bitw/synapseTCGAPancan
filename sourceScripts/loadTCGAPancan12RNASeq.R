require(synapseClient)
require(DreamBox7)

data(attractome.minimalist)
genes <- unlist(sapply(attractome.minimalist, function(x){x[,1]}))

# pancan 12 tags
pancan <- c("BLCA","BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")

synIDs <- c(
	"syn1571504", #BLCA
	"syn417812", #BRCA
	"syn1446197", #COAD
	"syn1446214", #GBM
	"syn1571420", #HNSC
	"syn417925", #KIRC
	"syn1571538", #LAML
	"syn1571468", #LUAD
	"syn418033", #LUSC
	"syn1446264", #OV
	"syn1446276", #READ
	"syn1446289" #UCEC
)

nf <- length(synIDs)
data <- list()

for(i in 1:nf){
	synid <- synIDs[i]
	cat("Processing",synid,"...\n" );flush.console()
	syn <- loadEntity(synid)
	ge <- load.exp(file.path(syn$cacheDir, syn$files[[1]]))
	oo <- getGeneSymbols(rownames(ge))
	rownames(ge) <- oo

	ge <- log2(ge+1)

	data[[i]] <- ge[intersect(genes, rownames(ge)),]
	names(data)[i] <- pancan[i]
}
