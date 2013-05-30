require(synapseClient)
require(cafr)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

# Table S1
syn <- synGet("syn1899444", downloadFile=F)

used <- list(
	"syn1899392",
	"syn1899407",
	"syn1899396",
	"syn1899411"
)

activity <- Activity(name="Combine tables", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)

# Table S2

syn <- synGet("syn1899445", downloadFile=F)

used <- list(
	"syn1899420",
	"syn1899422",
	"syn1899424",
	"syn1899426"
)

activity <- Activity(name="Combine tables", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)

# Table S3

syn <- synGet("syn1899446", downloadFile=F)
used <- list("syn1899443")

activity <- Activity(name="xls conversion", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


