require(synapseClient)
require(cafr)

synapseLogin()

tmpDir <- tempdir()
env <- new.env()

syn <- synGet("syn1899359", downloadFile=T, downloadLocation=tmpDir)

used <- list(
	"syn1876791",
	"syn1876793",
	"syn1876789",
	"syn1876795",
	"syn1876797",
	"syn1876799",
	"syn1876801",
	"syn1876804",
	"syn1876805",
	"syn1876806",
	"syn1895926",
	"syn1895927",
	"syn1895928",
	"syn1876909",
	"syn1876910"
)

activity <- Activity(name="Combine plots", used=used)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)

