require(cafr)
require(synapseClient)
require(limma)
require(impute)

synapseLogin()

tags <- c("BLCA","BRCA", "COAD", "GBM", "HNSC", "KIRC", "LAML", "LUAD", "LUSC", "OV", "READ", "UCEC")

tmpDir <- tempdir()
env <- new.env()


# load RNASeq attractome
syn <- synGet("syn1876552", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
attractome.rnaseq <- env[[nm]]


# load methylation probe map
syn <- synGet("syn1875841", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
map.meth <- env[[nm]]
# load methylation attractome
syn <- synGet("syn1876325", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
attractome.meth <- env[[nm]]


# load miRNA stem name map
syn <- synGet("syn1875840", downloadFile=T, downloadLocation=tmpDir)
nm <- load(file.path(tmpDir, syn$properties$name), env)
map.mirna <- env[[nm]]
# load miRNA attractome
syn <- synGet("syn1876073", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
attractome.mirna <- env[[nm]]


# load RPPA attractome
syn <- synGet("syn1876317", downloadFile=T, downloadLocation=tmpDir)
nm <- load(getFileLocation(syn), env)
attractome.rppa <- env[[nm]]


allattractome = list(rnaseq = attractome.rnaseq, mirna = attractome.mirna, rppa = attractome.rppa, meth = attractome.meth)
pancanParentID <- "syn1876760"
used <- list(
  "syn1876552",
  "syn1876325",
  "syn1875840",
  "syn1876317"
)
attractomeFile <- file.path(tmpDir, "allattractome.rda")
save(allattractome, file=attractomeFile)
activity <- Activity(name="Attractome, ASSEMBLE!!", used=used)
syn <- File(attractomeFile, synapseStore=TRUE, parentId=pancanParentID)
generatedBy(syn) <- activity
syn <- storeEntity(syn)
activity <- generatedBy(syn)


# load pancan syn ID table
syn <- synGet("syn1875837", downloadFile=T, downloadLocation=tmpDir)
pancanTable <- loadClin(getFileLocation(syn))
nf <- nrow(pancanTable)

pancanTable['BLCA',"methylation_27k"] <- "syn1889358"
pancanTable['HNSC',"methylation_27k"] <- "syn1889356"


meta.pancan = list()

for(t in tags){

message("Processing ", t, "...")
temp=list()

#=== GE ===
synID <- pancanTable[t, "RNASeq"]
if(!is.na(synID)){
  syn <- loadEntity(synID)
  ge <- loadExpr(file.path(syn$cacheDir, syn$files[[1]]))
  nz <- apply(ge, 1, function(x){sum(x==0)})
  ge <- ge[nz < 0.5*(ncol(ge)), ]
	
  ge[ge==0] <- NA
  ge <- log2(ge)
  ge <- impute.knn(ge)$data
  ge <- normalizeBetweenArrays(ge)
  
  meta.rnaseq = createMetageneSpace(ge, attractome.rnaseq, rownamesMap=TRUE)$metaSpace
  colnames(meta.rnaseq) <- substr(colnames(meta.rnaseq), 1, 12)
  temp[['rnaseq']] = meta.rnaseq
}else{
  temp[['rnaseq']] = NULL
}

#=== miRNA ===
synID <- pancanTable[t, "miRNA"]
if(!is.na(synID)){
  syn <- loadEntity(synID)
  mirna <- loadExpr(file.path(syn$cacheDir, syn$files[[1]]))
  nz <- apply(mirna, 1, function(x){sum(x==0)})
  mirna <- mirna[nz < 0.5*(ncol(mirna)), ]
	
  mirna[mirna==0] <- NA
  mirna <- log2(mirna)
  mirna <- impute.knn(mirna)$data
  mirna <- normalizeBetweenArrays(mirna)
  mirna <- probeSummarization(ge=mirna, map=map.mirna, threshold=0.7, gene.colname="miR_stem")
  
  meta.mirna = createMetageneSpace(mirna, attractome.mirna, rownamesMap=TRUE)$metaSpace
  colnames(meta.mirna) <- substr(colnames(meta.mirna), 1, 12)
  temp[['mirna']] = meta.mirna
}else{
  temp[['mirna']] = NULL
}

#=== RPPA ===
synID <- pancanTable[t, "RPPA"]
if(!is.na(synID)){
  syn <- loadEntity(synID)
  rppa <- loadExpr(file.path(syn$cacheDir, syn$files[[1]]))
  meta.rppa = createMetageneSpace(rppa, attractome.rppa, rownamesMap=TRUE)$metaSpace
  colnames(meta.rppa) <- substr(colnames(meta.rppa), 1, 12)
  temp[['rppa']] = meta.rppa
}else{
  temp[['rppa']] = NULL
}

#=== meth ===
synID <- pancanTable[t, "methylation_27k"]
if(!is.na(synID)){
  syn <- synGet(synID, downloadFile=TRUE, downloadLocation=tmpDir)
  if(class(syn)=="File"){
    meth <- loadExpr(getFileLocation(syn))
  }else{
    syn <- loadEntity(synID)
    meth <- loadExpr(file.path(syn$cacheDir, syn$files[[1]]))
  }

  nz <- apply(meth, 1, function(x){sum(x==0)})
  meth <- meth[nz < 0.5*(ncol(meth)), ]
  meth <- impute.knn(meth)$data
  meta.meth = createMetageneSpace(meth, attractome.meth, map.meth, gene.colname="GeneSymbol")$metaSpace
  colnames(meta.meth) <- substr(colnames(meta.meth), 1, 12)
  temp[['meth']] = meta.meth
}else{
  temp[['meth']] = NULL
}

meta.pancan[[t]] = temp
}


