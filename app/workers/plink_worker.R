source("worker.common.R")
source("lookup.R")

GenomicRangesToPlinkRangeFile <- function(ranges,file) {
	df <- as.data.frame(ranges) 
	fields <- sapply(c("seqnames", "start", "end", "name"),function(x) { which(names(df)==x) })
	print(fields)
	write.table(df[,fields],file=file, quot=F, row.names=F, col.names=F,sep="\t")
	system(paste("cat", file))
}

plinkExtractRegions <- function(regions, plinkStemIn) {
	# FIXME: use random/unique filename
	rangesFile <- workerTempFile("plink-ranges") 
	GenomicRangesToPlinkRangeFile(regions,rangesFile)

	plinkStemOut <- workerTempFile("plink-stem-")

	cmd <- paste("plink --noweb --bfile", plinkStemIn, "--range --extract",  rangesFile,  "--make-bed --out", plinkStemOut)
	print(cmd)
	system(cmd)

	return(plinkStemOut)
}

plinkFrequencies <- function(plinkStemIn) {
	plink.frq <- workerTempFile("plink.frq.")
	cmd <- paste("plink --noweb --bfile", plinkStemIn, "--freq --out", plink.frq)
	print(cmd)
	system(cmd)
	return(paste0(plink.frq,".frq"))
}

worker({
	# Get worker configuration
	worker.config <- NULL
	for(w in config$workers) {
		if (w$name == "plink analysis") 
			worker.config <- w
	}

	print(worker.config)
	print(paste("SESSION: ", session$session.key))

	# Get input regions from session
	lookupIdentifiers <- strsplit(session$input.list,"\n")

	ranges <- rangesFromInputList(lookupIdentifiers)

	print("RANGES:")
	print(ranges)
	print(as.data.frame(ranges))

	stem <- plinkExtractRegions(ranges,worker.config$plinkstem)  
	plink.frq <- plinkFrequencies(stem)

	read.table(plink.frq, head=T)
})
