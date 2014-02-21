source("worker.common.R")
source("lookup.R")

GenomicRangesToPlinkRangeFile <- function(ranges,file) {
	df <- as.data.frame(ranges) 
	fields <- sapply(c("seqnames", "start", "end", "name"),function(x) { which(names(df)==x) })
	print(fields)
	write.table(df[,fields],file=file, quot=F, row.names=F, col.names=F,sep="\t")
	system(paste("cat", file))
}

# Extract selected input regions and individuals corresponding to selected cohorts from plink file
plinkExtract <- function(regions, allCohorts, selectedCohorts, phenotypes, plinkStemIn) {
	# FIXME: use random/unique filename
	rangesFile <- workerTempFile("plink-ranges") 
	print("----")
	system(paste("cat", rangesFile))
	GenomicRangesToPlinkRangeFile(regions,rangesFile)

	studyids <- allCohorts[selectedCohorts %in% allCohorts$name,1]
	print(paste0("studyids:", studyids))
	particids <- phenotypes[phenotypes$studyid %in% studyids, 1] 
	#print(paste0("particids:", particids))
	keep.list <- data.frame(pid=particids, fid=rep(1, length(particids)))
	keep.list.file <- workerTempFile("plink.keep") 
	write.table(keep.list,keep.list.file,quot=F,row.names=F, col.names=F)


	plinkStemOut <- workerTempFile("plink-stem-")

	cmd <- paste("plink --noweb --bfile", plinkStemIn, "--range --extract",  rangesFile, "--keep", keep.list.file, "--make-bed --out", plinkStemOut)
	print(cmd)
	system(cmd)

	return(plinkStemOut)
}

plinkFrequencies <- function(plinkStemIn) {
	plink.frq <- workerTempFile("plink.frq.")
	cmd <- paste("plink --noweb --bfile", plinkStemIn, "--freq --out", plink.frq)
	print(cmd)
	system(cmd)
	frq <- read.table(paste0(plink.frq,".frq"),head=T)
	frq[,-which(names(frq) %in% c("NCHROBS"))]
}

plinkHardyWeinberg <- function(plinkStemIn) {
	plink.hwe <- workerTempFile("plink.hwe.")
	cmd <- paste("plink --noweb --bfile", plinkStemIn, "--hardy --out", plink.hwe)
	print(cmd)
	system(cmd)
	hwe <- read.table(paste0(plink.hwe,".hwe"),head=T)
	hwe <- hwe[hwe$TEST == "ALL", which(names(hwe) %in% c("P", "GENO", "CHR", "SNP"))]
	names(hwe)[which(names(hwe)=="P")] <- "P.HWE"
	return(hwe)
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

	# Read list of cohorts which is expected to be in tab separated format with two colummns: <cohort.id> <cohort.name>
	cohorts <- read.table(config$cohorts, head=F)
	names(cohorts) <- c("id", "name")

	# FIXME: Small hack to get studyid as part of cohort name
	# This should be thrown out, but we do not have unique cohort names!
	cohorts$name <- paste(cohorts$name, "(", cohorts$id, ")")

	# Read phenotype database
	phenotypes <- read.table(config$phenotypes,head=T)

	# Restrict to cohorts for which we have phenotype information
	cohorts <- cohorts[cohorts$id %in% unique(phenotypes$studyid),]


	# Get input regions from session
	lookupIdentifiers <- strsplit(session$input.list,"\n")

	ranges <- rangesFromInputList(lookupIdentifiers)

	stem <- plinkExtract(ranges,cohorts,session$cohorts,phenotypes,worker.config$plinkstem)  
	frq <- plinkFrequencies(stem)
	hwe <- plinkHardyWeinberg(stem) 

	print(merge(hwe,frq))

	result <- merge(frq,hwe)
})
