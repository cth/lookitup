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
	rangesFile <- workerTempFile("plink-ranges-") 
	print("----")
#	system(paste("cat", rangesFile))
	GenomicRangesToPlinkRangeFile(regions,rangesFile)

	studyids <- allCohorts[selectedCohorts %in% allCohorts$name,1]
	print(paste0("studyids:", studyids))
	particids <- phenotypes[phenotypes$studyid %in% studyids, 1] 
	#print(paste0("particids:", particids))
	keep.list <- data.frame(pid=particids, fid=rep(1, length(particids)))
	keep.list.file <- workerTempFile("plink-keep-") 
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

extractSelectPhenotypes <- function(session.phenotypes,phenotypes) {
	print(session.phenotypes)
	select.phenotypes <- phenotypes[,which(colnames(phenotypes) %in% c("particid",session.phenotypes))]
	
	# NOTE: We asssume that the first column of phenotypes is an ID 
	# that uniquely identifies the individual. Additionally, plink 
	# needs a family id column 
	plink.pheno <- data.frame(cbind(select.phenotypes[,1], rep(1,nrow(select.phenotypes)), select.phenotypes[,(2:ncol(select.phenotypes))]))
	names(plink.pheno) <- c("FID", "IID", names(select.phenotypes)[2:ncol(select.phenotypes)])
	plink.pheno.file <- workerTempFile("plink.pheno") 
	write.table(plink.pheno,plink.pheno.file, quot=F, row.names=F, col.names=T) 

	# Test it
	print(head(read.table(plink.pheno.file, head=T)))

	return(plink.pheno.file)
}

plinkAssociationAnalysis <- function(plink.phenotypes.file,plink.covariates.file, covar, plinkstem) {
	assoc.tables <- list()

	plink.phenotypes <- read.table(plink.phenotypes.file,head=T,nrow=1)
	
	for (pheno in names(plink.phenotypes)[3:ncol(plink.phenotypes)]) {
		print(paste("Running PLINK analysis for phenotype", pheno))
		assoc.file <- workerTempFile("plink.analysis")

                if(!is.null(plink.covariates.file)){             
                    print(covar)
                    tmp <- as.vector(as.matrix((covar[which(rownames(covar)==pheno),])))
                    print(tmp)
                    covar.names <- names(covar)[tmp]
                }
                else{
                    covar.names <- c()
                }
		if (length(covar.names) > 0) {
			plink.cmd <- paste("plink --noweb --bfile", plinkstem, "--pheno", plink.phenotypes.file, "--pheno-name", pheno, "--covar", plink.covariates.file, "--covar-name", paste(covar.names,sep=","), "--linear --out", assoc.file)
		} else {
			plink.cmd <- paste("plink --noweb --bfile", plinkstem, "--pheno", plink.phenotypes.file, "--pheno-name", pheno, "--linear --out", assoc.file)
		}
		print(plink.cmd)
		system(plink.cmd)
		t <- read.table(paste0(assoc.file,".assoc.linear"), head=T,stringsAsFactors=F)
		if (nrow(t) > 0) {
			t$phenotype <- rep(pheno, nrow(t))
			assoc.tables[[length(assoc.tables)+1]] <- t
		}
	}
	assoc.table <- Reduce(rbind,assoc.tables)
	print(assoc.table)
	assoc.table <- assoc.table[assoc.table[["TEST"]] == "ADD",]

	print(dim(assoc.table))
	print(assoc.table)

	return(assoc.table)
}

worker({
	# Get worker configuration
	worker.config <- NULL
	for(w in config$workers) {
		if (w$name == session$select.analysis) 
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
        
	plink.pheno <- extractSelectPhenotypes(session$phenotypes,phenotypes)
        if(!is.null(session$covariate.matrix)){
            plink.covar <- extractSelectPhenotypes(names(session$covariate.matrix),phenotypes)
        }
        else{
            plink.covar <- NULL
        }
            
	ressult <- list(
		snp.table = merge(frq,hwe),
		assoc.table =  plinkAssociationAnalysis(plink.pheno,plink.covar,session$covariate.matrix,stem)
	)	
})
