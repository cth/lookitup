source("worker.common.R")

# FIXME: 
temporaryPlinkStem <- function() {
	"/tmp/plink"
}

plinkExtractRegion(region, plinkfile) {
	splitRegion <- strsplit(region, ":") 	
	chr <- substring(splitRegion[[1]],4,nchar(splitRegion[[1]]))

	# FIXME: This is human specific and should be configurable
	if (chr == "X") 
		chr <- 22
	else if (chr == "Y") 
		chr <- 23

	splitRange <- strsplit(splitRegion[[2]], "-")

	newPlink <- temporaryPlinkStem()

	cmd <- paste("plink --noweb", "--chr", chr, "--from-kb", splitRange[[1]] / 1000, "--to-kb", splitRange[[2]], "--make-bed --out", newPlink))
	system(cmd)
}

worker({
	# Get worker configuration
	worker.config <- NULL
	for(w in config$workers) {
		if (w$name == "plink analysis") 
			worker.config <- w
	}

	print(worker.config)

	print(session$session.key)
	Sys.sleep(10)
	read.table("workers/test.annotation.txt",head=T)
	#data.frame(a=c(1,2,3),b=c("a", "b", "c"))
})
