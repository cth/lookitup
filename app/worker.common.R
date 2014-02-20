# Common stuff for all types of workers
library("rjson")

## Read config file
config <- fromJSON(file="config.json")

# FIXME: I eventually want this to be configurable per worker
workerTempDir <- function() {
	tempdir()
} 

workerTempFile <- function(...) {
	tmpdir=workerTempDir()
	if (!file.exists(tmpdir)) { 
		system(paste("mkdir -p", tmpdir))
		print("MKDIR")
		print(tmpdir)
	}
	tempfile(..., tmpdir=tmpdir)
}

worker <- function(worker.block) {
	args <- commandArgs(trailingOnly = TRUE)
	# Load session:
	load(args[1],.GlobalEnv) 
	result <- worker.block
	#result.json <- toJSON(result)
	#cat(result.json,file=args[2])
	save(result, file=args[2])
}
