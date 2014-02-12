# Common stuff for all types of workers

worker <- function(worker.block) {
	args <- commandArgs(trailingOnly = TRUE)
	# Load session:
	load(args[1],.GlobalEnv) 
	result <- worker.block
	save(result, file=args[2])
}
