# Common stuff for all types of workers


options(echo=TRUE) # if you want see commands in output file

worker <- function(worker.code) {
	args <- commandArgs(trailingOnly = TRUE)
	print(args)
	# Load session:
	load(args[1]) 
	input.file <- args[1]
	output.file <- args[2]
	result <- do.call(worker.code)
	save(result, output.file)
}
