library("rjson")

## Read config file
config <- fromJSON(file="config.json")

### Managing file names in the shared environment between server, dispatcher and workers
session.key <- function(sess.run.file) { 
	strsplit(basename(sess.run.file),".",fixed=T)[[1]][1]
}

session.file <- function(sess.key) {
	paste0(config$session.directory,sess.key,".session")
}


worker.script <- function(name) {
	for(w in config$workers) {
		if (w$name == name)
			return(w$script)
	} 
}

worker.extension <- function(name) {
	for(w in config$workers) {
		if (w$name == name)
			return(w$extension)
	} 
}

run.file <- function(sess.key) {
	paste0(config$session.directory, sess.key, ".run")
}

result.file <- function(sess.key, worker) {
	paste0(config$session.directory, sess.key, ".", worker.extension(worker))
} 

