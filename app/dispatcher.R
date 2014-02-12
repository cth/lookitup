#
# The dispatcher runs in its own process. It is responsible for starting analysis jobs.

library("rjson")

## Read config file
config <- fromJSON(file="config.json")

session.key <- function(sess.run.file) { 
	strsplit(basename(sess.run.file),".",fixed=T)[[1]][1]
}

session.file <- function(sess.key) {
	paste0("session/",sess.key,".session.Rdata")
}


# Main loop: Look for sessions that are ready to be run and assign them to a worker
#repeat {
	ready.sessions <- sapply(Sys.glob("session/*run"),session.key)
	for(sk in ready.sessions) {
		print(sk)
		load(sessfile(sk))
	} 
#}
