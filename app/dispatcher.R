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

output.file <- function(sess.key, worker) {
	paste0("session/", sess.key, ".", worker.extension(worker))
} 

rscript.qsub <- function(script,name) {
        qsub.header <- c( "#$ -S /bin/sh", paste("#$ -N", name), "#$ -cwd")
        write(c(qsub.header, read.delim(file="test.R", header=F, sep="\n",stringsAsFactors=F)$V1), file="tmp.sge", sep="\n")
        system("qsub tmp.sge")
}

launch.worker <- function(worker,session) {
	cmdline=paste('Rscript ', worker.script(worker),session.file(session), output.file(session, worker), "&" )
	if (identical(config$dispatcher$worker.type, "local")) {
		system(cmdline,intern=F)
	} else if (identical(config$dispatcher$worker.type, "grid")) {
        qsub <- c( "#$ -S /bin/sh", paste("#$ -N", session), "#$ -cwd", cmdline)
        write(qsub, file="tmp.sge", sep="\n")
        system("qsub tmp.sge && rm -f tmp.sge")
	} else { 
		print(paste("Unknown worker type: ", type))
	}
}

# Main loop: Look for sessions that are ready to be run and assign them to a worker
repeat {
	for(run in Sys.glob("session/*run")) {
		file.remove(run)
		load(session.file(session.key(run)))
		for(w in names(session$workers)) {
			if (session$workers[[w]] == T) {
				launch.worker(w, session.key(run))
			}
		}
	} 
	Sys.sleep(1)
	cat(".")
}
