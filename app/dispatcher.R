#
# The dispatcher runs in its own process. It is responsible for starting analysis jobs.

source("env.R")

launch.worker <- function(worker,session) {
	cmdline=paste('Rscript ', worker.script(worker),session.file(session), result.file(session, worker), "&" )
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
