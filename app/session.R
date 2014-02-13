
sessionDefault <- function(session, session.key, value) {
	if (is.null(session[[session.key]])) {
		print(paste0("session.default",session.key))
		session[[session.key]] <- value
	}
}

sessionUpdate <- function(session, session.key, value) {
	if (!is.null(value)) {
		print(paste0("session.update: ", session.key))
		session[[session.key]] <- value
	}
}

uniqueSessionKey <- function() {
	repeat {
		key <- Reduce(function(x,y) { paste0(x,y) }, sample(seq(0,9),20,replace=T))
		file <- paste0(key,".session.Rdata")
		if (!file.exists(session.file(key)))
			break
	}
	print(key)
	return(key)
}


