# TODO: I would like to be able to store rownames and column names in the widget, e.g., using hidden inputs, 
# but I think that some javascript is need to bind these to shiny inputs. 


renderEditableBooleanDataframe <- function(dataframe,label="df") {
	table.rows <- list()
	row.td <- lapply(c("",names(dataframe)), function(x) { tags$td(tags$strong(x)) })
	table.rows[[length(table.rows)+1]] <- row.td

	hidden.tags = list() 
	#for(r in 1:length(names(dataframe))) {
	#	hidden.tags[[length(hidden.tags)+1]] <- tags$input(type="hidden",id=paste0(label,".col.",r),names(dataframe)[r])
	#}
	#for(c in rownames(dataframe)) {
	#	hidden.tags[[length(hidden.tags)+1]] <- tags$input(type="hidden",id=paste0(label,".row.",c),rownames(dataframe)[c])
	#}

	for(x in 1:nrow(dataframe)) {
		row.td <- list()
		row.td[[length(row.td)+1]] <- tags$td(rownames(dataframe)[x])
		for (y in 1:ncol(dataframe)) {
			print(isolate({ paste0(label,".",x,".",y) }))
			row.td[[length(row.td)+1]] <- tags$td(checkboxInput(paste0(label,".",x,".",y), "", value=dataframe[x,y]))
		}
		table.rows[[length(table.rows)+1]] <- tags$tr(row.td)
	}
	list(hidden.tags,tags$table(border="1",table.rows))
}

booleanMatrixFromEditTable <- function(input,label,row.names,col.names) {
	maxsize <- 1000
	data <- c()
	max.x <- NULL
	max.y <- NULL
	for(x in 1:maxsize) {
		if (is.null(input[[paste0(label,".",x,".",1)]])) {
			max.x <- x-1
			break;
		}
		for(y in 1:maxsize) {
			if (is.null(input[[paste0(label,".",x,".",y)]])) {
				max.y <- y-1
				break
			} else {
				data[length(data)+1] <- c(input[[paste0(label,".",x,".",y)]])
			}
		}
	}
	if (is.null(max.y) || is.null(max.x)) {
		NULL
	} else {
		matrix(nrow=max.x, ncol=max.y, data=data,byrow=T)
		#row.names <- sapply(1:max.x, function(x) { input[[paste0(label,".row.",x)]] }) 
		#col.names <- sapply(1:max.y, function(x) { input[[paste0(label,".col.",x)]] }) 
	}
} 

