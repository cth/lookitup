
# TODO: I would like to be able to store rownames and column names in the widget, e.g., using hidden inputs, 
# but I think that some javascript is need to bind these to shiny inputs. 

# This function renders a table with of check
renderBooleanDataframe <- function(dataframe,label="df",row.name.name="",editable=T) {
	table.rows <- list()
	row.td <- lapply(c(row.name.name,names(dataframe)), function(x) { tags$td(tags$strong(x)) })
	table.rows[[length(table.rows)+1]] <- row.td

	for(x in 1:nrow(dataframe)) {
		row.td <- list()
		row <- rownames(dataframe)[x]
		row.td[[length(row.td)+1]] <- tags$td(row)
		for (y in 1:ncol(dataframe)) {
			col <- names(dataframe)[y]
			if (editable) {
				row.td[[length(row.td)+1]] <- tags$td(checkboxInput(paste0(label,".",row,".",col), "", value=dataframe[x,y]))
			} else {
				row.td[[length(row.td)+1]] <- tags$td(dataframe[x,y])
			}
		}
		table.rows[[length(table.rows)+1]] <- tags$tr(row.td)
	}
	list(tags$table(border="1",table.rows))
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

