library(shiny)
library(rjson)

source("env.R")
source("session.R")
source("boolEditTable.R")
source("lookup.R")

# Read list of cohorts which is expected to be in tab separated format with two colummns: <cohort.id> <cohort.name>
cohorts <- read.table(config$cohorts, head=F)
names(cohorts) <- c("id", "name")

# FIXME: Small hack to get studyid as part of cohort name
# This should be thrown out, but we do not have unique cohort names!
cohorts$name <- paste(cohorts$name, "(", cohorts$id, ")")

# Read phenotype database
phenotypes <- read.table(config$phenotypes,head=T)

# Restrict to cohorts for which we have phenotype information
cohorts <- cohorts[cohorts$id %in% unique(phenotypes$studyid),]


## 
# Step 3: Select cohorts to perfom analysis in
renderCohortsPanel <- function(input,session) {
	if (is.null(input$cohorts)) {
		wellPanel(checkboxGroupInput("cohorts", "Select cohorts:",choices=cohorts$name, selected=session$cohorts))
	} else {
		wellPanel(checkboxGroupInput("cohorts", "Select cohorts:",choices=cohorts$name, selected=input$cohorts))
	}
}


##
# Step 4: Select phenotypes relevant for analysis
renderSelectPhenotypes <- function(select,select.all.none="") {
	wellPanel(
		selectInput("phenotypes.all.none", c("","All","None"), choices=c("","all","none"), selected=select.all.none),
		checkboxGroupInput("phenotypes", "Select phenotypes:", 
			choices=names(phenotypes),
			selected=select)
	)
}

renderPhenotypePanel <- function(input,session) {
	if (!is.null(input$phenotypes.all.none) && input$phenotypes.all.none == "all") {
		renderSelectPhenotypes(names(phenotypes))
	} else if (is.null(input$phenotypes) && is.null(input$phenotypes.all.none)) {
		renderSelectPhenotypes(session$phenotypes)
	} else if (!is.null(input$phenotypes.all.none) && input$phenotypes.all.none == "none") {
		renderSelectPhenotypes(c())
	} else {
		renderSelectPhenotypes(input$phenotypes)
	}
}

createCovariateMatrix <- function(input,session) {
	if (is.null(input$phenotypes)) {
		use.phenotypes <- session$phenotypes
	} else {
		use.phenotypes <- input$phenotypes
	}

	if (!is.null(input$add.covar) && input$add.covar > 0 && !identical(input$add.covar,session$add.covar)) {
		session$add.covar <- input$add.covar
		if (!is.null(session$selected.covariates)) {
			session$selected.covariates <- sort(union(input$select.covar,session$selected.covariates))
		} else {
			session$selected.covariates <- input$select.covar
		}
	}

	if (!is.null(input$remove.covar) && input$remove.covar > 0 && !identical(input$remove.covar, session$remove.covar)) {
		session$remove.covar <- input$remove.covar
		if (!is.null(session$selected.covariates) && input$select.covar %in% session$selected.covariates) {
			session$selected.covariates <- setdiff(session$selected.covariates,input$select.covar)
		}
	}

	if (length(use.phenotypes) > 0 && length(session$selected.covariates) > 0) {
		# Create dataframe representing selected phenotypes x selected covariates
		covar.matrix <- matrix(nrow=length(use.phenotypes),ncol=length(session$selected.covariates),rep(F, length(use.phenotypes) * length(session$selected.covariates)))
		for(row in use.phenotypes) {
			for(col in session$selected.covariates) {
				list.key <- paste0("covar.matrix.",row,".",col)
				# If this covariate is currently selected in (input) table, then it should keep on being selected
				# On the other hand, if the input table have not been displayed yet, then we need to get values 
				# from session
				if (is.null(input[[list.key]]) && identical(session[[list.key]],T)) {
					covar.matrix[which(use.phenotypes==row),which(session$selected.covariates==col)] <- T
				} else if (identical(input[[list.key]],T)) {
					covar.matrix[which(use.phenotypes==row),which(session$selected.covariates==col)] <- T
					session[[list.key]] <- T 
				}
			}
		}

		df.covar.matrix <- data.frame(covar.matrix, row.names=use.phenotypes)
		colnames(df.covar.matrix) <- session$selected.covariates
		session$covariate.matrix <- df.covar.matrix
		return(session$covariate.matrix)
	} else {
		return(NULL)
	}
}

##
# Step 5: Select covariates
renderCovariatesPanel <- function(input,session) {
		df.covar.matrix <- createCovariateMatrix(input,session) 
		span(
			sidebarPanel(
	  			actionButton("add.covar", "Add"),
				actionButton("remove.covar", "Remove"),
				selectInput("select.covar","Covariates:",choices=names(phenotypes))
			),
			if (is.null(df.covar.matrix)) {
				p("No covariates selected")
			} else {
				renderBooleanDataframe(df.covar.matrix,"covar.matrix","trait")
			}
		)
}

renderAnalysisInput <- function(input,session) {
		displayElems = list()
#		for(worker in config$workers) {
##				displayElems[[length(displayElems)+1]] <- checkboxInput(worker$name, worker$name, F)
#		}
		displayElems[[length(displayElems)+1]] <- radioButtons("select.analysis", "Select analysis:", sapply(config$workers, function(x) { x$name }))
		displayElems[[length(displayElems)+1]] <- actionButton("run.analysis", "Run")

		return(displayElems)
}

renderAnalysisStatus <- function(input,session) {
	displayElems = list()

	displayElems[[length(displayElems)+1]] <- list(tags$strong(session$session.key))

	if (file.exists(result.file(session$session.key))) {
		displayElems[[length(displayElems)+1]] <- 
			list(tags$strong(session$select.analysis), tags$ul(tags$li("Finished")))
	} else { # Show "running" progress bar
		displayElems[[length(displayElems)+1]] <-
			list(tags$strong(session$select.analysis), tags$ul(tags$li(tags$progress(""))))
		displayElems[[length(displayElems)+1]] <- tags$script(paste0("delayedrefreshsession('", session$session.key, "');"))
	}

	return(displayElems)
}

renderAnalysisPanel <- function(input,session) {
	displayElems = list()
	if (!is.null(session$run.analysis)) {
		displayElems <- renderAnalysisStatus(input,session)
	} else if (!is.null(input$run.analysis) && input$run.analysis > 0) {
		session$session.key <- uniqueSessionKey()
		session$run.analysis <- input$run.analysis
		session$select.analysis <- input$select.analysis
		print(paste0("Saving session before analyses: ",session$session.key))
		save(session,file=session.file(session$session.key))
		write(c(), file=run.file(session$session.key))
		displayElems <- renderAnalysisStatus(input,session)
	} else {
		displayElems <- renderAnalysisInput(input,session)
	}
	wellPanel(displayElems)
}

renderSummaryPanel <- function(input,session) {
	span(
    	if (!is.null(session$run.analysis)) {
        	wellPanel(renderAnalysisStatus(input,session))
		} else {
			span("")
		},
		#renderInputSummary(session),
		renderCohortSummary(session),
		renderPhenotypeSummary(session),
		renderCovariateSummary(input,session)
	)
}

renderInputPanel <- function(input,session) {
    output<-span(h6("List of genes, SNPs (rsnumbers) and genomic ranges to by analysed (one per line)"),
            tags$textarea(id="input.list", rows="20", cols="60", ""))
	session$input.list <- input$input.list
	print(paste("input.list:", session$input.list))
	output
}



### Summary functions
renderGeneSummary <- function(geneInputted,gi) {
    		span(
			h6("Gene name: ",geneInputted),
			h6("Chromosome: ", gi$chr),
			h6("Range: ", gi$GeneLowPoint,"-",gi$GeneHighPoint),
			h6("Strand: ", gi$ori),
			h6("Description: ", gi$genesummary)
                    )
}

itemListSummary <- function(items,caption) {
	wellPanel(span(h4(caption),
		if (is.null(items)) {	
			p("None selected")
		} else {
			list(tags$ul(lapply(sort(items),function(i) { tags$li(i) } )))
		}
	))
}

renderInputSummary <- function(session) { 
	if (!is.null(session$input.list)) {
		print("rendering new summary")
		strsplit(session$input.list,"\n")[[1]]
		itemListSummary(, "Analysis regions:") 
	} else
		span()
}
renderCohortSummary <- function(session) { itemListSummary(session$cohorts, "Cohorts:") }
renderPhenotypeSummary <- function(session) { itemListSummary(session$phenotypes, "Phenotypes:") }
renderCovariateSummary <- function(input,session) { 
#	span(
#		itemListSummary(session$covariates, "Covariates:"),
		covar.matrix <- createCovariateMatrix(input,session)
		if (is.null(covar.matrix)) {
			wellPanel(h4("Covariates"), p("none selected"))
		} else {
			wellPanel(h4("Covariates"),renderBooleanDataframe(covar.matrix,"covar.summary","trait",F))
		}
#	)
	
}


renderExploratoriumPanel <- function(input,session){
    #Function that calls the right draw functions at the right time
    if(!is.null(session$lookupButton) && session$lookupButton!=input$lookupButton){
        #update session if it different from input (button is pressed). Update session input and button
        session$lookup <- input$lookup
        session$lookupButton <- input$lookupButton
        session$ExploratoriumSummaryPanel <- summaryExploratorium(input,session)
    }
    #Draw the session summary and a new search field if button have been pressed and updated
    if(!is.null(session$lookupButton) && session$lookupButton==input$lookupButton){
        session$ExploratoriumRangePanel <- rangeExploratorium(input,session)
        displayPanel <- list(drawExploratorium(input,session),session$ExploratoriumRangePanel,session$ExploratoriumSummaryPanel)
    }else{
    #Draw the search field without a summary and update session button
        session$lookupButton <- input$lookupButton
        displayPanel <- drawExploratorium(input,session)
    }
        
    return(displayPanel)
}


drawExploratorium <- function(input,session){
    #Function that draws the search field for lookup
    displayPanel <- list(
        mainPanel(
            textInput("lookup","Search:"),
            helpText(list("Submit a lookup in one of three formats:",tags$ul(tags$li("Range: chr1:12345-67890"),tags$li("SNP-name: rs1234567890"), tags$li("Gene-name: GENE"))))
            )
        )
    return(displayPanel)
}

rangeExploratorium <- function(input,session){
    displayPanel <- list()
    rangeChosen <- session$range

    if(session$lookupType != 'try-error'){
        #Original session start and end is used in hack
        sessionStart <- as.integer(strsplit((strsplit(session$range,":")[[1]][2]),"-")[[1]][1])
        sessionEnd <- as.integer(strsplit((strsplit(session$range,":")[[1]][2]),"-")[[1]][2])
        
        if(!is.null(input$sliderRange) && paste0('chr',session$chr,':',as.character(input$sliderRange[1]),'-',as.character(input$sliderRange[2]))!=rangeChosen){
            rangeChosen <- paste0('chr',session$chr,':',as.character(input$sliderRange[1]),'-',as.character(input$sliderRange[2]))
        }
        
        rangeStart <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][1])
        rangeEnd  <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][2])
        
        if(rangeStart-1000>0){
            rangeStartMinus1000 <- rangeStart-1000
        }else{
            rangeStartMinus1000 <- 0
        }

        #Hack to overwrite input if session$range starts or ends outside previous choseable area
        if( (sessionStart <= rangeStartMinus1000 || sessionStart > rangeEnd ) && (sessionEnd > rangeEnd+1000 || sessionEnd < rangeStart) )
        {
            rangeChosen <- session$range
            rangeStart <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][1])
            rangeEnd  <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][2])
            
            if(rangeStart-1000>0){
                rangeStartMinus1000 <- rangeStart-1000
            }else{
                rangeStartMinus1000 <- 0
            }
        }
      

        if(rangeStart>rangeEnd){
            session$type <- 'try-error'
            rangeChosen <- NA
            displayPanel[[length(displayPanel)+1]] <- list(
                h4('Impossible range (start>end)!')
                )
        }else{
            displayPanel[[length(displayPanel)+1]] <- list(
                sliderInput("sliderRange","Modify Range:",min=rangeStartMinus1000, max=rangeEnd+1000,value = c(rangeStart,rangeEnd))
                )
            displayPanel[[length(displayPanel)+1]] <- list(
                h4("Range chosen: ",tags$div(id="range.chosen", rangeChosen))
                )
        }
    }
  
    if(!is.na(rangeChosen) && session$range != rangeChosen){
        session$range <- rangeChosen
    }
   
    return(displayPanel)
}




summaryExploratorium <- function(input,session) {
    #Function that looks up session lookup input and report the summary and draws the modify range slider
    lookupTranslated <- lookupInterpreter(session$lookup)
    session$range <- lookupTranslated$range
    session$lookupType <- lookupTranslated$type
    session$chr <- lookupTranslated$chr
    displayPanel <- list()

    displayPanel[[length(displayPanel)+1]] <- list(
        h4("Summary:"),
        h6("Name: ",lookupTranslated$name),
        h6("Chromosome: ", lookupTranslated$chr),
        h6("Range: ", lookupTranslated$range),
        h6("Strand: ", lookupTranslated$strand),
        h6("Description: ", lookupTranslated$summary)
        )
    
    return(displayPanel)
}




# Define server logic
shinyServer(function(input, output, session) {
	# Load session if URL indicates that we should 
	# When the session.restore button is pushed, a JavaScript is invoked
	# to reload page with session identifier as part of query URL
	session <- isolate({
		queryString <- session$clientData$url_search
		query <- parseQueryString(queryString)
		if (!is.null(query$session) && !identical(query$session, session$session.key)) { 
			print(paste("load session:", query$session))
			load(session.file(query$session))
		}
		session
	})

	# Initialize session defaults:
	sessionDefault(session,"cohorts", cohorts$name)
	sessionDefault(session,"phenotypes", names(phenotypes))
            
	# Render tabs   
	output$input.tab <- renderUI({ renderInputPanel(input,session) })
    output$exploratorium.tab <- renderUI({ renderExploratoriumPanel(input,session) })
	output$cohorts.tab <- renderUI({ renderCohortsPanel(input,session) })
	output$phenotypes.tab <- renderUI({ renderPhenotypePanel(input,session) })
	output$covariates.tab <- renderUI({ renderCovariatesPanel(input,session) })
	output$analysis.tab <- renderUI({ renderAnalysisPanel(input,session) })

	output$result.table <- renderDataTable({
		result <- data.frame(a=c("default", "option")) 
		if (!is.null(session$select.analysis) && file.exists(result.file(session$session.key))) {
			load(result.file(session$session.key))
		}
		result
	})

	output$results.tab <- renderUI({
		print("Render result.tab:")
		tabPanels <- list(id="result.tabs") 
		tabPanels[[length(tabPanels)+1]] <-	tabPanel("Table", dataTableOutput("result.table"))
		do.call(tabsetPanel,tabPanels)
	})


	output$summary.tab <- renderUI({
		# update session with values from input sessionUpdate(session,"cohorts", input$cohorts)
		sessionUpdate(session,"phenotypes", input$phenotypes)
		renderSummaryPanel(input,session) })

	# Produce session.key and save session as side-effect
	output$session.key <- renderUI({

		if (!is.null(input$save.session) && input$save.session > 0 && !identical(input$save.session,session$save.session)) {

			# Make sure that any running analyses are not carried over to new session
			# update session with values from input
			sessionUpdate(session,"cohorts", input$cohorts)
			sessionUpdate(session,"phenotypes", input$phenotypes)
			
			# Make sure that any running analyses are not carried over to new session
			session$run.analysis <- NULL
			session$workers <- NULL

			print(session$cohorts)
			session$session.key <- uniqueSessionKey() 
			print(paste0("Saving session: ",session$session.key))
			save(session,file=session.file(session$session.key))
		}
		if (is.null(session$session.key))
			textInput("session.key","Session identifier:", "")
		else
			textInput("session.key","Session identifier:", session$session.key)
	})
}
)
