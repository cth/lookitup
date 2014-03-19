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

renderCohortsPanel <- function(input,session) {
	if (is.null(session$run.analysis)) {
		if (is.null(input$cohorts)) {
			wellPanel(checkboxGroupInput("cohorts", "Select cohorts:",choices=cohorts$name, selected=session$cohorts))
		} else {
			wellPanel(checkboxGroupInput("cohorts", "Select cohorts:",choices=cohorts$name, selected=input$cohorts))
		}
	} else {
		renderCohortSummary(session)
	}
}


renderSelectPhenotypes <- function(input,session,select) {
	controls <- list(selectInput("phenotypes.all.none", c("","All","None"), choices=c("","all","none")))
	session$stratificationProfiles <- session$stratificationProfiles

	for(p in names(phenotypes)) {
		strat.pheno <- paste0("stratification.",p) 
		controls[[length(controls)+1]] <-
			tags$table(cellspacing=5,
				tags$tr(tags$td(colspan=2, strong(p))), 
				tags$tr(
					if (p %in% select)
						tags$td(checkboxInput(p,"",value=T))
					else
						tags$td(checkboxInput(p,"",value=F))
					,
					tags$td(selectInput(paste0("transformation.", p),"Transformation", choices=c("none", "rank", "log"), selected="none")),
					if (!is.null(session$stratificationProfiles) && !identical(input[[strat.pheno]],names(session$stratificationProfiles))) {
						tags$td(selectInput(paste0("stratification.", p),"Stratification", choices=c("none", names(session$stratificationProfiles)), selected="none"))
					} else {
						tags$td(selectInput(paste0("stratification.", p),"Stratification", choices=c("none"), selected="none"))
					}
				))
	}

	return(controls)
}

renderPhenotypePanel <- function(input,session) {
	print("renderPhenotypePanel")
	if (is.null(session$run.analysis)) {
		# Case 1
		if (is.null(session$phenotypes) && is.null(session$phenotypes.none)) {
			print("case 1")
			session$phenotypes <- names(phenotypes) 
		} 
		# Case 2
		if (!is.null(input$phenotypes.all.none) && input$phenotypes.all.none == "all") {
			print("case 2")
			session$phenotypes <- names(phenotypes)
		}
		# Case 3
		if (!is.null(input$phenotypes.all.none) && input$phenotypes.all.none == "none") {
			print("case 3")
			# This is a hack to avoid returning to case 1 
			session$phenotypes <- c() 
			session$phenotypes.none <- T
		}
		if (!is.null(input$phenotypes.all.none) && input$phenotypes.all.none == "") {
			session$phenotypes <- c()
			session$phenotype.stratification <- list() 
			for(p in names(phenotypes)) {
				if(!is.null(input[[p]]) && input[[p]]==T) 
					session$phenotypes <- c(p,session$phenotypes)

				p.strat <- input[[paste0("stratification.", p)]]

				if(!is.null(p.strat) && p.strat != "none") 
					session$phenotype.stratification[[p]]  <- p.strat
			}

		}
		renderSelectPhenotypes(input,session,session$phenotypes)
	} else {
		renderPhenotypeSummary(session)
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
	if (is.null(session$run.analysis)) {
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
	} else {
		renderCovariateSummary(session)
	}
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

renderStratificationPanel <- function(input,session) {
	print("renderStratificationPanel")

	# If users has filled in the tab-name input field and pressed "Create", 
	# go ahead and create.. 
	if (!is.null(input$newStratTabButton) && input$newStratTabButton != 0 && !is.null(input$newStratTabName) && input$newStratTabName != "") {
		newStratTab <- tabPanel(input$newStratTabName, renderStratificationProfile(input,session,input$newStratTabName))
		
		if (is.null(session$stratificationProfiles)) 
			session$stratificationProfiles <- list()
		session$stratificationProfiles[[input$newStratTabName]] <- newStratTab
	}

	if (!is.null(session$stratificationProfiles) && input$newStratTabButton == 0 && !is.null(input$delStratTabButton) && input$delStratTabButton != 0 && !is.null(input$delStratTabSelect)) {
		tmp.stratificationProfiles <- session$stratificationProfiles
		if (length(tmp.stratificationProfiles) > 1) {
			session$stratificationProfiles <- list()
			keep.profiles <- setdiff(names(tmp.stratificationProfiles),input$delStratTabSelect)
			for(p in names(tmp.stratificationProfiles)) 
				if (p %in% keep.profiles)
					session$stratificationProfiles[[p]] <- tmp.stratificationProfiles[[p]]
		} else {
			session$stratificationProfiles <- NULL
		}
	}

	tabPanels <- list(
		tabPanel("Statification Profiles",
			sidebarPanel(
				textInput("newStratTabName", "Create stratification profile", ""),
				actionButton("newStratTabButton", "Create"),
				if (!is.null(session$stratificationProfiles))
					span(
						selectInput("delStratTabSelect","Remove stratification profile", choices=names(session$stratificationProfiles)),
						actionButton("delStratTabButton", "Remove"))
				else
					span("")
			)
		)
	)

	print("probably goes haywire here")

	if (!is.null(session$stratificationProfiles)) {
		for(tab in names(session$stratificationProfiles))
			if (is.null(session$stratificationProfiles[[tab]]))
				tabPanels[[tab]] <- tabPanel(input$newStratTabName, renderStratificationProfile(input,session,input$newStratTabName))
			else
				tabPanels[[tab]] <- session$stratificationProfiles[[tab]] 
	}


	do.call(tabsetPanel,tabPanels)
}

renderStratificationProfile <- function(input,session,name) {
	controls=list()
	print(paste0("renderStraticitionProfile:", name))

	for(p in names(phenotypes)) {
		if(!is.factor( phenotypes[,which(colnames(phenotypes)==p)] )) {
			values <- phenotypes[,which(colnames(phenotypes)==p)]
			uniq.values <- setdiff(unique(values),NA)

			# If the phenotype is represented by a few (upto 10) discrete values
			# then, make choice using checkbox, otherwise assume continuous range
			# and use a slider for input
			input.key <- paste("stratification",name,p,sep=".")
			if (length(uniq.values) <= 10) {
				controls[[length(controls)+1]] <- wellPanel(checkboxGroupInput(input.key, p, choices=uniq.values,selected=uniq.values))
			} else {
				minRange=min(uniq.values,na.rm=T)
				maxRange=max(uniq.values,na.rm=T)
				if (!is.null(session[[input.key]])) {
					print(paste("session ->",session[[input.key]]))
					select.choices <- session[[input.key]]
				} else {
					select.choices <- c(minRange,maxRange)
				}
				print(input.key)
				controls[[length(controls)+1]] <- wellPanel(sliderInput(input.key, p,min=minRange,max=maxRange,value=select.choices))
			}
		}
	}
	tags$span(controls)
}

saveStratificationProfiles <- function(input,session) {
	for (tab in names(session$stratificationProfiles))
		session$stratificationProfiles[[tab]] <- tabPanel(tab, renderStratificationProfile(input,session,tab))

	if (!is.null(session$stratificationProfiles)) {
		for (profile in names(session$stratificationProfiles)) {
			for(pheno in names(phenotypes)) {
				if(!is.factor( phenotypes[,which(colnames(phenotypes)==pheno)] )) {
					values <- phenotypes[,which(colnames(phenotypes)==pheno)]
					uniq.values <- setdiff(unique(values),NA)
					minRange=min(uniq.values,na.rm=T)
					maxRange=max(uniq.values,na.rm=T)
					# Determine if we should extract value from slider or from 
					input.key <- paste("stratification",profile,pheno,sep=".")
					if (!is.null(input[[input.key]])) {
						session[[input.key]] <- input[[input.key]]
						print(paste("Saving", input.key, input[[input.key]]))
					}
				}
			}
		}
	}
	if (!is.null(session$stratificationProfiles)) {
		for(tab in names(session$stratificationProfiles))
			session$stratificationProfiles[[tab]] <- NULL
	}

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
		saveStratificationProfiles(input,session)
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
		renderCohortSummary(session),
		renderPhenotypeSummary(session),
		renderCovariateSummary(session)
	)
}

renderInputPanel <- function(input,session) {
	if (is.null(session$run.analysis)) {
		if (!is.null(session$input.list) && is.null(input$input.list)) {
	    	output<-span(h6("List of genes, SNPs (rsnumbers) and genomic ranges to by analysed (one per line)"),
   	         	 tags$textarea(id="input.list", rows="20", cols="60", session$input.list)) 
		} else {
	    	output<-span(h6("List of genes, SNPs (rsnumbers) and genomic ranges to by analysed (one per line)"),
   	         	 tags$textarea(id="input.list", rows="20", cols="60", ""))
		}
		if (!is.null(input$input.list))
			session$input.list <- input$input.list
		print(paste("input.list:", session$input.list))
		output
	} else {
		renderInputSummary(session)
	}
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
renderCovariateSummary <- function(session) { 
	# We use list() as placeholder for "empty" input. This is a hack to sure only session is used
	covar.matrix <- createCovariateMatrix(list(),session)
	if (is.null(covar.matrix)) {
		wellPanel(h4("Covariates"), p("none selected"))
	} else {
		wellPanel(h4("Covariates"),renderBooleanDataframe(covar.matrix,"covar.summary","trait",F))
	}
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
	print("call rangeExploratorium")
    displayPanel <- list()
    rangeChosen <- session$range

    if(!is.null(session$lookupType) && session$lookupType != 'try-error'){
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
			print("hack 1")
            rangeChosen <- session$range
            rangeStart <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][1])
            rangeEnd  <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][2])
			if(rangeStart-1000>0){
            	rangeStartMinus1000 <- rangeStart-1000
        	}else{
            	rangeStartMinus1000 <- 0
        	}
        }

		# Hack two: If User changed both ends of slider, then he didn't and it was a bug
		if ( (sessionStart != rangeStart) && (sessionEnd != rangeEnd) ) {
			print("hack 2")
			print(rangeChosen)
			print(paste(sessionStart,rangeStart,sessionEnd,rangeEnd))
            rangeChosen <- session$range
			print(rangeChosen)
			rangeStart <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][1])
            rangeEnd  <- as.integer(strsplit((strsplit(rangeChosen,":")[[1]][2]),"-")[[1]][2])
			if(rangeStart-1000>0){
            	rangeStartMinus1000 <- rangeStart-1000
        	}else{
            	rangeStartMinus1000 <- 0
        	}


		}

		print(rangeStart)
		print(rangeEnd)
      

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
  
    if(!is.null(rangeChosen) && !is.na(rangeChosen) && session$range != rangeChosen){
		print("update session range")
        session$range <- rangeChosen
    }

	Sys.sleep(1)
   
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
	#sessionDefault(session,"phenotypes", names(phenotypes))
            
	# Render tabs   
	output$input.tab <- renderUI({ renderInputPanel(input,session) })
    output$exploratorium.tab <- renderUI({ renderExploratoriumPanel(input,session) })
	output$cohorts.tab <- renderUI({ renderCohortsPanel(input,session) })
	output$stratification.tab <- renderUI({ renderStratificationPanel(input,session) })
	output$phenotypes.tab <- renderUI({ renderPhenotypePanel(input,session) })
	output$covariates.tab <- renderUI({ renderCovariatesPanel(input,session) })
	output$analysis.tab <- renderUI({ renderAnalysisPanel(input,session) })

	# FIXME: Memoize loading to avoid double loading
	output$snp.table <- renderDataTable({
		result <- data.frame(a=c("default", "option")) 
		if (!is.null(session$select.analysis) && file.exists(result.file(session$session.key))) {
			load(result.file(session$session.key))
		}
		result$snp.table
	})

	output$assoc.table <- renderDataTable({
		result <- data.frame(a=c("default", "option")) 
		if (!is.null(session$select.analysis) && file.exists(result.file(session$session.key))) {
			load(result.file(session$session.key))
		}
		result$assoc.table
	})


	output$results.tab <- renderUI({
		print("Render result.tab:")
		tabPanels <- list(id="result.tabs") 
		tabPanels[[length(tabPanels)+1]] <-	tabPanel("SNP Table", dataTableOutput("snp.table"))
		tabPanels[[length(tabPanels)+1]] <-	tabPanel("Association Table", dataTableOutput("assoc.table"))
		do.call(tabsetPanel,tabPanels)
	})


	output$summary.tab <- renderUI({
		# update session with values from input sessionUpdate(session,"cohorts", input$cohorts)
		sessionUpdate(session,"phenotypes", input$phenotypes)
		renderSummaryPanel(input,session) })

	# Produce session.key and save session as side-effect
	output$session.key <- renderUI({
		if (!is.null(input$save.session) && input$save.session > 0 && !identical(input$save.session,session$save.session)) {

			# update session with values from input
			sessionUpdate(session,"cohorts", input$cohorts)
			sessionUpdate(session,"phenotypes", input$phenotypes)

			saveStratificationProfiles(input,session)
			
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
