library(shiny)
library(NCBI2R)
library(memoise)

source("env.R")
source("session.R")
source("boolEditTable.R")

test.genes = list(
	"FTO" = list("chromosome" = 16, start=53737875, end=54148379,strand="+",description="blah blah blah")
)

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
		for(worker in config$workers) {
				displayElems[[length(displayElems)+1]] <- checkboxInput(worker$name, worker$name, F)
		}
		displayElems[[length(displayElems)+1]] <- actionButton("run.analysis", "Run")

		return(displayElems)
}

renderAnalysisStatus <- function(input,session) {
	displayElems = list()
	displayElems[[length(displayElems)+1]] <- tags$span(tags$strong("Session ID:"), tags$p(session$session.key))

	for(worker in config$workers) {
#j		if ((!is.null(input[[worker$name]]) && identical(input[[worker$name]],T)) {
# ||
		if (session$workers[[worker$name]]) {
			print(paste("result: ",result.file(session$session.key,worker$name)))
			if (file.exists(result.file(session$session.key,worker$name))) {
				displayElems[[length(displayElems)+1]] <- 
					list(tags$strong(worker$name), tags$ul(tags$li("Finished")))
			} else { # Show "running" progress bar
				displayElems[[length(displayElems)+1]] <- 
					list(tags$strong(worker$name), tags$ul(tags$li(tags$progress(""))))
			}
		} else {
			displayElems[[length(displayElems)+1]] <- list(tags$strong(worker$name), tags$ul(tags$li("Skipped.")))
		}
	}
	
	return(displayElems)
}

renderAnalysisPanel <- function(input,session) {
	displayElems = list()
	print(config$workers$name)
	if (!is.null(session$run.analysis)) {
		displayElems <- renderAnalysisStatus(input,session)
	} else if (!is.null(input$run.analysis) && input$run.analysis > 0) {
		session$session.key <- uniqueSessionKey()
		session$run.analysis <- input$run.analysis
		session$workers <- list()
		for(worker in config$workers) {
			session$workers[[worker$name]] <- input[[worker$name]]
		}
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
#		renderGeneSummary(input),
		renderCohortSummary(session),
		renderPhenotypeSummary(session),
		renderCovariateSummary(input,session)
	)
}

##
# Step 6: Select stratification
# TODO


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

gene.info <- function(gene) {
	if (is.null(gene)) {
		return(try(a.little.bit.harder))
	}	
	cat(paste("gene.info: ", gene, "\n"))
	lookup.str <- paste0(gene, "[sym]")
	cat(paste("lookup str: ", lookup.str, "\n"))
	id <- try(GetIDs(lookup.str))
	if (class(id) == "try-error") {
		return(id)
	} else {
		cat(paste("got id:", id[1], "\n"))
		return(try(GetGeneInfo(id[1])))
	}
}

gene.info.persist <- memoise(gene.info)


gene.strand <- function(input) {
		if (input$gene %in% labels(test.genes)) {
			test.genes[[input$gene]]$strand
		} else if (!is.null(input$gene.lookup) && input$gene.lookup != 0) {
			gi <- gene.info.persist(input$gene)
			if (class(gi)=="try-error") {
				"NA"
			} else {
				gi$ori
			}
		}
}

gene.description <- function(input) {
		if (input$gene %in% labels(test.genes)) {
			test.genes[[input$gene]]$description
		} else if (!is.null(input$gene.lookup) && input$gene.lookup != 0) {
			gi <- gene.info.persist(input$gene)
			if (class(gi)=="try-error") {
				"NA"
			} else {
				gi$genesummary
			}
		}
}

#Lookup logical interpreter
lookupInterpreter <- function(inputLookup){
    returnOutput <- data.frame(NA,NA,NA,NA,NA)
    colnames(returnOutput) <- c('range','summary','name','strand','chr')
       #range as chr1-22,X,Y:123-456
    if(grepl(paste0(Reduce(function(...) {paste(...,sep="|") },paste0('chr',1:22)),'chrX|chrY'),inputLookup)){       
        returnOutput$range <- inputLookup
        returnOutput$name <- paste0("Range: ",inputLookup)
    } else if(grepl('rs[0-9]',inputLookup)){ #snp with rs name rs[0-9]
        returnOutput$range <- "rsnumber looked up and return as range"
        returnOutput$summary <- "summary of SNP"
        returnOutput$name <- inputLookup
        returnOutput$strand <- "strand?!"
        returnOutput$chr <- 'lolwut'
    } else {
            gi <- isolate({gene.info.persist(inputLookup)})
            if (class(gi)=="try-error") {
                print("Input not understanded, gene is not found or servers are down")
                returnOutput$name <- inputLookup
                returnOutput$summary <- "Input not understanded, gene is not found or servers are down"
            }
            else{
                returnOutput$range <- paste0(gi$chr,gi$GeneLowPoint,'-',gi$GeneHighPoint)
                returnOutput$name <- inputLookup
                returnOutput$strand <- gi$ori
                returnOutput$summary <- gi$genesummary
                returnOutput$chr <- gi$chr
            }
        }
    return(returnOutput)
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
			print(session$cohorts)
		}
		session
	})

	# Initialize session defaults:
	sessionDefault(session,"cohorts", cohorts$name)
	sessionDefault(session,"phenotypes", names(phenotypes))

       

        
            
	# Render tabs   
        # Render gene tab - fundemental different design by Vincent

      
        output$lookupPrint <- renderUI({
            lookupInputted <- input$lookup
            lookupTranslated <- lookupInterpreter(lookupInputted)
            print(lookupInputted)
            
            span(
                h6("Name: ",lookupTranslated$name),
                h6("Chromosome: ", lookupTranslated$chr),
                h6("Range: ", lookupTranslated$range),
                h6("Strand: ", lookupTranslated$strand),
                h6("Description: ", lookupTranslated$summary)
                )


            
#            print(lookupInputted)
#            gi <- isolate({gene.info.persist(lookupInputted)})
#            if (class(gi)=="try-error") {
#                print("Gene not found or servers are down")
#            }
#            else {
#                renderGeneSummary(lookupInputted,gi)
#            }
        })
        
	output$cohorts.tab <- renderUI({ renderCohortsPanel(input,session) })
	output$phenotypes.tab <- renderUI({ renderPhenotypePanel(input,session) })
	output$covariates.tab <- renderUI({ renderCovariatesPanel(input,session) })
	output$analysis.tab <- renderUI({ renderAnalysisPanel(input,session) })


	output$summary.tab <- renderUI({
		# update session with values from input
		sessionUpdate(session,"cohorts", input$cohorts)
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
