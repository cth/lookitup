library(shiny)
library(NCBI2R)
library(memoise)
library(rjson)

source("boolEditTable.R")
source("session.R")

test.genes = list(
	"FTO" = list("chromosome" = 16, start=53737875, end=54148379,strand="+",description="blah blah blah")
)

## Read config file
config <- fromJSON(file="config.json")

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
# Step 1: Select a gene or a genomic range to be analyzed
renderGeneSelect <- function(set.gene="", set.chromosome="",set.start=1, set.end=1000000) {
  sidebarPanel(
	  textInput("gene", "Gene:"),
	  actionButton("gene.lookup", "Lookup gene"),
      selectInput("chromosome", "Chromosome:", 
                choices = c(as.character(seq(1,22)),"X","Y","MT"),set.chromosome),
	  sliderInput("range", "Range:",
                min = set.start, max = set.end, value = c(set.start,set.end))
	)
}

##
# Step 2: Adjust selection after choice of gene
renderGeneAdjust <- function(set.gene="", set.chromosome,set.start, set.end, strand, description, input) {
  gene.size <-  set.end-set.start
  min.slider <- ifelse(set.start-gene.size > 0, set.start-gene.size, 0)  

  span(
  sidebarPanel(
	  textInput("gene.selected", "Gene:", set.gene),
#	  actionButton(input.gene.lookup, "Lookup gene"),      
      selectInput("chromosome", "Chromosome:", 
                choices = set.chromosome, set.chromosome),
	  sliderInput("range", "Range (gene):",
                min = min.slider, max = set.end+gene.size, value = c(set.start,set.end))
#	  sliderInput("range", "Range (analysis):",
#                min = min.slider, max = set.end+gene.size, value = c(set.start,set.end))
  )
      ,
        renderGeneSummary(input)
      )
}

renderGeneLookup <- function(input) {
    if(!is.null(input$gene.lookup) && input$gene.lookup > 0){
        gi <- isolate({gene.info.persist(input$gene)})
        if (class(gi)=="try-error") {
            renderGeneSelect()
        }
        else {
            renderGeneAdjust(paste(input$gene),gi$chr,gi$GeneLowPoint,gi$GeneHighPoint,gi$ori,gi$genesummary,input)
        }
        
    }
    else{
        span('')
    }
}



    
renderGenePanel <- function(input,session) {
    span(
        if (!is.null(input$gene.lookup) && input$gene.lookup != 0) {
            gi <- isolate({gene.info.persist(input$gene)})
            if (class(gi)=="try-error") {
                renderGeneSelect(input$gene)
            }
            else {
		renderGeneAdjust(paste(input$gene),gi$chr,gi$GeneLowPoint,gi$GeneHighPoint,gi$ori,gi$genesummary)
            }
	}

        else if (is.null(input$gene)) {
            renderGeneSelect()
	}
        else if (input$gene %in% labels(test.genes)) {
		renderGeneAdjust(paste(input$gene),test.genes[[input$gene]]$chromosome,test.genes[[input$gene]]$start,test.genes[[input$gene]]$end,test.genes[[input$gene]]$strand,test.genes[[input$gene]]$description)
	}
        else {
		renderGeneSelect(input$gene)
	}
        ,
        renderGeneSummary(input)
        )
}

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

renderAnalysisPanel <- function(input,session) {
	workerList = list()
	print(config$workers$name)
	workers <- c()
	if ((!is.null(input$run.analysis) && input$run.analysis > 0) || !is.null(session$run.analysis)) {
		
		session$session.key <- uniqueSessionKey()
		workerList[[length(workerList)+1]] <- tags$span(tags$strong("Session ID:"), tags$p(session$session.key))

		session$workers <- list()
		for(worker in config$workers) {
				session$workers[[worker$name]] <- input[[worker$name]]

				if (!is.null(input[[worker$name]]) && identical(input[[worker$name]],T))
					workerList[[length(workerList)+1]] <- list(tags$strong(worker$name), tags$progress(""))	
				else 
					workerList[[length(workerList)+1]] <- list(tags$strong(worker$name), tags$ul(tags$li("Skipped.")))
		}

		session$run.analysis <- input$run.analysis
		print(paste0("Saving session before analyses: ",session$session.key))
		file <- paste0(session$session.key, ".session.Rdata")
		save(session,file=file)
		file.run <- paste0(session$session.key,".run")
		write(c(), file=file.run)

	} else {
		for(worker in config$workers) {
				workerList[[length(workerList)+1]] <- checkboxInput(worker$name, worker$name, F)
		}
		workerList[[length(workerList)+1]] <- actionButton("run.analysis", "Run")

	}
	#checkboxGroup("workers", choices=workers, selected  
	wellPanel(workerList)
}

renderSummaryPanel <- function(input,session) {
	span(
		renderGeneSummary(input),
		renderCohortSummary(session),
		renderPhenotypeSummary(session),
		renderCovariateSummary(input,session)
	)
}

##
# Step 6: Select stratification
# TODO


### Summary functions
renderGeneSummary <- function(input) {
	if(!is.null(input$gene.lookup) && input$gene.lookup > 0) {
		span(
			h4("Gene:"),
			h6("Gene name: ",ifelse(is.null(input$gene),"Unknown", input$gene)),
			h6("Chromosome: ", input$chromosome),
			h6("Range: ", input$range[[1]],"-",input$range[[2]]),
			h6("Strand: ", isolate({gene.strand(input)})),
			h6("Description: ", isolate({gene.description(input)})))
	} else { span("") }
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
			load(paste0(query$session,".session.Rdata"))
			print(session$cohorts)
		}
		session
	})

	# Initialize session defaults:
	sessionDefault(session,"cohorts", cohorts$name)
	sessionDefault(session,"phenotypes", names(phenotypes))
	
	# Render tabs
	output$gene.tab <- reactive({
    	if(is.null(input$gene.lookup)){
        	renderGeneSelect(input,session)
    	} else if(input$gene.lookup == 0){
        	renderGeneSelect(input,session)
    	} else {
        	renderGeneLookup(input,session)
		}
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
			file <- paste0(session$session.key, ".session.Rdata")
			save(session,file=file)
		}
		if (is.null(session$session.key))
			textInput("session.key","Session identifier:", "")
		else
			textInput("session.key","Session identifier:", session$session.key)
	})
}
)
