library(shiny)
library(NCBI2R)
library(memoise)

source("boolEditTable.R")
source("session.R")

test.genes = list(
	"FTO" = list("chromosome" = 16, start=53737875, end=54148379,strand="+",description="blah blah blah")
)

# Read list of cohorts which is expected to be in tab separated format with two colummns: <cohort.id> <cohort.name>
cohorts <- read.table("data/cohorts.txt", head=F)
names(cohorts) <- c("id", "name")
cohorts$name <- paste(cohorts$name, "(", cohorts$id, ")")

# Read phenotype database
phenotypes <- read.table("data/MergedPhenotypes_13dec2013.txt",head=T)

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

##
# Step 5: Select covariates
renderCovariatesPanel <- function(input,session) {
	if (is.null(input$phenotypes)) {
		#wellPanel(h6("Please select phenotypes in the Phenotypes tab first"))
		use.phenotypes <- session$phenotypes
	} else {
		use.phenotypes <- input$phenotypes
	}
		# Add covar button pressed:
		#print(paste("input:",input$add.covar))
		#print(paste("session:", session$add.covar))
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

		span(
		sidebarPanel(
	  		actionButton("add.covar", "Add"),
			actionButton("remove.covar", "Remove"),
			selectInput("select.covar","Covariates:",choices=names(phenotypes))
		),
		span(
			if (length(use.phenotypes) > 0 && length(session$selected.covariates) > 0) {
				# Create dataframe representing selected phenotypes x selected covariates
				covar.matrix <- matrix(nrow=length(use.phenotypes),ncol=length(session$selected.covariates),rep(F, length(use.phenotypes) * length(session$selected.covariates)))
				print(session$selected.covariates)
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
				renderEditableBooleanDataframe(df.covar.matrix,"covar.matrix","trait")
			} else {
				span("")
			}
		)
		)
	#}
}

renderSummaryPanel <- function(input,session) {
	span(
		renderGeneSummary(input),
		renderCohortSummary(session),
		renderPhenotypeSummary(session),
		renderCovariateSummary(session)
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
	wellPanel(
		if (is.null(items)) {	
			span(h4(caption), p("None selected"))
		} else {
			span(h4(caption), p(paste(items,sep=",")))
		}
	)
}

renderCohortSummary <- function(session) { itemListSummary(session$cohorts, "Cohorts:") }
renderPhenotypeSummary <- function(session) { itemListSummary(session$phenotypes, "Phenotypes:") }
renderCovariateSummary <- function(session) { itemListSummary(session$covariates, "Covariates:") }

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
	output$summary.tab <- renderUI({ renderSummaryPanel(input,session) })

	# Produce session.key and save session as side-effect
	output$session.key <- renderUI({

		if (!is.null(input$save.session) && input$save.session > 0 && !identical(input$save.session,session$save.session)) {

			# update session with values from input
			sessionUpdate(session,"cohorts", input$cohorts)
			sessionUpdate(session,"phenotypes", input$phenotypes)

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
