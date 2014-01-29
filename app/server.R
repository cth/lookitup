library(shiny)
library(NCBI2R)
library(memoise)

test.genes = list(
	"FTO" = list("chromosome" = 16, start=53737875, end=54148379,strand="+",description="blah blah blah")
)

# Read list of cohorts which is expected to be in tab separated format with two colummns: <cohort.id> <cohort.name>
cohorts <- read.table("cohorts.txt", head=F)
names(cohorts) <- c("id", "name")
cohorts$name <- paste(cohorts$name, "(", cohorts$id, ")")

# Read phenotype database
phenotypes <- read.table("MergedPhenotypes_13dec2013.txt",head=T)

# Restrict to cohorts for which we have phenotype information
cohorts <- cohorts[cohorts$id %in% unique(phenotypes$studyid),]

##
# Step 1: Select a gene or a genomic range to be analyzed
renderGeneSelect <- function(set.gene="", set.chromosome="",set.start=1, set.end=1000000) {
sidebarPanel(
  sidebarPanel(
	  textInput("gene", "Gene:", set.gene),
	  actionButton("gene.lookup", "Lookup gene"),
      selectInput("chromosome", "Chromosome:", 
                choices = c(as.character(seq(1,22)),"X","Y","MT"),set.chromosome),
	  sliderInput("range", "Range:",
                min = set.start, max = set.end, value = c(set.start,set.end)),
	actionButton("gene.select.continue", "Continue")
  )
)
}

##
# Step 2: Adjust selection after choice of gene
renderGeneAdjust <- function(set.gene="", set.chromosome,set.start, set.end, strand, description) {
  gene.size <-  set.end-set.start
  min.slider <- ifelse(set.start-gene.size > 0, set.start-gene.size, 0)  

  sidebarPanel(
	  textInput("gene.selected", "Gene:", set.gene),
      selectInput("chromosome", "Chromosome:", 
                choices = set.chromosome, set.chromosome),
	  sliderInput("range", "Range:",
                min = min.slider, max = set.end+gene.size, value = c(set.start,set.end)),

	actionButton("gene.select.continue", "Continue")
  )
}


## 
# Step 3: Select cohorts to perfom analysis in
renderSelectCohorts <- function(cohorts,cohorts.selected,gene) {
	sidebarPanel(
		checkboxGroupInput("cohorts", "Select cohorts:", 
			choices=cohorts,
			selected=cohorts.selected),
		actionButton("select.cohorts.continue", "Continue")
	)
}

##
# Step 4: Select phenotypes relevant for analysis
renderSelectPhenotypes <- function(select,select.all.none="") {
	sidebarPanel(
		actionButton("select.phenotypes.continue", "Continue"),
		selectInput("phenotypes.all.none", c("","All","None"), choices=c("","all","none"), selected=select.all.none),
		checkboxGroupInput("phenotypes", "Select phenotypes:", 
			choices=names(phenotypes),
			selected=select)
	)
}


##
# Step 5: Select covariates
renderSelectCovariates <- function(set) {
	sidebarPanel(
		actionButton("select.covariates.skip", "Skip (use predefined covariates)"),
		actionButton("select.covariates.continue", "Continue"),
		checkboxGroupInput("covariates", "Select covariates:", 
			choices=names(phenotypes),
			selected=c())
	)
}

##
# Step 6: Select stratification
# TODO

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
shinyServer(function(input, output) {

	###########################
	# Display Stage management
	###########################

	output$panel.controls <- renderUI({
		if (!is.null(input$select.phenotypes.continue) && input$select.phenotypes.continue != 0) {
		 	 renderSelectCovariates("")
		# Stage 3: Selection of phenotypes
		} else if (!is.null(input$select.cohorts.continue) && input$select.cohorts.continue != 0) {
			if (input$phenotypes.all.none == "all" || (is.null(input$phenotypes) && !identical(input$phenotypes.all.none,"none"))) {
				renderSelectPhenotypes(names(phenotypes),"all")
			} else if (input$phenotypes.all.none == "none") {
				renderSelectPhenotypes(c(),"none")
			} else {
				renderSelectPhenotypes(input$phenotypes,"")
			}
		## Stage 2: Selection of cohorts
		} else if (!is.null(input$gene.select.continue) && input$gene.select.continue != 0) {
			renderSelectCohorts(cohorts$name,cohorts$name,input$gene)
		## Stage 1: Selection of a gene or region
		} else if (is.null(input$gene)) {
			renderGeneSelect()
		} else if (input$gene %in% labels(test.genes)) {
				renderGeneAdjust(paste(input$gene),test.genes[[input$gene]]$chromosome,test.genes[[input$gene]]$start,test.genes[[input$gene]]$end,test.genes[[input$gene]]$strand,test.genes[[input$gene]]$description)
		} else if (!is.null(input$gene.lookup) && input$gene.lookup != 0) {
			gi <- isolate({gene.info.persist(input$gene)})
			if (class(gi)=="try-error") {
				renderGeneSelect(input$gene)
			} else {
				renderGeneAdjust(paste(input$gene),gi$chr,gi$GeneLowPoint,gi$GeneHighPoint,gi$ori,gi$genesummary)
			}
		} else {
			renderGeneSelect(input$gene)
		}
	})

	#output$panel.controls <- renderUI({panel})

	output$mainframe <- renderUI({

		mainPanel(
				if(!is.null(input$gene.selected)) {
					span(
						h1("Summary"),
						sidebarPanel(span(
							h4("Gene:"),
							h6("Gene name: ",ifelse(is.null(input$gene.selected),"Unknown", input$gene.selected)),
							h6("Chromosome: ", input$chromosome),
							h6("Range: ", input$range[[1]],"-",input$range[[2]]),
							h6("Strand: ", isolate({gene.strand(input)})),
							h6("Description: ", isolate({gene.description(input)}))
))
					)
				} else { span("") }
				,
				if (!is.null(input$select.cohorts.continue)) {  
					sidebarPanel(span(h4("Cohorts:"), p(paste(input$cohorts,sep=","))))
				} else { span("") }
				,
				if (!is.null(input$select.phenotypes.continue)) {
					if (!is.null(input$phenotypes)) {
						sidebarPanel(span(h4("Phenotypes:"), p(paste(sort(input$phenotypes),sep=","))))
					} else {
						sidebarPanel(span(h4("Phenotypes:"), p("None selected")))
					}
					#sidebarPanel(tableOutput("pheno.cohort.table"))
					#verbatimTextOutput("pheno.cohort.table")
				} else { span("") }
				,
				if (!is.null(input$select.covariates.continue)) {
					if (!is.null(input$covariates)) {
						sidebarPanel(span(h4("Covariates:"), p(paste(sort(input$covariates),sep=","))))
					} else {
						sidebarPanel(span(h4("Covariates:"), p("None selected")))
					}
				} else { span("") }
				, 
				if (!is.null(input$phenotypes) && !is.null(input$covariates)) {
					tableOutput("pheno.covar.table")
				} else { span("") }

		)
	})
})
