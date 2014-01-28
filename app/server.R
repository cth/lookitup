library(shiny)
library(NCBI2R)

known.genes = list(
	"FTO" = list("chromosome" = 16, start=53737875, end=54148379)
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
	  textInput("gene", "Gene:", set.gene),
      selectInput("chromosome", "Chromosome:", 
                choices = c(as.character(seq(1,22)),"X","Y","MT"),set.chromosome),
	  sliderInput("range", "Range:",
                min = set.start, max = set.end, value = c(set.start,set.end)),
	actionButton("gene.select.continue", "Continue")
  )
}

##
# Step 2: Adjust selection after choice of gene
renderGeneAdjust <- function(set.gene="", set.chromosome,set.start, set.end) {
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
		selectInput("phenotypes.all.none", c("","All","None"), choices=c("","all","none"), selected=select.all.none),
		checkboxGroupInput("phenotypes", "Select phenotypes:", 
			choices=names(phenotypes),
			selected=select),
		actionButton("select.phenotypes.continue", "Continue")
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

# Define server logic
shinyServer(function(input, output) {

	# This function handles the main control logic of the application
	# It displays sidebars relevant to the particular stage of input
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
			print(summary)
			renderSelectCohorts(cohorts$name,cohorts$name,input$gene)
		## Stage 1: Selection of a gene or region
		} else if (is.null(input$gene)) {
			renderGeneSelect()
		} else if (input$gene %in% labels(known.genes)) {
			renderGeneAdjust(paste(input$gene),known.genes[[input$gene]]$chromosome,known.genes[[input$gene]]$start,known.genes[[input$gene]]$end)
		} else {
			renderGeneSelect(input$gene)
		}
	}) 


	# Create a data frame that contains the number of individuals by phenotype and cohorte 
	pheno.cohort.table <- isolate({
		if (!is.null(input$select.phenotypes.continue)) {  
			# FIXME: Precalculate instead and subset as needed subset as needed subset as needed 
			df <- head(sapply(input$cohorts, function(x) {
				cohort.id <- cohorts[[input$cohorts[x]]]
				sapply(names(phenotypes), function(y) { 
					sum(!is.na(phenotypes[phenotypes$studyid==cohort.id,which(colnames(phenotypes)==y)])) 
				})
			}))
			cat(dim(df))
			renderTable(df)
		} else {
			NULL
		}
	})

	pheno.covar.table <- reactive({
			renderTable(data.frame(matrix(nrow=length(input$phenotypes), ncol=length(input$covariates), data=rep(T, length(input$phenotypes) * length(input$covariates)))))
	})

	output$mainframe <- renderUI({
		mainPanel(
				if(!is.null(input$gene.select.continue)) {
					span(
						h1("Summary"),
						sidebarPanel(span(
							h4("Gene:"),
							h6("Gene name: ",ifelse(is.null(input$gene.selected),"Unknown", input$gene.selected)),
							h6("Chromosome: ", input$chromosome),
							h6("Range: ", input$range[[1]],"-",input$range[[2]])))
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
