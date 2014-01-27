library(shiny)
library(NCBI2R)

known.genes = list(
	"FTO" = list("chromosome" = 16, start=53737875, end=54148379)
)

cohorts <- c("inter99", "Helbred06", "Steno", "Vejle")

phenotypes <- read.table("MergedPhenotypes_13dec2013.txt",head=T)

renderSelectPhenotypes <- function() {
	sidebarPanel(
		checkboxGroupInput("phenotypes", "Select phenotypes:", 
			choices=names(phenotypes),
			selected=names(phenotypes)),
		actionButton("select.phenotypes.continue", "Continue")
	)
}

renderSelectCohorts <- function(cohorts,cohorts.selected,gene) {
	sidebarPanel(
		checkboxGroupInput("cohorts", "Select cohorts:", 
			choices=cohorts,
			selected=cohorts.selected),
		actionButton("select.cohorts.continue", "Continue")
	)
}

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


# Define server logic required to generate and plot a random distribution
shinyServer(function(input, output) {

#	observe({
#		if (!is.null(input$lookup)) {
#			if (input$lookup != 0) {
#				cat(paste("LOOKUP", input$lookup, "\n"))
#			}
#		}
#	})

	output$panel.controls <- renderUI({
		if (!is.null(input$select.cohorts.continue) && input$select.cohorts.continue != 0) {
			renderSelectPhenotypes()
		} else if (!is.null(input$gene.select.continue) && input$gene.select.continue != 0) {
			print(summary)
			renderSelectCohorts(cohorts,cohorts,input$gene)
		} else if (is.null(input$gene)) {
			renderGeneSelect()
		} else if (input$gene %in% labels(known.genes)) {
			renderGeneAdjust(paste(input$gene),known.genes[[input$gene]]$chromosome,known.genes[[input$gene]]$start,known.genes[[input$gene]]$end)
		} else {
			renderGeneSelect(input$gene)
		}
	}) 

	output$mainframe <- renderUI({
		
		mainPanel(
			span(
				if(!is.null(input$gene.select.continue)) {
					span(
						h1("Summary"),
						h2("Selected gene: ",ifelse(is.null(input$gene.selected),"Unknown", input$gene.selected)),
						h3("Chromosome: ", input$chromosome),
						h3("Range: ", input$range[[1]],"-",input$range[[2]])
					)
				} else { span("") }
				,
				if (!is.null(input$select.cohorts.continue)) {  
					span(h4("Cohorts:"), p(paste(input$cohorts,sep=",")))
				} else { span("") }
				,
				if (!is.null(input$select.phenotypes.continue)) {  
					span(h4("Phenotypes:"), p(paste(sort(input$phenotypes),sep=",")))
				} else { span("") }

			)
		)
	})
})
