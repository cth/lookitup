library(shiny)
library(NCBI2R)
library(memoise)

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



    
renderGenePanel <- function(input) {
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
renderCohortsPanel <- function(input) {
	wellPanel(
		if(is.null(input$cohorts)) {
			checkboxGroupInput("cohorts", "Select cohorts:",choices=cohorts$name, selected=cohorts$name)
		} else {
			checkboxGroupInput("cohorts", "Select cohorts:",choices=cohorts$name, selected=input$cohorts)
		})
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

renderPhenotypePanel <- function(input) {
	if (input$phenotypes.all.none == "all" || (is.null(input$phenotypes) && !identical(input$phenotypes.all.none,"none"))) {
		renderSelectPhenotypes(names(phenotypes),"")
	} else if (input$phenotypes.all.none == "none") {
		renderSelectPhenotypes(c(),"")
	} else {
		renderSelectPhenotypes(input$phenotypes,"")
	}
}


##
# Step 5: Select covariates
renderSelectCovariates <- function(set) {
	wellPanel(
		checkboxGroupInput("covariates", "Select covariates:", 
			choices=names(phenotypes),
			selected=set)
	)
}

renderCovariatesPanel <- function(input) {
	if (is.null(input$covariates)) {
		renderSelectCovariates("")
	} else {
		renderSelectCovariates(input$covariates)
	}
}

##
# Step 6: Select stratification
# TODO


### Summary functions
renderGeneSummary <- function(input) {
	if(!is.null(input$gene.lookup) && input$gene.lookup > 0) {
		sidebarPanel(span(
			h4("Gene:"),
			h6("Gene name: ",ifelse(is.null(input$gene),"Unknown", input$gene)),
			h6("Chromosome: ", input$chromosome),
			h6("Range: ", input$range[[1]],"-",input$range[[2]]),
			h6("Strand: ", isolate({gene.strand(input)})),
			h6("Description: ", isolate({gene.description(input)}))))
	} else { span("") }
}

itemListSummary <- function(items,caption) {
	if (is.null(items)) {	
		sidebarPanel(span(h4(caption), p("None selected")))
	} else {
		sidebarPanel(span(h4(caption), p(paste(items,sep=","))))
	}
}

renderCohortSummary <- function(input) { itemListSummary(input$cohorts, "Cohorts:") }
renderPhenotypeSummary <- function(input) { itemListSummary(input$phenotypes, "Phenotypes:") }
renderCovariateSummary <- function(input) { itemListSummary(input$covariates, "Covariates:") }

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
	output$mainframe <- renderUI({
		print(isolate({ifelse(is.null(input$top), "NULL", input$top)}))
		mainPanel(
				tabsetPanel(id="top",selected=ifelse(is.null(input$top),"Gene",input$top),
					tabPanel("Gene",
                                                 if(is.null(input$gene.lookup)){
                                                     renderGeneSelect(input)
                                             }else if(input$gene.lookup == 0){
                                                 renderGeneSelect(input)
                                             }else{
                                                 renderGeneLookup(input)
                                             }
#                                                 , renderGenePanel(input)
                                                 ),             
					tabPanel("Cohorts", renderCohortsPanel(input)),
					tabPanel("Phenotypes", renderPhenotypePanel(input)),
					tabPanel("Covariates", renderCovariatesPanel(input)),
					tabPanel("Summary",  
						renderGeneSummary(input),
						renderCohortSummary(input),
						renderPhenotypeSummary(input),	
						renderCovariateSummary(input))
				)
		)
	})
})
