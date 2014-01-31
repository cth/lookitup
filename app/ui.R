library(shiny)

# The UI is mostly dynamically configured from the server-side.
# In this UI definition, we only setup an empty page to be filled
# with stuff. All this happens in server.R

shinyUI(bootstrapPage(
	mainPanel(
		tabsetPanel(id="top",
			tabPanel("Gene", uiOutput("gene.tab")),
			tabPanel("Cohorts", uiOutput("cohorts.tab")),
			tabPanel("Phenotypes", uiOutput("phenotypes.tab")),
			tabPanel("Covariates", uiOutput("covariates.tab")),
			tabPanel("Summary",  uiOutput("summary.tab"))))))

