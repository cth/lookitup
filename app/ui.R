library(shiny)

# The UI is mostly dynamically configured from the server-side.
# In this UI definition, we only setup an empty page to be filled
# with stuff. All this happens in server.R

shinyUI(bootstrapPage(
	tags$head(
		tags$script(src = "js/session.js")
	),
	mainPanel(
		tabsetPanel(id="top",
			tabPanel("Session", 
				mainPanel(
				wellPanel(
					#p(uiOutput("session.status")),
					uiOutput("session.key"),
					actionButton("save.session", "Save session"),
					actionButton("restore.session","Restore session")),
					span(h3("Summary of selections:"),uiOutput("summary.tab")))
			),
			tabPanel("Exploratorium", 
                                 mainPanel(
                                     textInput("lookup","Search:"),
                                     submitButton("Lookup")
                                     ),
                                 mainPanel(
                                     h4("Summary:"),
                                     uiOutput("lookupPrint")
                                     ) 
                                 ),
			tabPanel("Cohorts", uiOutput("cohorts.tab")),
			tabPanel("Phenotypes", uiOutput("phenotypes.tab")),
			tabPanel("Covariates", uiOutput("covariates.tab")),
			tabPanel("Results", 
				tabsetPanel(id="results.tabset",
					tabPanel("test1", h1("test1")),
					tabPanel("test2", h1("test2"))))))))
