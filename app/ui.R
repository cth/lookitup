library(shiny)

# The UI is mostly dynamically configured from the server-side.
# In this UI definition, we only setup an skeleton page with tabs to be filled
# with from the server side.. 
# All this happens from server.R

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
			tabPanel("Exploratorium", mainPanel(
                            actionButton("lookupButton","Lookup"),
                            uiOutput("exploratorium.tab"))),
			tabPanel("Cohorts", uiOutput("cohorts.tab")),
			tabPanel("Phenotypes", uiOutput("phenotypes.tab")),
			tabPanel("Covariates", uiOutput("covariates.tab")),
			tabPanel("Analysis", uiOutput("analysis.tab")),
			tabPanel("Results", uiOutput("results.tab"))))))
		#		tabsetPanel(id="results.tabset",
		#			tabPanel("test1", h1("test1")),
		#			tabPanel("test2", h1("test2"))))))))
