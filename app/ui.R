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
                                     uiOutput("session.key"),
                                     actionButton("save.session", "Save session"),
                                     actionButton("restore.session","Restore session")),
                                 span(h3("Summary of selections:"),uiOutput("summary.tab")))
                             ),

                    tabPanel("Input", 
						span(
					  	  div(class="span4",
					    	tags$form(class="well",
					        actionButton("lookupButton","Lookup"),
					        actionButton("add.input","Add"),
					        uiOutput("exploratorium.tab"),
					        actionButton("addRangeButton","Add Range"))),
					  	  div(class="span3",
							uiOutput("input.tab")))),
                    tabPanel("Cohorts", uiOutput("cohorts.tab")),
                    tabPanel("Phenotypes", uiOutput("phenotypes.tab")),
                    tabPanel("Covariates", uiOutput("covariates.tab")),
                    tabPanel("Analysis", uiOutput("analysis.tab")),
                    tabPanel("Results", uiOutput("results.tab"))))))

