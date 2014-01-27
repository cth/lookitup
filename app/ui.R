library(shiny)

shinyUI(pageWithSidebar(
  	headerPanel("Lookup"),
	uiOutput("panel.controls"),
	uiOutput("mainframe")
  #	mainPanel(
#		h3(textOutput("summary"))
#  	)
))
