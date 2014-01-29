library(shiny)

# The UI is mostly dynamically configured from the server-side.
# In this UI definition, we only setup an empty page to be filled
# with stuff. All this happens in server.R

shinyUI(bootstrapPage(uiOutput("mainframe")))
