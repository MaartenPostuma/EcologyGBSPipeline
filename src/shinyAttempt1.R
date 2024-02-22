rm(list=ls())
library(shiny)
dataSet<-read.table("results/test.tsv",h=T)
vars=names(dataSet)

dataSet

shinyApp(
  ui = fluidPage(
    fluidRow(style = "padding-bottom: 20px;",
             column(4, selectInput('xcol', 'X Variable', vars[3:5])),
             column(4, selectInput('ycol', 'Y Variable', vars[3:5],
                                   selected = vars[3])),
             column(4, selectInput('maf', 'maf',unique(dataSet$maf))),
             column(4, selectInput('max_missing', 'max_missing',unique(dataSet$max_missing)))
    ),
    fluidRow(
      plotOutput('pca', height = "400px")  
    )
  ),
  
  server = function(input, output, session) {
    
    # Combine the selected variables into a new data frame
    selectedData = reactive({
      dataSet[, c(input$xcol, input$ycol)]
    })
    
    
    output$pca = renderPlot(height = 400, {
      par(mar = c(5.1, 4.1, 0, 1))
      plot(selectedData())
    })
  },
  
  options = list(height = 500)
)
)
