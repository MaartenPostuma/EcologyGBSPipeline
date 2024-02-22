rm(list=ls())
library(shiny)
library(ggplot2)
dataSet<-read.table("results/pcaAll.tsv",h=T)
vars=names(dataSet)
dataSet$maf
dataSet$max_missing
dataSet

shinyApp(
  ui = fluidPage(
    fluidRow(style = "padding-bottom: 20px;",
             column(4, selectInput('xcol', 'X Variable', vars[3:5],selected = vars[3])),
             column(4, selectInput('ycol', 'Y Variable', vars[3:5],selected = vars[4])),
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
      dataSet[dataSet$maf==input$maf&dataSet$max_missing==input$max_missing, c(input$xcol, input$ycol,"metaPop")]
    })
    
    
    output$pca = renderPlot(height = 400, {
      ggplot(data=selectedData(),aes(x=selectedData()[,input$xcol],y=selectedData()[,input$ycol]))+
        geom_point()
    })
  },
  
  options = list(height = 500)
)
# 
# ggplot(dataSet,aes(x=EV1,y=EV2))+
#   geom_point()+
#   facet_grid(max_missing~maf)
