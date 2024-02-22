rm(list=ls())

#You can use this script to create a shiny app of the output of the SNPFiltering, it uses the "filters/pcaAll.tsv" file.
#Either you download this file from the server and run it on your local Rstudio or you do some magic that I do not understand.
dataSet<-read.table("results/pcaAll.tsv",h=T) #This is important to fill in. It needs to point to the pcaAll.tsv file in filters (after downloading it to your local drive)
#



library(shiny)
library(ggplot2)
library(ggrepel)
library(ggpubr)
vars=names(dataSet)

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
      dataSet[dataSet$maf==input$maf&dataSet$max_missing==input$max_missing, c(input$xcol, input$ycol,"metaPop","pop","PCA1Var","PCA2Var","PCA3Var")]
    })
    
    
    output$pca = renderPlot(height = 400, {
      ggplot(data=selectedData(),aes(x=selectedData()[,input$xcol],y=selectedData()[,input$ycol],
                                     col=selectedData()[,"metaPop"],label=selectedData()[,"pop"]))+
        theme_pubclean()+
        geom_point()+
        geom_text_repel(col="black")+
        xlab(paste0(input$xcol," (",round(unique(selectedData()[,paste0(input$xcol,"Var")],1)*100),"%)"))+
        ylab(paste0(input$ycol," (",round(unique(selectedData()[,paste0(input$ycol,"Var")],1)*100),"%)"))
      
    })
  },
  
  options = list(height = 500)
)
# 
# ggplot(dataSet,aes(x=EV1,y=EV2))+
#   geom_point()+
#   facet_grid(max_missing~maf)

