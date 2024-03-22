rm(list=ls())

#You can use this script to create a shiny app of the output of the SNPFiltering, it uses the "filters/pcaAll.tsv" file.
#Either you download this file from the server and run it on your local Rstudio or you do some magic that I do not understand.
treeSegments<-read.table("results/treeSegmentsAll.tsv",h=T) #This is important to fill in. It needs to point to the pcaAll.tsv file in filters (after downloading it to your local drive)
treeLabels<-read.table("results/treeLabelsAll.tsv",h=T)
#



library(shiny)
library(ggplot2)
library(ggrepel)
library(ggpubr)
vars=names(treeSegments)

shinyApp(
  ui = fluidPage(
    fluidRow(style = "padding-bottom: 20px;",
             column(4, selectInput('maf', 'maf',unique(treeSegments$maf))),
             column(4, selectInput('max_missing', 'max_missing',unique(treeSegments$max_missing)))
    ),
    fluidRow(
      plotOutput('cluster', height = "400px")  
    )
  ),
  
  server = function(input, output, session) {
    
    # Combine the selected variables into a new data frame
    selectedSegment= reactive({
      treeSegments[treeSegments$maf==input$maf&treeSegments$max_missing==input$max_missing,c("x","y", "xend","yend")]
    })
    
    selectedLabel= reactive({
      treeLabels[treeLabels$maf==input$maf&treeLabels$max_missing==input$max_missing,c("x","y","metaPop","nSNPs","metaPop","pop","label")]
    })
    
    output$cluster = renderPlot(height = 400, {
      ggplot(selectedSegment()) +
              geom_segment(aes(x = selectedSegment()[,"x"], y = selectedSegment()[,"y"], 
                               xend = selectedSegment()[,"xend"], yend=selectedSegment()[,"yend"]))+
              theme_pubr()+coord_cartesian(ylim=c(-3,max(selectedSegment()[,"y"])*1.15))+
            geom_text(data = selectedLabel(), aes(selectedLabel()[,"x"], 
                                                  selectedLabel()[,"y"], 
                                                  label = selectedLabel()[,"label"],
                                                  col=selectedLabel()[,"metaPop"]),
                                                  hjust = 1, angle = 90, size = 5)+
              xlab("")+
              ylab("genetic distance")+
              ggtitle(paste("number of SNPs =",unique(selectedLabel()[,"nSNPs"])))
      
    })
  },
  
  options = list(height = 500)
)
