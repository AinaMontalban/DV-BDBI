library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)

library(readxl)
data <- read_excel("nature22402-s3.xlsx", skip = 1)

colnames(data)[1] <- "Position"
colnames(data)[2] <- "Ancestral"
colnames(data)[3] <- "Derived"
colnames(data)[4] <- "Allele.Freq"
colnames(data)[5] <- "Derived.Freq"
colnames(data)[8] <- "Codon.position"
colnames(data)[10] <- "Codon.change"
colnames(data)[11] <- "Amino.acid.change"


data$Degeneracy <- as.factor(data$Degeneracy)
data$Protein <- as.factor(data$Protein)


ui <- fluidPage(
  verticalLayout(
    h1("MAF of nonsynonymous mutations"),
    sliderInput("Pos", "Position:",
                min = 0, max = 11000, value = c(3994,7064)
    ),
    plotOutput("plot1", brush = "plot_brush", click = "plot_click"),
    br(),
    
    h4("Brushed cols"),
    fluidRow(
      column(width = 10, textOutput("info")),
    column(width = 6, tableOutput("infoBrush"))
    
    ),
    column(width=6, h4("Click cols")),
    
    column(width = 6, tableOutput("infoClick"))
  )
)

server <- function(input, output) {
  output$plot1 <- renderPlot({
    
    filtered <- data %>% 
      filter(Position >= input$Pos[1], 
             Position <= input$Pos[2])
    
    
    p <- ggplot(filtered, aes(x=as.factor(Position), y = Derived.Freq, fill = as.factor(Degeneracy)))  + 
      geom_col() + 
      scale_fill_viridis_d(name = "Degeneracy") + 
      labs(x= "Position of the genomes", y = "Minor Allele Frequency (MAF)", title = "MAF of nonsynonymous mutations related with degeneracy") + 
      theme(axis.text.x = element_blank(), axis.ticks = element_blank())  
    p
    
  })
  
  
  
  output$infoBrush <- renderTable({
    # With base graphics, need to tell it what the x and y variables are.
    tbl_brush <- brushedPoints(data, input$plot_brush, xvar="Position", yvar = "Derived.Freq")
    tbl_brush[,c("Position", "Ancestral", "Derived", "Allele.Freq","Derived.Freq", "Codon.position", "Degeneracy", "Codon.change","Amino.acid.change", "Protein")]
  })
  
  output$info <- renderText({
    tbl_brush <- brushedPoints(data, input$plot_brush, xvar="Position", yvar = "Derived.Freq")
    
    paste("Number of results in the table: ", nrow(tbl_brush))
  })
  
  output$infoClick <- renderTable({
    # With base graphics, need to tell it what the x and y variables are.
    tbl_click <- nearPoints(data, input$plot_click, xvar="Position", yvar = "Derived.Freq")
    tbl_click[,c("Position", "Ancestral", "Derived", "Allele.Freq","Derived.Freq", "Codon.position", "Degeneracy", "Codon.change","Amino.acid.change", "Protein")]
    
  })
  
  
  
  
}
# Run the application 
shinyApp(ui = ui, server = server)


