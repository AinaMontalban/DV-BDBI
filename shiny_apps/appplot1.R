library(shiny)
library(dplyr)
library(ggplot2)
library(plotly)
library(readxl)
data <- read_excel("nature22402-s3.xlsx", skip = 1)

  
colnames(data)[1] <- "Position"
colnames(data)[2] <- "Ancestral.Allele"
colnames(data)[3] <- "Derived.Allele"
colnames(data)[4] <- "Ancestral.Freq"
colnames(data)[5] <- "Derived.Freq"
colnames(data)[8] <- "Codon.position"
colnames(data)[10] <- "Codon.change"
colnames(data)[11] <- "Amino.Acid.Change"

data$Degeneracy <- as.factor(data$Degeneracy)
class(data$Protein) #can be a factor?
data$Protein <- as.factor(data$Protein)

ui <- fluidPage(
  titlePanel(HTML('Protein related with degeneracy and MAF')),
  sidebarLayout(
    
    sidebarPanel(
      checkboxGroupInput("protInput", "Protein: ", choices = c("2K" ,"capsid", "envelope" ,"membrane", "NS1", "NS2A", "NS2B", "NS3","NS4A", "NS4B",
                                                               "NS5", "propeptide"), selected = c("NS1", "NS5", "capsid")),
      checkboxGroupInput("degInput", "Degeneracy",
                         choices = c("1", "2" ,"3"),
                         selected = c("1")),
    ),
    mainPanel(plotlyOutput("plot"),
              tableOutput("table")
              
    )
    
  ))

# Define server logic required to draw a histogram
server <- function(input, output) {
  
  output$plot <- renderPlotly({
    
    filtered <- data %>% 
      filter(Degeneracy %in% input$degInput, 
             Protein %in% input$protInput)
    
    p1 <- ggplot(filtered, aes(as.factor(Protein),Derived.Freq, fill=as.factor(Degeneracy), 
                               text=paste("Protein: ", Protein,"\n", "Derived Freq: ", Derived.Freq,"\n", "Degeneracy: ",
                                          Degeneracy,"\n" ,"Genome Position: ", Position,"\n" ,"AA Change: ", Amino.Acid.Change))) +
      geom_col(position = "dodge") +
      labs(x= "Protein", y = "Minor Allele Frequency (MAF)", title = "Protein related with degeneracy and Minor Allele Frequency") +  
      theme_minimal() + theme(axis.text.x = element_text(angle=60, hjust =1)) + scale_fill_viridis_d()
    
    ggplotly(p1, tooltip = "text")
  })
  
  
  output$table <- renderTable({
    filtered <- data %>% 
      filter(Degeneracy %in% input$degInput, 
             Protein %in% input$protInput)
    filtered
    
  })
  
}



# Run the application 
shinyApp(ui = ui, server = server)