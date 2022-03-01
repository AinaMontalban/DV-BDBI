## ----setup, include=FALSE------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, fig.width = 6, fig.height = 4)
  #install.packages("prettydoc")


## ----message=FALSE, warning=FALSE----------------------------------------
library(readxl)
data <- read_excel("nature22402-s3.xlsx", skip = 1)


## ------------------------------------------------------------------------
colnames(data)[1] <- "Position"
colnames(data)[2] <- "Ancestral.Allele"
colnames(data)[3] <- "Derived.Allele"
colnames(data)[4] <- "Ancestral.Freq"
colnames(data)[5] <- "Derived.Freq"
colnames(data)[8] <- "Codon.position"
colnames(data)[10] <- "Codon.change"
colnames(data)[11] <- "Amino.Acid.Change"


## ------------------------------------------------------------------------
#Allele
class(data$Ancestral.Allele)
class(data$Derived.Allele)

#Allele frequency within outbreak
class(data$Ancestral.Freq)
class(data$Derived.Freq)

#Number of alleles
class(data$Minor)
class(data$Total)

#Impact on ZIKV proteins
class(data$Codon.position)
class(data$Degeneracy) #factor?  
data$Degeneracy <- as.factor(data$Degeneracy)
class(data$Codon.change)
class(data$Amino.Acid.Change)
class(data$Protein) #can be a factor?
data$Protein <- as.factor(data$Protein)


## ------------------------------------------------------------------------
subset_data_alleles <- data[,c(1,4,5)]

library(reshape2)
long_subset_data_alleles <- melt(subset_data_alleles, 
                                 id.vars ="Position", 
                                 variable.name = "A_D",
                                 value.name = "Frequency")



## ------------------------------------------------------------------------
subset_data_alleles_nuc <- data[,c(1,2,3)]

long_subset_data_alleles_nuc <- melt(subset_data_alleles_nuc, 
                                 id.vars ="Position", 
                                 variable.name = "A_D",
                                 value.name = "Nucleotide")



## ------------------------------------------------------------------------
library(ggplot2)


## ------------------------------------------------------------------------
ggplot(long_subset_data_alleles, aes(Frequency, fill=A_D)) + 
  geom_histogram(position = "dodge") +  xlab("Allele Frequency") + scale_x_continuous(breaks = seq(0,1, by = 0.1)) + 
  theme_minimal() + 
  scale_fill_viridis_d(name="Allele Type", labels=c("Ancestral", "Derived")) + ylab("Count") + labs(title="Allele Frequency Distribution")


## ------------------------------------------------------------------------
ggplot(long_subset_data_alleles_nuc, aes(Nucleotide, fill=A_D)) + 
  geom_bar(position = "dodge") + theme_minimal() + 
  scale_fill_viridis_d(name="Allele Type", labels=c("Ancestral", "Derived")) + labs(title = "Nucleotide Distribution", y="Count")


## ------------------------------------------------------------------------
ggplot(data, aes(x=factor(1), fill=Degeneracy)) + geom_bar(width = 1) + coord_polar("y", start = 0) + scale_fill_viridis_d() + theme_minimal() + theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        panel.grid  = element_blank(), legend.position = "left") + geom_text(stat='count', aes(label=paste(100*round((..count..)/sum(..count..),3), "%")), hjust=1, color="red", size=3) + labs(x="", y="", title="Degeneracy Distribution") 




## ------------------------------------------------------------------------

ggplot(data, aes(Protein, fill=Protein)) + geom_bar(show.legend = FALSE) + 
  theme_minimal() + scale_fill_viridis_d() + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_text(stat='count', aes(label=..count..), hjust=-0.3)


## ------------------------------------------------------------------------
AA <- c("A","R","N","D","C","Q","E","G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V")
length(AA)
M <- matrix(0, nrow = 20, ncol = 20)
colnames(M) <- AA
rownames(M) <- AA
data$Amino.Acid.Change <- gsub("->", "-", data$Amino.Acid.Change)


## ------------------------------------------------------------------------
for (i in c(1:20)){
  for(j in c(1:20)){
    s <- paste(rownames(M)[i], colnames(M)[j],sep=" - ")
    if (s %in% data$Amino.Acid.Change){
    for (ele in data$Amino.Acid.Change){
      if (ele == s){
        M[i,j] <- M[i,j] +  1
      }else{
        next;}}}}
}


## ------------------------------------------------------------------------
long_data_M <- melt(M)
colnames(long_data_M) <- c("Original", "Change", "Count")
ggplot(long_data_M, aes(x=Change, y=Original, fill=Count)) + geom_tile() + 
    scale_fill_gradient(low="white", high="purple4")   + 
    geom_text(aes(label = (ifelse(Count > 0, Count, ""))), color="white") + labs(title=" Distribution of the AA changes")


## ------------------------------------------------------------------------
ggplot(data, aes(x=as.factor(Position), y = Derived.Freq, fill = as.factor(Degeneracy)))  + geom_col() + 
  scale_fill_viridis_d(name = "Degeneracy") + 
  labs(x= "Position of the genomes", y = "Minor Allele Frequency (MAF)", title = "MAF of nonsynonymous mutations related with degeneracy") + 
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) + 
  annotate(geom="text", x=59, y=1, label="2853") + annotate(geom="text", x=70, y=0.4, label="3534") + annotate(geom="text", x=170, y=0.37, label="8631")   + annotate(geom="text", x=193, y=0.38, label="10301")  


## ------------------------------------------------------------------------
p1 <- ggplot(data, aes(Protein,y=Derived.Freq, group=Derived.Freq, fill=Degeneracy, 
                       text=paste("Protein: ", Protein,"\n", "Derived Freq: ", Derived.Freq,"\n", 
                                  "Degeneracy: ", Degeneracy,"\n" ,"Genome Position: ", Position))) + 
  geom_col(position = "dodge") +
  labs(x= "Protein", y = "Minor Allele Frequency (MAF)", title = "Protein related with degeneracy and Minor Allele Frequency") +  
  theme_minimal() + theme(axis.text.x = element_text(angle=60, hjust =1)) + scale_fill_viridis_d()


p1


## ------------------------------------------------------------------------
data$mutation <- paste(data$Ancestral.Allele, data$Derived.Allele)
ggplot(data, aes(mutation, fill = as.factor(Degeneracy))) +geom_bar() + 
  scale_fill_viridis_d( name = "Degeneracy") + scale_x_discrete(labels=c("A C" = "A -> C", "A G" = "A -> G", "A T" = "A -> T", "C A" = "C -> A", "C G" = "C -> G", "C T" = "C -> T", "G A" = "G -> A", "G C" = "G -> C", "G T" = "G -> T", "T A" = "T -> A", "T C" = "T -> C","T G" = "T -> G" )) + 
  theme_minimal() + theme(axis.text.x = element_text(angle=60, hjust =1)) +
  labs(x="Mutation", y="Frequency", title = "Mutation distribution related with degeneracy")


## ----message=FALSE, warning=FALSE, paged.print=FALSE---------------------
library(plotly)


## ------------------------------------------------------------------------

p2 <- ggplot(long_data_M, aes(x=Change, y=Original, fill=Count,  text = paste('Original: ', Original,"\n", 'Change: ', Change,"\n", 'Count: ', Count))) + geom_tile() + 
    scale_fill_gradient(low="white", high="purple4")   + 
    geom_text(aes(label = (ifelse(Count > 0, Count, ""))), color="white") + labs(title=" Distribution of the AA changes")

ggplotly(p2, tooltip = "text")


## ----eval=TRUE-----------------------------------------------------------

library(shiny)
library(dplyr)
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
                                          Degeneracy,"\n" ,"Genome Position: ", Position,"\n" ,"AA Change: ", Amino.acid.change))) +
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


## ----include=FALSE-------------------------------------------------------
colnames(data)[1] <- "Position"
colnames(data)[2] <- "Ancestral"
colnames(data)[3] <- "Derived"
colnames(data)[4] <- "Allele.Freq"
colnames(data)[5] <- "Derived.Freq"
colnames(data)[8] <- "Codon.position"
colnames(data)[10] <- "Codon.change"
colnames(data)[11] <- "Amino.acid.change"


## ------------------------------------------------------------------------
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

