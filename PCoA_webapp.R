# Author: Nikeeta S. Chavan
# Date: 21 November 2020
# Description: Webapp for normalization and Principle Component Analysis 
# of transcriptome data. 

# check and install required packages
cran.pkgs <- c('shiny','vegan','plotly','ape')
if (length(setdiff(cran.pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(cran.pkgs, rownames(installed.packages())))  
}
bioc.pkgs <- c('DESeq2')
if (length(setdiff(bioc.pkgs, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(bioc.pkgs, rownames(installed.packages())))  
}

# load packages
library(shiny)
library(DESeq2)
library(ape)
library(vegan)
library(plotly)

# UI
ui <- fluidPage(
  titlePanel("Principal Component Analysis (PCA)"),
  sidebarLayout(
    sidebarPanel(
      # upload file of counts data
      fileInput("countsFile", "Counts file (Comma separated)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",".txt")),
      # Horizontal line
      tags$hr(),
      # upload metadata file
      fileInput("metaFile","Metadata file (Comma separated)",
                multiple = FALSE,
                accept = c("text/csv",
                           "text/comma-separated-values,text/plain",
                           ".csv",".txt")),
      
      tags$hr()
      
    ),
    
    # Main panel
    mainPanel("This app performs Principal Component Analysis (PCoA) on transcriptome data and visualize ordination in 2D plot. The steps involved nomralization using DESeq2 package, calculation of Bray-Curtis distance matrix using vegan package and subsequently calculation of PCA using ape package.",
      tags$hr(),
      tags$hr(),
      uiOutput("columns"),
      plotlyOutput(outputId = "plot")
    )
    
  )
)

# Define server 
server <- function(input, output) {
  # PCA object
  pca_obj <- reactive({
    # read metadata
    inFileMeta <- input$metaFile
    if (is.null(inFileMeta))
      return(NULL)
    df_meta <- read.csv(inFileMeta$datapath, header = T, row.names = 1)
    # read counts file
    inFileCounts <- input$countsFile
    if (is.null(inFileCounts))
      return(NULL)
    df_counts <- read.csv(inFileCounts$datapath, header = T, row.names = 1)
    # create deseq2 object
    dds <- DESeqDataSetFromMatrix(countData = df_counts,colData = df_meta,design = ~Samples)
    # estimate size factors
    dds <- estimateSizeFactors(dds)
    # perform normalization
    normalized_counts <- counts(dds, normalized=TRUE)
    # Calculate distance
    d <- vegdist(x = t(normalized_counts),method = 'bray')
    # perform PCA analysis
    pca.res <- pcoa(D = d)
    # Percent variations explained
    axis1.pve <- round(pca.res$values[1,'Relative_eig']*100,2)
    axis2.pve <- round(pca.res$values[2,'Relative_eig']*100,2)
    # Formating dataframe
    V <- as.data.frame(pca.res$vectors)
    V <- merge(x = V,y = df_meta,by=0)
    write.table(x = V,file = 'temp.txt',quote = F,sep = '\t')
    return(list(axis1.pve,axis2.pve,V))
  })
  
  # Metadata option for dropdown choices
  output$columns <- renderUI({
    inFile <- input$metaFile
    if (is.null(inFile))
      return()
    dat = read.csv(inFile$datapath, header = T,row.names = 1)
    selectInput(inputId = "metadataValue", label = "Select metadata variable", choices = colnames(dat),selected = NULL)
  })
  
  # PCA plot 
  output$plot <- renderPlotly({
    if(is.null(input$metaFile))
      return(NULL)
    df <- as.data.frame(pca_obj()[[3]])
    f <- list(face='bold')
    x <- list(title=paste("Axis.1 (",pca_obj()[[1]],"%)"),font=f)
    y <- list(title=paste("Axis.2 (",pca_obj()[[2]],"%)"),font=f)
    plot_ly(data = df,x=~Axis.1,y=~Axis.2,color = df[,input$metadataValue],colors = 'Set1',mode="markers") %>% layout(xaxis=x,yaxis=y)
    }
  )
}

# Create app
shinyApp(ui, server)