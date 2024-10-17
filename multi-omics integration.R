#loading packages
library(TCGAbiolinks) #for accessing TCGA data
library(SummarizedExperiment)
library(DESeq2) #for DEA
library(tidyverse) #for data wrangling and visualization
library(org.Hs.eg.db) #for gene symbol mapping
library(biomaRt)
library(survival) #for survival anaylsis based on exp data
library(survminer)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) #for annotation of methylation data
library(sesame) #to create summarizedexperiment dataset of methylation data
library(sesameData)
library(GenomicRanges) #to get gene coordinates for CNV analysis

#getting TCGA-LUSC data
query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Transcriptome Profiling",
  experimental.strategy = "RNA-Seq",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts")

setwd("C:/Users/Spectre/Desktop")

#dowloading query and creating summarizedexperiment object
GDCdownload(query)
data <- GDCprepare(query)

#saving summarizedexperiment dataset as RDS file
saveRDS(data, file = "summarized_experiment", refhook = NULL)
data <- read_rds("summarized_experiment", refhook = NULL)

#normalization and DEA
dds <- DESeqDataSet(data, design = ~ definition)
dds <- DESeq(dds)

#DESeqDataSet alters duplicate names by adding numbers at the end
#removing those numbers for proper query to ENSMBL ID mapping
rownames(dds) <- sub("\\..*", "", rownames(dds))
rownames(res) <- sub("\\..*", "", rownames(res))

#saving dds as RDS file
saveRDS(dds, file = "dds_object", refhook = NULL)
dds <- read_rds("dds_object", refhook = NULL)

#getting DESeq results from dds
res <- results(dds)
#extracting normalized counts from dds
count <- counts(dds, normalized = TRUE)

#taking user input
input <- toupper(readline("Enter gene name: "))

#function to validate input and map it to ensbl IDs
validate_map <- function(gene_name, dds) {
  ID <- mapIds(org.Hs.eg.db,
               keys = gene_name,
               column = "ENSEMBL",
               keytype = "SYMBOL",
               multiVals = "first")
  
  ifelse (is.na(ID) || !ID %in% rownames(dds), stop("INVALID INPUT! make sure you are entering the correct gene symbol"), return(ID))
}

#function to fetch user defined gene for status evaluation
get_gene <- function(gene_name, dds, res) {
  
  ID <- validate_map(gene_name, dds)
  
  #get the DEA result of corresponding ID from DESeq result
  #returning DEA result
  return(res[ID, ])
}

#function to get counts data for user defined gene
get_expression_data <- function(gene_name, dds, count) {
  #mapping
  ID <- validate_map(gene_name, dds)
  gene_expression <- count[ID, ]
  
  if (is.null(gene_expression)) {
    stop(paste("No expression data found for gene:", gene_name))
  }
  
  gene_df <- data.frame(expression = as.numeric(gene_expression),
                        sample_type = colData(dds)$definition,
                        sample_id = colnames(count)
                        )
  return(gene_df)
}

####################################
#function call to get inputted gene
gene_result <- as.data.frame(get_gene(input, dds, res))
#Function call to get count data of inputted gene
expression_result <- get_expression_data(input, dds, count)

#drawing box plot to check expression levels of user defined gene
#in normal vs tumor samples
bplot <- ggplot(expression_result,
       aes(x = sample_type,
           y = expression,
           fill = sample_type)) +
  geom_boxplot(alpha =0.5) +
  geom_jitter(shape = 16,
              position = position_jitter(0.1),
              color = "black",
              size = 1) +
  labs(title = paste("Expression of", input, "in Normal vs Tumor samples"),
       x = "Sample Type",
       y = "Normalized Counts") +
    theme_minimal() +
    theme(legend.position = "none")

#Now
#performing survival analysis
survival_analysis <- function(gene_name, clinical_data, dds, count) {
   
  expression_result <- get_expression_data(gene_name, dds, count)
  
  #getting sample IDs from expression_result
  expression_samples <- expression_result$sample_id
  
  #matching clinical data to expression data based on sample ids
  clinical_data_f <- clinical_data[clinical_data$barcode %in% expression_samples, ]
  
  if (nrow(clinical_data_f) == 0) {
    stop("No matching clinical data found for the provided expression samples.")
  }
  
  #creating a df with clinical data and expression
  survival_df <- data.frame(
    expression = expression_result$expression,
    time = as.numeric(clinical_data_f$days_to_death),
    status = ifelse(clinical_data_f$vital_status == "Dead", 1, 0),
    sample_type = expression_result$sample_type 
  )
  
  #removing NAs
  survival_df <- survival_df[!is.na(survival_df$time) & !is.na(survival_df$status), ]
  
  #categorizing high and low groups based on median value of expression
  survival_df$expression_group <- ifelse(survival_df$expression > median(survival_df$expression), "high", "low")
  
  #performing kaplan-meier analysis
  surv_obj <- Surv(time = survival_df$time, event = survival_df$status)
  fit <- survfit(surv_obj ~ expression_group, data = survival_df)
  
  #survival plot
  survival_plot <- ggsurvplot(fit,
                              data = survival_df,
                              pval = T,
                              risk.table = T,
                              ggtheme = theme_minimal(),
                              title = paste("Survival analysis"),
                              legend.title = "Expression group",
                              palette = c("cyan", "pink"))
  return(survival_plot)
}

##################################
#getting clinical data
clinical_data <- as.data.frame(data@colData)


#function call
splot <- survival_analysis(input, clinical_data, dds, count)

###################################
# creating a query for methylation data
methylation <- GDCquery(project = 'TCGA-LUSC',
                     data.category = 'DNA Methylation',
                     platform = 'Illumina Human Methylation 450',
                     access = 'open',
                     data.type = 'Methylation Beta Value')

#downloading and preparing query
GDCdownload(methylation)
methylation_data <- GDCprepare(methylation, summarizedExperiment = T)

#saving RDS for later use
saveRDS(methylation_data, file = "methylation_data.rds")
methylation_data <- readRDS("methylation_data.rds")


#getting annotations from annotation package
annotation <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

#adjusting multivalued cells
split <- c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")
annotation <- annotation %>%
  separate_rows(!!!syms(split), sep = ";")

#extracting beta values
beta_values <- assay(methylation_data)

#converting beta to m-values
mval <- t(apply(beta_values, 1, function(x) log2(x/(1-x))))
saveRDS(mval, file = "mval.rds")
mval <- readRDS("mval.rds")

#Function to get gene associated probes and generate boxplots to visulize m-values by sample types
meth_status <- function(input, mval, annotation){
  
  #filtering out probes associated to the user inputted gene
  gene_probes <- annotation[annotation$UCSC_RefGene_Name == input, ]
  
  #getting probes associated to gene
  gene_mval <- mval[rownames(mval) %in% gene_probes$Name, ]
  
  #getting sample information
  sample_info <- as.data.frame(colData(methylation_data))
  definition <- sample_info[, c("barcode", "definition")]
  
  #adding Sample Types to mvals associated to inputted gene
  gene_mval_df <- as.data.frame(t(gene_mval))
  gene_mval_df$definition <- definition$definition
  
  #melting data to long format
  gene_mval_long <- gene_mval_df %>%
    gather(key = "CpG_Probe", value = "Mval", -definition)
  
  # Check for NA or infinite values in Mval
  sum(is.na(gene_mval_long$Mval))  # Count NAs
  sum(!is.finite(gene_mval_long$Mval))  # Count non-finite values
  
  # Remove rows with NA or non-finite values
  gene_mval_long <- gene_mval_long %>%
    filter(is.finite(Mval) & !is.na(Mval))
  
  
  # Plot the boxplots for each CpG probe
  meth_bplot <- ggplot(gene_mval_long, aes(x = definition, y = Mval, fill = definition)) +
    geom_boxplot(alpha = 0.5) +
    geom_jitter(shape = 16,
                alpha = 0.4,
                position = position_jitter(0.1),
                color = "black",
                size = 0.5) +
    facet_wrap(~ CpG_Probe, scales = "free") +
    labs(title = paste("Methylation status for", input, "associated Probes by Sample Type"),
         x = "Sample Type", y = "M-value") +
    theme_minimal() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          strip.text = element_text(size=10))
  
  return(meth_bplot)
}

################################
#function call to inspect methylation status of gene
<- meth_res <- meth_status(input, mval, annotation)
################################

# Creating a query for CNV data
cnv_query <- GDCquery(
  project = "TCGA-LUSC",
  data.category = "Copy Number Variation",
  data.type = "Masked Copy Number Segment",
  access = "open"
)

# Downloading and preparing the CNV data
GDCdownload(cnv_query)
cnv_data <- GDCprepare(cnv_query, summarizedExperiment = TRUE)

#RDS for later use
saveRDS(cnv_data, "cnv_data.rds")
cnv_data <- read_rds("cnv_data.rds")

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

#function for CNV analysis
cnv_analysis <- function(input, dds, cnv_data, mart) {
  
  #getting ensemble ID for the inputted symbol
  gene_ID <- validate_map(input, dds)
  
  #getting inputted gene coordinates
  gene_coords <- getBM(attributes = c("ensembl_gene_id", "chromosome_name", "start_position", "end_position", "strand"),
                       filters = "ensembl_gene_id",
                       values = gene_ID,
                       mart = mart)
  
  # Create GRanges object for the input gene using its coordinates
  gene_input <- GRanges(seqnames = paste0("chr", gene_coords$chromosome_name),
                        ranges = IRanges(start = gene_coords$start_position, end = gene_coords$end_position),
                        strand = ifelse(gene_coords$strand == 1, "+", "-"))
  
  # Creating GRanges object for CNV data
  cnv_granges <- GRanges(
    seqnames = paste0("chr", cnv_data$Chromosome),  # Ensure chromosomes are prefixed with 'chr'
    ranges = IRanges(start = cnv_data$Start, end = cnv_data$End)
  )
  
  # Find overlaps between input gene coordinates and CNVs
  overlaps <- findOverlaps(gene_input, cnv_granges)
  
  # Extracting CNV data that overlaps with your gene
  gene_cnv <- cnv_data[subjectHits(overlaps), ]
  
  gene_cnv$Segment_Mean <- as.numeric(gene_cnv$Segment_Mean)
  
  # Creating new column for amplification status
  gene_cnv$CNV_Status <- ifelse(gene_cnv$Segment_Mean > 0.03, "Amplified",
                                ifelse(gene_cnv$Segment_Mean < -0.03, "Deleted", "Neutral"))
  
  # Create the dot plot
  cnv_dplot <- ggplot(gene_cnv, aes(x = Sample, y = Segment_Mean, color = CNV_Status)) +
    geom_jitter(alpha =0.3, width = 0.2, height = 0) + 
    geom_hline(yintercept = 0, linetype = "solid", color = "black") + 
    theme(axis.text.x = element_blank()) + 
    labs(title = paste("Copy Number Variation (CNV) for", input),
         x = "Samples",
         y = "Segment Mean") +
    scale_y_continuous(breaks = seq(-1, 1, by = 0.2))
  
  return(cnv_dplot)
}

cnv_res <- cnv_analysis(input, dds, cnv_data, mart)




#####################################
#########  USER INTERFACE  #########
#####################################

library(shiny)
library(shinythemes)

ui <- fluidPage(
  #title panel
  titlePanel("Multi-omics integration app for gene evaluation in LUSC"),
  
  #sidebar layout for input
  sidebarLayout(
    sidebarPanel(
      #text input to get gene name
      textInput("gene_name", "Enter Gene Symbol:", value = "KRT5"),
      #action button
      actionButton("analyze", "Analyze"),
      
      #dropdown menu to select anaylsis
      selectInput("analysis_type", "Select Analysis Type:",
                  choices = c(
                    "Expression Analysis",
                    "Survival Analysis",
                    "Methylation Status",
                    "Copy Number Variation"
                  ),
                  selected = "Expression Analysis")
    ),
    
    #main panel to display results
    mainPanel(
      tabsetPanel(
        tabPanel("Results",
                 plotOutput("resultplot"),
                 verbatimTextOutput("resultText"))
      ),
      
      #Tab for instructions
      tabPanel("Instructions",
               h4("Instructions for Use"),
               p("1. Enter a gene symbol in the text box."),
               p("2. Select the analysis type from the dropdown."), 
               p("3. Click 'Analyze' to see the results.") 
            )
    )
    )
  )

#server logic
server <- function(input, output, session) {
  #load RDS files
  data <- read_rds("summarized_experiment", refhook = NULL)
  dds <- read_rds("dds_object", refhook = NULL)
  res <- results(dds)
  cnv_data <- read_rds("cnv_data.rds")
  mval <- readRDS("methylation_data.rds")
  
  #getting clinical data
  clinical_data <- as.data.frame(data@colData)
  
  #getting probe annotation
  annotation <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))
  
  #adjusting multivalued cells
  split <- c("UCSC_RefGene_Name", "UCSC_RefGene_Accession", "UCSC_RefGene_Group")
  annotation <- annotation %>%
    separate_rows(!!!syms(split), sep = ";")
  
  #getting gene annotations for CNVs
  mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  
  observeEvent(input$analyze, {
    req(input$gene_name)
    
    if(input$analysis_type == "Expression Analysis") { #if user selects exp analysis then perform DEA and display boxplots
      stats <- get_gene(input$gene_name, dds, res)
      expression_result <- get_expression_data(input$gene_name, dds, count)
      
      output$resultplot <- renderPlot({
        ggplot(expression_result,
               aes(x = sample_type,
                   y = expression,
                   fill = sample_type)) +
          geom_boxplot(alpha =0.5) +
          geom_jitter(shape = 16,
                      position = position_jitter(0.1),
                      color = "black",
                      size = 1) +
          labs(title = paste("Expression of", input, "in Normal vs Tumor samples"),
               x = "Sample Type",
               y = "Normalized Counts") +
          theme_minimal() +
          theme(legend.position = "none")
      })
      output$resulText <- renderText({stats})
    
    } else if (input$analysis_type == "Survival Analysis") { #if surv analysis, then perform survival analysis and display kap-mier plot
        splot <- survival_analysis(input$gene_name, clinical_data, dds, count)
        output$resultplot <- renderPlot(splot)
        output$resultText <- renderText({"Survival Analysis Complete!"})
      
    } else if (input$analysis_type == "Methylation Status") { #if meth status, display mvalue plots
       meth_res <- meth_status(input$gene_name, mval, annotation)
       output$resultplot <- renderPlot(meth_res)
       output$resultText <- renderText({"Methylation Status Displayed!"})
       
    } else if (input$analysis_type == "Copy Number Variation") { #if CNV, then display dotplots
          cnv_res <- cnv_analysis(input$gene_name, dds, cnv_data, mart)
          output$resultplot <- renderPlot(cnv_res)
          output$resultText <- renderText({"Copy Number Variations displayed!"})
        }
  })
  
}

#run the app
shinyApp(ui = ui, server = server)






