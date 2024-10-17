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

#getting clinical data
clinical_data <- as.data.frame(data@colData)

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





