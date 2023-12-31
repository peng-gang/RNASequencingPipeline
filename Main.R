library(ggplot2)
library(ggrepel)
library(DESeq2)
library(ComplexHeatmap)
library(ggtext)
library(RColorBrewer)
library(org.Mm.eg.db)
library(enrichplot)
library(ggsci)
library(clusterProfiler)
library(ggupset)
library(circlize)
library(RColorBrewer)
source("Analysis_DGE.R")
source("Quality_Control.R")

# Input variable name in sampleInfo as Variable_Of_Interest
# Select two variables in Variable_of_Interest as Groups_Selected, Groups_selected[1] as reference
# Clarify the sample name in sampleInfo as Samples_Column_Name
Variable_Of_Interest = "group"
Groups_Selected = c("KO_2","KO_8")
Samples_Column_Name = "ID"

# Filter the groups required by analysis 
idx <- sampleInfo[[Variable_Of_Interest]] %in% Groups_Selected
sampleInfoSel <- sampleInfo[idx,]
sampleInfoSel[[Variable_Of_Interest]] <- factor(sampleInfoSel[[Variable_Of_Interest]], levels = Groups_Selected)

Folder_Name <- "KO_2_vs_KO_8"

# function to calculate the normalize counts for all the samples
dds <- Normalize_Counts(cts = cts, sampleInfoSel = sampleInfoSel)

# function to plot the heatmap of sample distance
Distance_Clustering(dds = dds)

# function to perform the principle component analysis
PCA_Plots(dds = dds, sampleInfoSel = sampleInfoSel, 
          Variables_For_PCA = c("group"),
          Color_Choice = c("Blue","black"), Shape_Choice = c(12, 13, 14, 15))

# function to generate the plots for total and proportion reads on Y chromosome
Y_Reads(cts = cts, geneInfo = geneInfo, sampleInfo = sampleInfo, Folder_Name = "qcdata", Chromosome_CN = "Chr", gender_column_name = "sex", Samples_column_name = "id")

# function to generate the plots for total and proportion reads on X chromosome
X_Reads(cts = cts, geneInfo = geneInfo, sampleInfo = sampleInfo, Folder_Name = "qcdata", Chromosome_CN = "Chr", gender_column_name = "sex", Samples_column_name = "id")

# function to generate the plots for total and proportion reads on XIST gene
XIST_Counts(cts = cts, geneInfo = geneInfo, genes_column_name = "SYMBOL", sampleInfo = sampleInfo, Folder_Name = "qcdata", gender_column_name = "sex", Samples_column_name = "id")

# geneInfo giving names of input files containing read mapping results of each gene
# function to perform DESeq2 analysis using a categorical variable between two selected groups
twoGroupCompare(Feature_Counts = cts, Sample_Info = sampleInfo, 
                Gene_Info = geneInfo, Genes_Column_Name = "ENSEMBL",
                Covariates = NULL,
                pvalue_Cutoff = 0.05, log2Fold_Cutoff = log2(1.5))

# Uncomment this section if you want to check for age-dependent effect
# 
# # Addition of fictitious age data to sampleInfo
# Sample_Info <- data.frame(sampleInfo)  
# 
# # Set the number of rows in the dataframe
# num_rows <- nrow(Sample_Info)
# 
# # Set the range for the random age values
# min_age <- 18
# max_age <- 65
# 
# # Generate a vector of random age values
# age_values <- sample(min_age:max_age, num_rows, replace = TRUE)
# 
# # Add the "age" column to the dataframe
# Sample_Info$age <- age_values


# function to perform DESeq2 analysis using a continuous variable
# continuousCompare(Feature_Counts = cts, Sample_Info = Sample_Info, Samples_Column_Name = "id", Gene_Info = geneInfo, Genes_Column_Name = "SYMBOL", Variable_Of_Interest = "age", Folder_name = "age_compare", pvalue_Cutoff = 0.05)

