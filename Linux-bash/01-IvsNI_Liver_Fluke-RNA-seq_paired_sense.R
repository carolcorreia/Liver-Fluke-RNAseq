##########################################################################
# RNA-seq analysis: Peripheral blood  mononuclear cells and              #
#                     liver fluke in bovine paired-end reads             #
#                      Comparison Infected vs No Infected                #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 1                                 #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version: Andres Garcia-Campos
# DOI badge of current version:
# Last updated on 20/12/2017

# For analysis infected versus no infected, all variables are going to the preceded by prefix: "IvsNI_" see from point 3

##################################
# 01 Working directory and RData #
##################################

# Set working directory
setwd("/Users/andresgarciacampos/Dropbox/Liver-Fluke-RNAseq/R_analysis")
getwd()

# Define variables for working and file directories
workDir <- getwd()
workDir
fileDir <- "/Users/andresgarciacampos/Dropbox/Liver-Fluke-RNAseq/R_analysis/edgeR_sense/Counts"
imgDir <- paste0(workDir, "/Figures")
tabDir <- paste0(workDir, "/Tables")

# Load previously saved data
load("Liver_Fluke_BovinePBMC.RData")

############################################
# 02 Load and/or install required packages #
############################################

library(Biobase)
library(edgeR)
library(AnnotationFuncs)
library(org.Bt.eg.db)
library(devtools)
library(plyr)
library(tidyverse)
library(stringr)
library(magrittr)
library(biobroom)
library(ggridges)
library(ggrepel)
library(extrafont)

# Uncomment functions below to install packages in case you don't have them

# Bioconductor packages
source("https://bioconductor.org/biocLite.R")
biocLite("AnnotationFuncs")
biocLite("DESeq2") #I can delete it because I am not going to use it
biocLite("edgeR")
biocLite("Biobase")
biocLite("geneLenDataBase")
biocLite("org.Bt.eg.db")

# CRAN packages
install.packages("plyr")
install.packages("tidyverse")
install.packages("ggrepel")
install.packages("ggridges")
install.packages("devtools") # I added it because I could not do library devtools (But it is in tidyverse but its not required)
install.packages("extrafont")

# Registered fonts with R for the PDF output device
loadfonts()

#########################################
# 03 Import featureCounts sense counts  #     
#########################################

# Create a vector of file names
IvsNI_files <- list.files (path       = "/Users/andresgarciacampos/Dropbox/Liver-Fluke-RNAseq/Counts",
                     pattern     = "*S*_sense*",
                     all.files   = TRUE,
                     full.names  = FALSE,
                     recursive   = FALSE,
                     ignore.case = FALSE)

IvsNI_files
length(IvsNI_files)

# Create a dataframe with raw counts for all samples
IvsNI_rawCounts <- readDGE(path         = "/Users/andresgarciacampos/Dropbox/Liver-Fluke-RNAseq/Counts",
                     files        = IvsNI_files,
                     header       = TRUE,
                     comment.char = "#",
                     columns      = c(1, 7))
names(IvsNI_rawCounts)
head(IvsNI_rawCounts$samples)
head(IvsNI_rawCounts$counts)

#################################
# 04 Clean column and row names #
#################################

# Correct column names in counts
colnames(IvsNI_rawCounts$counts) %<>%
  str_replace("_S.*", "")

# Correct row names in counts
rownames(IvsNI_rawCounts$counts) %<>%
  str_replace("BGD.*,", "") %>%
  str_replace(",miRBase:.*", "") %>%
  str_replace("GeneID:", "")

# Correct row names in samples
rownames(IvsNI_rawCounts$samples) %<>%
  str_replace("_S.*", "")

# Check data frames
head(IvsNI_rawCounts$counts)
head(IvsNI_rawCounts$samples)

#################################
# 05 Get gene names and symbols #
#################################

# Create annotation table with counts information
IvsNI_annotCounts <- as.data.frame(IvsNI_rawCounts$counts)
columns(org.Bt.eg.db)

# Retrieve gene names using NCBI Ref Seq gene identifiers
IvsNI_annotCounts$gene_name <- mapIds(org.Bt.eg.db,
                                keys      = rownames(IvsNI_annotCounts),
                                column    = "GENENAME",
                                keytype   = "ENTREZID",
                                multiVals = "first")

# Retrieve gene symbols using NCBI Ref Seq gene identifiers
IvsNI_annotCounts$gene_symbol <- mapIds(org.Bt.eg.db,
                                  keys      = rownames(IvsNI_annotCounts),
                                  column    = "SYMBOL",
                                  keytype   = "ENTREZID",
                                  multiVals = "first")

head(IvsNI_annotCounts)
dim(IvsNI_annotCounts)

# Ouptut data
IvsNI_annotCounts %>%
  dplyr::select(gene_symbol, gene_name, everything()) %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(file.path(paste0(workDir, "/Liver_fluke-bovine_PBMC-sense_annot-rawcounts_infected-noinfected.csv")),
            col_names = TRUE)

#########################################
# 06 Add sample information for DGElist #
#########################################

# Treatment group (avoid using underscores)
IvsNI_rawCounts$samples$group <- rownames(IvsNI_rawCounts$samples)
IvsNI_rawCounts$samples$group %<>%
  str_replace("P14", "Infected") %>%
  str_replace("P05", "Infected") %>%
  str_replace("P01", "Infected") %>%
  str_replace("P21", "Infected") %>%
  str_replace("P16", "Infected") %>%
  str_replace("P04", "Infected") %>%
  str_replace("P12", "Infected") %>%
  str_replace("P19", "Infected") %>%
  str_replace("P02", "Infected") %>%
  str_replace("P15", "Infected") %>%
  str_replace("P17", "Infected") %>%
  str_replace("A07", "Infected") %>%
  str_replace("A10", "Infected") %>%
  str_replace("A02", "Infected") %>%
  str_replace("A21", "Infected") %>%
  str_replace("A12", "Infected") %>%
  str_replace("A08", "Infected") %>%
  str_replace("A11", "Infected") %>%
  str_replace("A17", "Infected") %>%
  str_replace("A01", "Infected") %>%
  str_replace("A16", "Infected") %>%
  str_replace("A19", "Infected") %>%
  str_replace("C03", "Infected") %>%
  str_replace("C16", "Infected") %>%
  str_replace("C01", "Infected") %>%
  str_replace("C17", "Infected") %>%
  str_replace("C13", "Infected") %>%
  str_replace("C08", "Infected") %>%
  str_replace("C15", "Infected") %>%
  str_replace("C22", "Infected") %>%
  str_replace("C02", "Infected") %>%
  str_replace("C12", "Infected") %>%
  str_replace("C18", "Infected") %>%
  str_replace("P13", "Control") %>%
  str_replace("P07", "Control") %>%
  str_replace("P09", "Control") %>%
  str_replace("P10", "Control") %>%
  str_replace("P06", "Control") %>%
  str_replace("P22", "Control") %>%
  str_replace("P11", "Control") %>%
  str_replace("P18", "Control") %>%
  str_replace("P08", "Control") %>%
  str_replace("A13", "Control") %>%
  str_replace("A06", "Control") %>%
  str_replace("A14", "Control") %>%
  str_replace("A09", "Control") %>%
  str_replace("A03", "Control") %>%
  str_replace("A22", "Control") %>%
  str_replace("A15", "Control") %>%
  str_replace("A18", "Control") %>%
  str_replace("A05", "Control") %>%
  str_replace("C10", "Control") %>%
  str_replace("C07", "Control") %>%
  str_replace("C14", "Control") %>%
  str_replace("C11", "Control") %>%
  str_replace("C04", "Control") %>%
  str_replace("C21", "Control") %>%
  str_replace("C09", "Control") %>%
  str_replace("C20", "Control") %>%
  str_replace("C05", "Control") %>%
  factor(levels = c("Control", "Infected"))

# Animal ID (cannot start with numbers)
IvsNI_rawCounts$samples$animal <- rownames(IvsNI_rawCounts$samples)
IvsNI_rawCounts$samples$animal %<>%
  str_replace("P05", "A2552") %>%
  str_replace("P14", "A2503") %>%
  str_replace("P01", "A2505") %>%
  str_replace("P21", "A2512") %>%
  str_replace("P16", "A2513") %>%
  str_replace("P04", "A2519") %>%
  str_replace("P12", "A2522") %>%
  str_replace("P19", "A2531") %>%
  str_replace("P02", "A2532") %>%
  str_replace("P15", "A2540") %>%
  str_replace("P17", "A2547") %>%
  str_replace("A07", "A2552") %>%
  str_replace("A10", "A2503") %>%
  str_replace("A02", "A2505") %>%
  str_replace("A21", "A2512") %>%
  str_replace("A12", "A2513") %>%
  str_replace("A08", "A2519") %>%
  str_replace("A11", "A2522") %>%
  str_replace("A17", "A2531") %>%
  str_replace("A01", "A2532") %>%
  str_replace("A16", "A2540") %>%
  str_replace("A19", "A2547") %>%
  str_replace("C03", "A2552") %>%
  str_replace("C16", "A2503") %>%
  str_replace("C01", "A2505") %>%
  str_replace("C17", "A2512") %>%
  str_replace("C13", "A2513") %>%
  str_replace("C08", "A2519") %>%
  str_replace("C15", "A2522") %>%
  str_replace("C22", "A2531") %>%
  str_replace("C02", "A2532") %>%
  str_replace("C12", "A2540") %>%
  str_replace("C18", "A2547") %>%
  str_replace("P13", "A2502") %>%
  str_replace("P07", "A2507") %>%
  str_replace("P09", "A2511") %>%
  str_replace("P10", "A2520") %>%
  str_replace("P06", "A2525") %>%
  str_replace("P22", "A2526") %>%
  str_replace("P11", "A2539") %>%
  str_replace("P18", "A2544") %>%
  str_replace("P08", "A2553") %>%
  str_replace("A13", "A2502") %>%
  str_replace("A06", "A2507") %>%
  str_replace("A14", "A2511") %>%
  str_replace("A09", "A2520") %>%
  str_replace("A03", "A2525") %>%
  str_replace("A22", "A2526") %>%
  str_replace("A15", "A2539") %>%
  str_replace("A18", "A2544") %>%
  str_replace("A05", "A2553") %>%
  str_replace("C10", "A2502") %>%
  str_replace("C07", "A2507") %>%
  str_replace("C14", "A2511") %>%
  str_replace("C11", "A2520") %>%
  str_replace("C04", "A2525") %>%
  str_replace("C21", "A2526") %>%
  str_replace("C09", "A2539") %>%
  str_replace("C20", "A2544") %>%
  str_replace("C05", "A2553") %>%
  factor()

# Time points (avoid using underscores)
IvsNI_rawCounts$samples$time.point <- rownames(IvsNI_rawCounts$samples)
IvsNI_rawCounts$samples$time.point %<>%
  str_replace("P\\d\\d", "PRE") %>%
  str_replace("A\\d\\d", "W1") %>%
  str_replace("C\\d\\d", "W14") %>%
  factor(levels = c("PRE", "W1", "W14"))

# Check data frame
head(IvsNI_rawCounts$samples)

#####################
# 07 Create DGElist #
#####################

# Assign required information to variables
IvsNI_gene_annotation <- dplyr::select(IvsNI_annotCounts,
                                 gene_name,
                                 gene_symbol)

IvsNI_raw_counts <- as.data.frame(IvsNI_rawCounts$counts)

IvsNI_group <- IvsNI_rawCounts$samples$group

# Use newly assigned variables to create DGElist
IvsNI_Liver_Fluke_dgelist <- DGEList(counts       = IvsNI_raw_counts,
                        group        = IvsNI_group,
                        genes        = IvsNI_gene_annotation,
                        lib.size     = NULL,
                        norm.factors = NULL,
                        remove.zeros = FALSE)

names(IvsNI_Liver_Fluke_dgelist)
dim(IvsNI_Liver_Fluke_dgelist)
head(IvsNI_Liver_Fluke_dgelist$counts)
head(IvsNI_Liver_Fluke_dgelist$samples)
head(IvsNI_Liver_Fluke_dgelist$genes)

# Include addtional experimental information into DGElist
identical(rownames(IvsNI_rawCounts$samples), rownames(IvsNI_Liver_Fluke_dgelist$samples))

IvsNI_Liver_Fluke_dgelist$samples$animal <- IvsNI_rawCounts$samples$animal
IvsNI_Liver_Fluke_dgelist$samples$time.point <- IvsNI_rawCounts$samples$time.point

head(IvsNI_Liver_Fluke_dgelist$samples)

################################################
# 08 Density plot: raw gene counts per library #
################################################

# Tidy DGElist and plot data
IvsNI_Liver_Fluke_dgelist %>%
  tidy() %>%
  ggplot() +
      geom_density(aes(x     = log10(count + 1),
                       group = sample)) +
      theme_bw(base_size = 12, base_family = "Arial") +
      ylab("Density of raw gene counts per sample") +
      xlab(expression(paste(log[10], "(counts + 1)"))) -> IvsNI_density_raw


IvsNI_density_raw

# Export image
ggsave("Liver_Fluke_Bovine-density_plot_raw_counts_infected_vs_noinfected.png",
       plot      = IvsNI_density_raw,
       device    = "png",
       limitsize = FALSE,
       dpi       = 300,
       path      = workDir)

###########################################
# 09 Remove zero and lowly expressed tags #
###########################################

# Filter non expressed tags (all genes that have zero counts in all samples)
IvsNI_Liver_Fluke_no_zeros <- IvsNI_Liver_Fluke_dgelist[rowSums(IvsNI_Liver_Fluke_dgelist$counts) > 0, ]
dim(IvsNI_Liver_Fluke_no_zeros$counts)
head(IvsNI_Liver_Fluke_no_zeros$counts)
colnames(IvsNI_Liver_Fluke_no_zeros$counts)

# Filter lowly expressed tags, retaining only tags with
# more than 1 count per million in 10 or more libraries
# (9 libraries correspond to 9 biological replicates and represent
# the smallest group (control group)at any given time point)
IvsNI_Liver_Fluke_filt <- IvsNI_Liver_Fluke_no_zeros[rowSums(cpm(IvsNI_Liver_Fluke_no_zeros) > 1) >= 9, ] #(9 because is the smallest group)
dim(IvsNI_Liver_Fluke_filt$counts)
head(IvsNI_Liver_Fluke_filt$counts)

# Ouptut filtered counts
IvsNI_Liver_Fluke_filt$counts %>%
  as.data.frame() %>%
  rownames_to_column(var = "EntrezID") %>%
  write_csv(file.path(paste0(workDir, "/Liver_Fluke-sense_filt_counts-infected_vs_noinfected.csv")),
            col_names = TRUE)

##############################
# 10 Recompute library sizes #
##############################

IvsNI_Liver_Fluke_filt$samples$lib.size <- colSums(IvsNI_Liver_Fluke_filt$counts)
head(IvsNI_Liver_Fluke_filt$samples)
head(IvsNI_Liver_Fluke_dgelist$samples)

###########################################################################
# 11 Calculate normalisation factors using Trimmed Mean of M-values (TMM) #
###########################################################################

# With edgeR, counts are not transformed in any way after
# calculating normalisation factors
IvsNI_Liver_Fluke_norm <- calcNormFactors(IvsNI_Liver_Fluke_filt, method = "TMM")
head(IvsNI_Liver_Fluke_norm$samples)


#######################
# 12 Save .RData file # ----
#######################

save.image(file = "Liver_Fluke_BovinePBMC.RData")

##########################
# 13 Save R session info # ----
##########################

devtools::session_info()

######################################
# Proceed to Part 2 of this analysis #
######################################

# File: 02Liver_Fluke-RNA-seq_paired_sense.R