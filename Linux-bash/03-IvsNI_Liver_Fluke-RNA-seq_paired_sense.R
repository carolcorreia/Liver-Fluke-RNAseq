##########################################################################
# RNA-seq analysis: Peripheral blood  mononuclear cells and              #
#                     liver fluke in bovine paired-end reads             #
#                      Comparison Infected vs No Infected                #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 3                                 #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version (4.0.0): Garcia-Campos A. 
# DOI badge of current version:
# Last updated on 20/12/2017

##################################
# 30 Working directory and RData #
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
# 31 Load and/or install required packages #
############################################

library(statmod)
library(edgeR)
library(devtools)
library(plyr)
library(tidyverse)
library(stringr)
library(magrittr)
library(purrr)
library(forcats)
library(biobroom)
library(ggjoy)
library(ggrepel)
library(Cairo)
library(cowplot)
library(extrafont)
library(VennDiagram)
library(treemap)

# Uncomment functions below to install packages in case you don't have them

#install.packages("cowplot")
#install.packages("extrafont")
#install.packages("statmod")
install.packages("VennDiagram")
install.packages("treemap")

###################
# 32 Set up fonts #
###################

# Registered fonts with R for the PDF output device
loadfonts()


# For analysis infected versus no infected, all variables are going to the preceded by prefix: "IvsNI_"



#####################################################
# 33 Determine differential expression by fitting a #
# negative binomial GLM with Quasi-likelihood Tests #
#####################################################

# Test for differential expression between the different time points/treatments,
# using the coefficients from IvsNI_Liver_Fluke_fit$design

# 0 wk
IvsNI_W0.QL <- glmQLFTest(IvsNI_Liver_Fluke_fitQL, coef = "Infected")
IvsNI_testDE.W0 <- topTags(object        = IvsNI_W0.QL,
                     n             = "inf",
                     adjust.method = "BH")

head(IvsNI_testDE.W0$table)

# +1 wk
IvsNI_W1.QL <- glmQLFTest(IvsNI_Liver_Fluke_fitQL, coef = "Infected:time.pointW1")
IvsNI_testDE.W1 <- topTags(object        = IvsNI_W1.QL,
                     n             = "inf",
                     adjust.method = "BH")

head(IvsNI_testDE.W1$table)


# +14 wk
IvsNI_W14.QL <- glmQLFTest(IvsNI_Liver_Fluke_fitQL, coef = "Infected:time.pointW14")
IvsNI_testDE.W14 <- topTags(object        = IvsNI_W14.QL,
                      n             = "inf",
                      adjust.method = "BH")

head(IvsNI_testDE.W14$table)

### Output all genes tested for DE

# 0 wk
IvsNI_testDE.W0$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/W0_AllGenes-infected-vs-noinfected.csv")),
            col_names = TRUE)

# +1 wk
IvsNI_testDE.W1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/W1_AllGenes-infected-vs-noinfected.csv")),
            col_names = TRUE)


# +14 wk
IvsNI_testDE.W14$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  dplyr::arrange(FDR) %>%
  write_csv(file.path(paste0(workDir, "/W14_AllGenes-infected-vs-noinfected.csv")),
            col_names = TRUE)


### Filter genes considered DE (FDR < 0.05)

# 0 wk
IvsNI_testDE.W0$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> IvsNI_W0_FDR_05

# +1 wk
IvsNI_testDE.W1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> IvsNI_W1_FDR_05


# +14 wk
IvsNI_testDE.W14$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::filter(FDR < 0.05) %>%
  dplyr::arrange(desc(logFC)) %>%
  as.tibble() -> IvsNI_W14_FDR_05

### Output genes considered DE (FDR < 0.05)
IvsNI_DElists <- list(IvsNI_W0_FDR_05, IvsNI_W1_FDR_05, IvsNI_W14_FDR_05)
IvsNI_DEfiles <- c(paste0(c("IvsNI_W0_FDR_05", "IvsNI_W1_FDR_05", "IvsNI_W14_FDR_05"),
                    "_genes.csv"))
IvsNI_DEpaths <- file.path(tabDir, IvsNI_DEfiles)

pwalk(list(IvsNI_DElists, IvsNI_DEpaths),
      write_csv,
      col_names = TRUE)

###############################################
# 34 Plot: Treemaps of DE genes (FDR < 0.05) #
###############################################

# Get numbers of up and down regulated genes
# at each time point
IvsNI_list_DE <- list(IvsNI_W0_FDR_05, IvsNI_W1_FDR_05, IvsNI_W14_FDR_05)
names(IvsNI_list_DE) <- c("0 wk", "+1 wk", "+14 wk")

IvsNI_Up_Down <- map_df(IvsNI_list_DE,
                  ~ dplyr::count(.x,
                                 up = sum(logFC > 0),
                                 down = sum(logFC < 0),
                                 zero = sum(logFC == 0)),
                  .id = "time_point")

# Time point as factor
IvsNI_Up_Down$time_point %<>%
  factor() %>%
  fct_inorder()

# Plotting labels
IvsNI_Up_Down %<>% dplyr::mutate(labelsUp = paste(time_point, up, sep = ' '))
IvsNI_Up_Down %<>% dplyr::mutate(labelsDown = paste(time_point, down, sep = ' '))

# Check data frame
IvsNI_Up_Down

# Plot chart increased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/tree_up-infected-vs-noinfected.pdf")),
          width    = 8,
          height   = 4,
          family   = "Arial",
          fallback_resolution = 300)
treemap(IvsNI_Up_Down,
        index             = "labelsUp",
        vSize             = "up",
        type              = "index",
        palette           = "PRGn",
        title             = "Increased expression",
        fontsize.title    = 14,
        fontfamily.title  = "Arial",
        fontfamily.labels = "Arial",
        fontsize.labels   = 16)
dev.off()

# Plot chart decreased expression
# Run this chunk together
cairo_pdf(filename = file.path(paste0(imgDir, "/tree_down-infected-vs-noinfected.pdf")),
          width    = 8,
          height   = 4,
          family   = "Arial",
          fallback_resolution = 300)
treemap(IvsNI_Up_Down,
        index             = "labelsDown",
        vSize             = "down",
        type              = "index",
        palette           = "-PRGn",
        title             = "Decreased expression",
        fontsize.title    = 14,
        fontfamily.title  = "Arial",
        fontfamily.labels = "Arial",
        fontsize.labels   = 16)
dev.off()

###################################################
# 35 Plot: Venn diagram of DE genes (FDR < 0.05) #
###################################################

# Turn gene IDs into vectors
IvsNI_W0.vector <- IvsNI_W0_FDR_05$EntrezID
IvsNI_W1.vector <- IvsNI_W1_FDR_05$EntrezID
IvsNI_W14.vector <- IvsNI_W14_FDR_05$EntrezID

# Turn log files off from VennDiagram package
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# Plot Venn diagram for all time points
IvsNI_venn.plot <- venn.diagram(list(A = IvsNI_W0.vector,
                               B = IvsNI_W14.vector,
                               C = IvsNI_W1.vector),
                          filename        = file.path(paste0(imgDir,
                                                             "/Venn_DE_FDR_05.png")),
                          imagetype       = "png",
                          col             = "transparent",
                          fill            = c("#ffffcc",
                                              "#225ea8",
                                              "#a1dab4"),
                          alpha           = 0.50,
                          label.col       = "#003333",
                          cex             = 3,
                          fontfamily      = "Arial",
                          category.names  = c("0 wk",
                                              "+14 wk",
                                              "+1 wk"),
                          cat.col         = "black",
                          cat.cex         = 3,
                          cat.pos         = c(-11, 11, 0),
                          cat.dist        = c(0.21, 0.21, 0.1),
                          cat.fontfamily  = "Arial",
                          rotation.degree = 360,
                          margin          = 0,
                          height          = 11,
                          width           = 11,
                          units           = 'in',
                          compression     = 'lzw',
                          resolution      = 300)



#######################################END  ############################
######################
# 36 Common DE genes #
######################

# Join common DE genes to all time points into single data frame
IvsNI_W0_FDR_05 %>%
  dplyr::inner_join(IvsNI_W1_FDR_05,
                    by = "EntrezID",
                    suffix = c("_W0", "_W1")) %>%
  dplyr::inner_join(IvsNI_W1_FDR_05,
                    by = "EntrezID") %>%
  dplyr::inner_join(IvsNI_W14_FDR_05,
                    by = "EntrezID",
                    suffix = c("_W1", "_W14")) -> IvsNI_common_DE

# Join common DE genes to infected samples into single data frame    I have not done this
IvsNI_W1_FDR_05 %>%
  dplyr::inner_join(IvsNI_W14_FDR_05,
                    by = "EntrezID",
                    suffix = c("_W1", "_W14")) %>%
  dplyr::anti_join(IvsNI_W0_FDR_05,
                   by = "EntrezID") -> IvsNI_common_DE_infected

# Select DE genes unique to healthy samples
IvsNI_W0_FDR_05 %>%
  dplyr::anti_join(IvsNI_W1_FDR_05,
                    by = "EntrezID") %>%
  dplyr::anti_join(IvsNI_W14_FDR_05,
                    by = "EntrezID") -> IvsNI_healthy_DE

# Check data frames
IvsNI_common_DE
IvsNI_common_DE_infected
IvsNI_healthy_DE

# Output data
IvsNI_VennDE <- list(IvsNI_common_DE, IvsNI_common_DE_infected, IvsNI_healthy_DE)
IvsNI_VennDEfiles <- c(paste0(c("IvsNI_common_DE", "IvsNI_common_DE_infected", "IvsNI_healthy_DE"),
                        "_genes.csv"))
IvsNI_VennDEpaths <- file.path(tabDir, IvsNI_VennDEfiles)

pwalk(list(IvsNI_VennDE, IvsNI_VennDEpaths),
      write_csv,
      col_names = TRUE)

####################################################
# 37 Plot: Volcano of DE genes at each time point  #
####################################################

# 0 wk
IvsNI_testDE.W0$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 3 & FDR < 0.05,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Arial",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-4, 4)) +
  ylim(c(0, 10)) +
  theme_bw(base_size = 14, base_family = "Arial") +
  ggtitle("PRE") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.05",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> IvsNI_W0_Volcano

IvsNI_W0_Volcano

# +1 wk
IvsNI_testDE.W1$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 3 & FDR < 0.05,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Arial",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-4, 4)) +
  ylim(c(0, 10)) +
  theme_bw(base_size = 14, base_family = "Arial") +
  ggtitle("1wpi") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.05",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> IvsNI_W1_Volcano

IvsNI_W1_Volcano

# +14 wk
IvsNI_testDE.W14$table %>%
  rownames_to_column(var = "EntrezID") %>%
  dplyr::mutate(is.DE = if_else(FDR < 0.05, "true", "false")) %>%
  ggplot(aes(x = logFC,
             y = -log10(FDR),
             colour = is.DE)) +
  geom_point(alpha = 0.4, size  = 1.75) +
  geom_label_repel(aes(label = if_else(abs(logFC) > 5 & FDR < 0.05,
                                       gene_symbol, NULL),
                       fontface = "italic",
                       family = "Arial",
                       size = 4,
                       fill = "gray",
                       alpha = 0.5),
                   show.legend = FALSE) +
  xlim(c(-4, 4)) +
  ylim(c(0, 10)) +
  theme_bw(base_size = 14, base_family = "Arial") +
  ggtitle("14wpi") +
  xlab(expression(paste(log[2], " fold-change"))) +
  ylab(expression(paste(-log[10], " FDR"))) +
  scale_colour_manual("FDR \n< 0.05",
                      labels = c("False", "True"),
                      values = c('#41b6c4', '#225ea8')) -> IvsNI_W14_Volcano

IvsNI_W14_Volcano


### Export high quality image for all volcano plots
IvsNI_files_V <- paste0(c("IvsNI_W0_Volcano", "IvsNI_W1_Volcano",
                    "IvsNI_W14_Volcano"), ".pdf")
IvsNI_plots_V <- list(IvsNI_W0_Volcano, IvsNI_W1_Volcano, IvsNI_W14_Volcano)

purrr::pwalk(list(IvsNI_files_V, IvsNI_plots_V),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 12,
             width     = 14,
             units     = "in")

####################################################
# 38 Plot: Combine all volcanos into single figure #
####################################################

# Set grids
IvsNI_V_grid <- plot_grid(IvsNI_W0_Volcano,
                           IvsNI_W1_Volcano, IvsNI_W14_Volcano,
                     labels = c("A", "B", "C"),
                     nrow = 3)

# Check plots
IvsNI_V_grid

# Export high quality image for both grids
IvsNI_files_Vgrid <- paste0(c("IvsNI_V_grid"), ".pdf")
IvsNI_plots_Vgrid <- list(IvsNI_V_grid)

purrr::pwalk(list(IvsNI_files_Vgrid, IvsNI_plots_Vgrid),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 15,
             width     = 10,
             units     = "in")

#######################
# 39 Save .RData file #
#######################

save.image(file = "Liver_Fluke_BovinePBMC.RData")

##########################
# 40 Save R session info #
##########################

devtools::session_info()

############################
# Proceed to RNA-seq stats #
############################