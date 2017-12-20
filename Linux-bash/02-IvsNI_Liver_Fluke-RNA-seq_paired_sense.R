##########################################################################
# RNA-seq analysis: Peripheral blood  mononuclear cells and              #
#                     liver fluke in bovine paired-end reads             #
#                      Comparison Infected vs No Infected                #
#                                                                        #
#           --- R workflow for analyses of known sense genes ---         #
#                                 Part 2                                 #
##########################################################################

# Based on the workflow created by Nalpas, N.C. (2014)
# DOI badge: http://dx.doi.org/10.5281/zenodo.12474

# Author of current version (4.0.0): Garcia-Campos A.
# DOI badge of current version:
# Last updated on 20/12/2017

# For analysis infected versus no infected, all variables are going to the preceded by prefix: "IvsNI_"

##################################
# 14 Working directory and RData #
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
# 15 Load and/or install required packages #
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
library(ggridges)

# Uncomment functions below to install packages in case you don't have them

install.packages("cowplot")
#install.packages("extrafont")
install.packages("statmod")
install.packages("Cairo")

###################
# 16 Set up fonts #
###################

# Registered fonts with R for the PDF output device
loadfonts()

##########################################################
# 17 Tidy DGElist for exploratory data analysis plotting #
##########################################################

IvsNI_tidy_Liver_Fluke_norm <-tidy(IvsNI_Liver_Fluke_norm, addSamples = TRUE)

# Check factors
levels(IvsNI_tidy_Liver_Fluke_norm$group)

# Clean animal IDs
IvsNI_tidy_Liver_Fluke_norm$animal %<>%
  stringr::str_replace("A", "") %>%
  fct_inorder()

# Check factors
levels(IvsNI_tidy_Liver_Fluke_norm$animal)

# Correct time point info
IvsNI_tidy_Liver_Fluke_norm$time.point %<>%
  factor(levels = c("PRE", "W1", "W14"))

# Check factors
levels(IvsNI_tidy_Liver_Fluke_norm$time.point)

# Combine animal and time point info for
# plotting labels
IvsNI_tidy_Liver_Fluke_norm %<>%
  dplyr::mutate(labels = paste0(time.point, "_", animal))

IvsNI_tidy_Liver_Fluke_norm$labels %<>%
  factor(levels = c("PRE_2532", "PRE_2505", "PRE_2525", "PRE_2553", "PRE_2507", "PRE_2519",
                    "PRE_2520", "PRE_2503", "PRE_2522", "PRE_2513", "PRE_2502", "PRE_2511",
                    "PRE_2539", "PRE_2540", "PRE_2531", "PRE_2544", "PRE_2547", "PRE_2512","PRE_2526","PRE_2552",
                    "W1_2532", "W1_2505", "W1_2525", "W1_2553", "W1_2507", "W1_2519",
                    "W1_2520", "W1_2503", "W1_2522", "W1_2513", "W1_2502", "W1_2511","W1_2552",
                    "W1_2539", "W1_2540", "W1_2531", "W1_2544", "W1_2547", "W1_2512","W1_2526",
                    "W14_2532", "W14_2505", "W14_2525", "W14_2553", "W14_2507", "W14_2519",
                    "W14_2520", "W14_2503", "W14_2522", "W14_2513", "W14_2502", "W14_2511","W14_2552",
                    "W14_2539", "W14_2540", "W14_2531", "W14_2544", "W14_2547", "W14_2512","W14_2526"))

# Check factors
levels(IvsNI_tidy_Liver_Fluke_norm$labels)

# Check data frame
IvsNI_tidy_Liver_Fluke_norm

########################################################
# 18 Plot: density of filtered gene counts per library #
########################################################

ggplot(IvsNI_tidy_Liver_Fluke_norm, aes(x = log10(count + 1),
                          y = labels)) +
  scale_y_discrete(limits = rev(levels(IvsNI_tidy_Liver_Fluke_norm$labels))) +
  geom_density_ridges(aes(fill = group), alpha = 0.5) +
  scale_fill_manual("Treatment", values = c("#b2b2b2", "#e06377")) +
  theme_bw(base_size = 12, base_family = "Arial") +
  ggtitle("Density of filtered gene counts per sample") +
  facet_grid(. ~ group, scales = "free") +
  ylab("Time point_Animal number") +
  xlab(expression(paste(log[10], "(counts + 1)"))) -> IvsNI_density_norm

IvsNI_density_norm

# Export high quality image       
ggsave("Liver_Fluke-density-filt-infected-vs-noinfected.pdf",
       plot      = IvsNI_density_norm,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 9,
       units     = "in")

####################################
# 19 Plots: MDS at each time point #
####################################

### Plot MDS of Pre-infection (PRE)
IvsNI_mds_PRE <- plotMDS.DGEList(IvsNI_Liver_Fluke_norm [, grep(pattern = "P\\d\\d",
                                                 x = colnames(IvsNI_Liver_Fluke_norm))],
                           plot = FALSE,
                           method = "bcv")

names(IvsNI_mds_PRE)
IvsNI_PRE_coord <- IvsNI_mds_PRE$cmdscale.out # Get coords to plot with ggplot2
# Tidy coords
IvsNI_PRE_coord %<>% 
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)
# Clean animal IDs for plotting
IvsNI_PRE_coord$animal %<>% 
  str_replace("P14", "A2503") %>%
  str_replace("P05", "A2552") %>%
  str_replace("P01", "A2505") %>%
  str_replace("P21", "A2512") %>%
  str_replace("P16", "A2513") %>%
  str_replace("P04", "A2519") %>%
  str_replace("P12", "A2522") %>%
  str_replace("P19", "A2531") %>%
  str_replace("P02", "A2532") %>%
  str_replace("P15", "A2540") %>%
  str_replace("P17", "A2547") %>%
  str_replace("P07", "A2507") %>%
  str_replace("P09", "A2511") %>%
  str_replace("P10", "A2520") %>%
  str_replace("P06", "A2525") %>%
  str_replace("P22", "A2526") %>%
  str_replace("P11", "A2539") %>%
  str_replace("P18", "A2544") %>%
  str_replace("P08", "A2553") %>%
  str_replace("P13", "A2502") %>%
  str_replace("A", "") %>%
  factor()
# Clean group info for plotting
IvsNI_PRE_coord$group %<>% 
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
  str_replace("P13", "Control") %>%
  str_replace("P07", "Control") %>%
  str_replace("P09", "Control") %>%
  str_replace("P10", "Control") %>%
  str_replace("P06", "Control") %>%
  str_replace("P22", "Control") %>%
  str_replace("P11", "Control") %>%
  str_replace("P18", "Control") %>%
  str_replace("P08", "Control") %>%
  factor(levels = c("Control", "Infected"))

# Check tidy coords
IvsNI_PRE_coord

# Plot MDS with ggplot2
IvsNI_MDS_PRE <- ggplot(IvsNI_PRE_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  expand_limits(x=c(-0.6, 0.6), y=c(-0.6, 0.6)) +
  scale_colour_manual("Treatment",
                      values = c("#e06377", "#b2b2b2")) +
  scale_shape_manual("Treatment",
                     values = c(17, 19)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Arial") +
  ggtitle("PRE INFECTION") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")


# Check PRE MDS plot
IvsNI_MDS_PRE


### Plot MDS of Week +1 time point
IvsNI_mds_W1 <- plotMDS.DGEList(x = IvsNI_Liver_Fluke_norm[, grep(pattern = "A\\d\\d",
                                               x = colnames(IvsNI_Liver_Fluke_norm))],
                          plot = FALSE,
                          method = "bcv")

names(IvsNI_mds_W1)
IvsNI_W1_coord <- IvsNI_mds_W1$cmdscale.out # Get coords to plot with ggplot2
# Tidy coords
IvsNI_W1_coord %<>% 
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)
# Clean animal IDs for plotting
IvsNI_W1_coord$animal %<>% 
  str_replace("A10", "A2503") %>%
  str_replace("A07", "A2552") %>%
  str_replace("A02", "A2505") %>%
  str_replace("A21", "A2512") %>%
  str_replace("A12", "A2513") %>%
  str_replace("A08", "A2519") %>%
  str_replace("A11", "A2522") %>%
  str_replace("A17", "A2531") %>%
  str_replace("A01", "A2532") %>%
  str_replace("A16", "A2540") %>%
  str_replace("A19", "A2547") %>%
  str_replace("A13", "A2502") %>%
  str_replace("A06", "A2507") %>%
  str_replace("A14", "A2511") %>%
  str_replace("A09", "A2520") %>%
  str_replace("A03", "A2525") %>%
  str_replace("A22", "A2526") %>%
  str_replace("A15", "A2539") %>%
  str_replace("A18", "A2544") %>%
  str_replace("A05", "A2553") %>%
  str_replace("A", "") %>%
  factor()

IvsNI_W1_coord$group %<>% # Clean group info for plotting
  str_replace("A10", "Infected") %>%
  str_replace("A07", "Infected") %>%
  str_replace("A02", "Infected") %>%
  str_replace("A21", "Infected") %>%
  str_replace("A12", "Infected") %>%
  str_replace("A08", "Infected") %>%
  str_replace("A11", "Infected") %>%
  str_replace("A17", "Infected") %>%
  str_replace("A01", "Infected") %>%
  str_replace("A16", "Infected") %>%
  str_replace("A19", "Infected") %>%
  str_replace("A13", "Control") %>%
  str_replace("A06", "Control") %>%
  str_replace("A14", "Control") %>%
  str_replace("A09", "Control") %>%
  str_replace("A03", "Control") %>%
  str_replace("A22", "Control") %>%
  str_replace("A15", "Control") %>%
  str_replace("A18", "Control") %>%
  str_replace("A05", "Control") %>%
  factor(levels = c("Control", "Infected"))

# Check tidy coords
IvsNI_W1_coord

# Plot MDS with ggplot2
IvsNI_MDS_W1 <- ggplot(IvsNI_W1_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  expand_limits(x=c(-0.6, 0.6), y=c(-0.6, 0.6)) +
  scale_colour_manual("Treatment",
                      values = c("#e06377", "#b2b2b2")) +
  scale_shape_manual("Treatment",
                     values = c(17, 19)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Arial") +
  ggtitle("1 WPI") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Check W1 MDS plot
IvsNI_MDS_W1

### Plot MDS of Week +14 time point
IvsNI_mds_W14 <- plotMDS.DGEList(x = IvsNI_Liver_Fluke_norm[, grep(pattern = "C\\d\\d",
                                        x = colnames(IvsNI_Liver_Fluke_norm))],
                           plot = FALSE,
                           method = "bcv")

names(IvsNI_mds_W14)
IvsNI_W14_coord <- IvsNI_mds_W14$cmdscale.out # Get coords to plot with ggplot2
# Tidy coords
IvsNI_W14_coord %<>% 
  tidy() %>%
  dplyr::rename(sample = .rownames, x = X1, y = X2) %>%
  dplyr::mutate(animal = sample, group = sample)
# Clean animal IDs for plotting
IvsNI_W14_coord$animal %<>% 
  str_replace("C16", "A2503") %>%
  str_replace("C03", "A2552") %>%
  str_replace("C01", "A2505") %>%
  str_replace("C17", "A2512") %>%
  str_replace("C13", "A2513") %>%
  str_replace("C08", "A2519") %>%
  str_replace("C15", "A2522") %>%
  str_replace("C22", "A2531") %>%
  str_replace("C02", "A2532") %>%
  str_replace("C12", "A2540") %>%
  str_replace("C18", "A2547") %>%
  str_replace("C10", "A2502") %>%
  str_replace("C07", "A2507") %>%
  str_replace("C14", "A2511") %>%
  str_replace("C11", "A2520") %>%
  str_replace("C04", "A2525") %>%
  str_replace("C21", "A2526") %>%
  str_replace("C09", "A2539") %>%
  str_replace("C20", "A2544") %>%
  str_replace("C05", "A2553") %>%
  str_replace("A", "") %>%
  factor()

IvsNI_W14_coord$group %<>% # Clean group info for plotting
  str_replace("C16", "Infected") %>%
  str_replace("C03", "Infected") %>%
  str_replace("C01", "Infected") %>%
  str_replace("C17", "Infected") %>%
  str_replace("C13", "Infected") %>%
  str_replace("C08", "Infected") %>%
  str_replace("C15", "Infected") %>%
  str_replace("C22", "Infected") %>%
  str_replace("C02", "Infected") %>%
  str_replace("C12", "Infected") %>%
  str_replace("C18", "Infected") %>%
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

# Check tidy coords
IvsNI_W14_coord

# Plot MDS with ggplot2
IvsNI_MDS_W14 <- ggplot(IvsNI_W14_coord) +
  geom_point(aes(x = x, y = y,
                 colour = group,
                 shape  = group),
             size = 3) +
  expand_limits(x=c(-0.6, 0.6), y=c(-0.6, 0.6)) +
  scale_colour_manual("Treatment",
                      values = c("#e06377", "#b2b2b2")) +
  scale_shape_manual("Treatment",
                     values = c(17, 19)) +
  geom_text_repel(aes(x = x, y = y, label = animal)) +
  theme_bw(base_size = 12, base_family = "Arial") +
  ggtitle("14 WPI") +
  xlab("BCV distance 1") +
  ylab("BCV distance 2")

# Check W10 MDS plot
IvsNI_MDS_W14


### Export high quality image for all MDS plots 
IvsNI_files_MDS <- paste0(c("IvsNI_MDS_PRE", "IvsNI_MDS_W1", "IvsNI_MDS_W14"), ".pdf")
IvsNI_plots_MDS <- list(IvsNI_MDS_PRE, IvsNI_MDS_W1, IvsNI_MDS_W14)

purrr::pwalk(list(IvsNI_files_MDS, IvsNI_plots_MDS),
      ggsave,
      device    = cairo_pdf,
      path      = imgDir,
      limitsize = FALSE,
      dpi       = 300,
      height    = 9,
      width     = 10,
      units     = "in")

#####################################################
# 20 Plot: Combine all MDS plots into single figure #
#####################################################

# Set grid
IvsNI_MDS_grid <- plot_grid((IvsNI_MDS_PRE + theme(legend.position = "none")),
                      (IvsNI_MDS_W1  + theme(legend.position = "none")),
                      (IvsNI_MDS_W14 + theme(legend.position = "none")),
                      labels = c("A", "B", "C"),
                      ncol = 1,
                      scale = .96)

# Check plot
IvsNI_MDS_grid

# Get legend from one MDS plot
IvsNI_legend <- get_legend(IvsNI_MDS_PRE)

# Add legend to grid
IvsNI_MDS_grid_leg <- plot_grid(IvsNI_MDS_grid, legend, ncol = 2,
                          rel_widths = c(1, .3),
                          rel_heights = c(3, .3))

# Check plot
IvsNI_MDS_grid_leg

# Export high quality image for both plots  
IvsNI_files_MDSgrid <- paste0(c("IvsNI_MDS_grid", "IvsNI_MDS_grid_leg"), ".pdf")
IvsNI_plots_MDSgrid <- list(IvsNI_MDS_grid, IvsNI_MDS_grid_leg)

purrr::pwalk(list(IvsNI_files_MDSgrid, IvsNI_plots_MDSgrid),
             ggsave,
             device    = cairo_pdf,
             path      = imgDir,
             limitsize = FALSE,
             dpi       = 300,
             height    = 15,
             width     = 10,
             units     = "in")

##################################
# 21 Define experimental factors #
##################################

head(IvsNI_Liver_Fluke_norm$samples)

# Liver Fluke Infection factor
IvsNI_condition <- IvsNI_Liver_Fluke_norm$samples$group

# Time point factor
IvsNI_time.point <- IvsNI_Liver_Fluke_norm$samples$time.point

# Animal factor
IvsNI_animal <- IvsNI_Liver_Fluke_norm$samples$animal

# Combine Liver Fluke infection and time point into one factor
# to simplify contrats
IvsNI_cond.time <- factor(paste(IvsNI_Liver_Fluke_norm$samples$group,
                                IvsNI_Liver_Fluke_norm$samples$time.point,
                          sep="."),
                    levels = c("Control.PRE", "Control.W1", "Control.W14", "Infected.PRE", "Infected.W1", "Infected.W14"))

#################################################
# 22 Create a design matrix for paired analysis #    This is not a PAIRED analysis so we need to delete intercept column (0)
#################################################  Look at page 37 in EdgeR manual

# Create a design matrix with animal as a blocking factor     ####### I remove animal to make matrix
IvsNI_matrix_group_time.point <- model.matrix(~group + group:time.point,
                             data = IvsNI_Liver_Fluke_norm$samples)

dim(IvsNI_matrix_group_time.point)
dim(IvsNI_Liver_Fluke_norm$samples)
head(IvsNI_matrix_group_time.point)

# Rename design matrix columns for simplicity
colnames(IvsNI_matrix_group_time.point) %<>%
  
  str_replace("group", "")

head(IvsNI_matrix_group_time.point)

# Output the design matrix info
write_csv(as.data.frame(IvsNI_matrix_group_time.point),
          path = file.path(paste0(workDir, "/Liver_Fluke_design-matrix_only_cond.time-infected-vs-noinfected.csv")),
          col_names = TRUE)

#########################################
# 23 Estimate the dispersion parameters #
#########################################

# Common and trended dispersions are estimated with the
# Cox-Reid method and tagwise dispersions with the
# empirical Bayes method
IvsNI_Liver_Fluke_disp <- estimateDisp.DGEList(y       = IvsNI_Liver_Fluke_norm,
                                  design  = IvsNI_matrix_group_time.point,
                                  robust  = TRUE,
                                  verbose = TRUE)

names(IvsNI_Liver_Fluke_disp)

# Check the calculated dispersion
IvsNI_Liver_Fluke_disp$common.dispersion

# Check the calculated dispersion's square root,
# which corresponds to the biological coefficient of variation (BCV)
sqrt(IvsNI_Liver_Fluke_disp$common.dispersion)
sqrt(IvsNI_Liver_Fluke_disp$tagwise.dispersion)

# Create a matrix of the tagwise dispersion associated with gene information
IvsNI_Tagwisedisp <- cbind(IvsNI_Liver_Fluke_disp$genes, IvsNI_Liver_Fluke_disp$tagwise.dispersion)
head(IvsNI_Tagwisedisp)
dim(IvsNI_Tagwisedisp)

# Output tagwise dispersion values with gene info
IvsNI_Tagwisedisp <- as.data.frame(cbind(IvsNI_Liver_Fluke_disp$genes,
                                         IvsNI_Liver_Fluke_disp$tagwise.dispersion))
write_csv(IvsNI_Tagwisedisp,
          path = file.path(paste0(workDir, "/Liver_Fluke_Tagwise_dispersion-infected-vs-noinfected.csv")),
          col_names = TRUE)

################################
# 24 Plot: BCV and dispersions #
################################

# Create a dataframe with the dispersion values
names(IvsNI_Liver_Fluke_disp)

IvsNI_Disp <- as.data.frame(cbind(IvsNI_Liver_Fluke_disp$genes,
                                  IvsNI_Liver_Fluke_disp$tagwise.dispersion,
                                  IvsNI_Liver_Fluke_disp$common.dispersion,
                                  IvsNI_Liver_Fluke_disp$trended.dispersion,
                                  IvsNI_Liver_Fluke_disp$AveLogCPM))

colnames(IvsNI_Disp) %<>%
  str_replace("Liver_Fluke_disp\\$", "")

IvsNI_Disp %<>%
  dplyr::mutate(type_point = "Tagwise dispersion") %>%
  dplyr::mutate(type_hline = "Common dispersion") %>%
  dplyr::mutate(type_smooth = "Trended dispersion")

# Plot all dispersions
IvsNI_Liver_Fluke_BCV <- ggplot(IvsNI_Disp) +
              geom_point(aes(x = IvsNI_AveLogCPM,
                             y = sqrt(IvsNI_tagwise.dispersion),
                             fill = type_point),
                         alpha = 0.5) +
              geom_hline(aes(yintercept = sqrt(IvsNI_common.dispersion),
                             colour = type_hline)) +
              geom_smooth(aes(x = IvsNI_AveLogCPM,
                              y = sqrt(IvsNI_trended.dispersion),
                              colour = type_smooth),
                              linetype = 2) +
              scale_fill_manual("", values = c("black")) +
              scale_colour_manual("", values = c("red", "blue")) +
              theme_bw(base_size = 14, base_family = "Arial") +
              ggtitle("Estimated dispersions (NB model)") +
              xlab(expression(paste("Average ", log[2],"CPM"))) +
              ylab("Biological Coefficient of Variation")

IvsNI_Liver_Fluke_BCV

# Output high resolution plot
ggsave("Liver_Fluke_BCV-infected-vs-noinfected.pdf",
       plot = IvsNI_Liver_Fluke_BCV,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")

####################
# 25 Fit GLM model #
####################

# Fit a quasi-likelihood negative binomial generalized log-linear model
# for each tag using the design matrix and calculated dispersion
IvsNI_Liver_Fluke_fitQL <- glmQLFit(y = IvsNI_Liver_Fluke_disp,
                       design = IvsNI_matrix_group_time.point,
                       robust = TRUE)

names(IvsNI_Liver_Fluke_fitQL)
colnames(IvsNI_Liver_Fluke_fitQL$design)

###########################
# 26 Plot: QL dispersions #
###########################

# Create a dataframe with the dispersion values
names(IvsNI_Liver_Fluke_fitQL)

IvsNI_DispQL <- as.data.frame(cbind(AveLogCPM = IvsNI_Liver_Fluke_fitQL$AveLogCPM,
                              deviance = IvsNI_Liver_Fluke_fitQL$deviance,
                              df.residual.zeros = IvsNI_Liver_Fluke_fitQL$df.residual.zeros,
                              var.prior = IvsNI_Liver_Fluke_fitQL$var.prior,
                              var.post = IvsNI_Liver_Fluke_fitQL$var.post))

IvsNI_DispQL %<>%
  dplyr::mutate(type_point = "Raw dispersion (NB)") %>%
  dplyr::mutate(type_point2 = "Squeezed EB dispersion") %>%
  dplyr::mutate(type_smooth = "Trended EB dispersion")

head(IvsNI_DispQL)

# Plot all dispersions
IvsNI_Liver_Fluke_BCVQL <- ggplot(IvsNI_DispQL) +
                geom_point(aes(x = AveLogCPM,
                               y = sqrt(sqrt(deviance/df.residual.zeros)),
                               fill = type_point),
                           alpha = 0.5) +
                geom_point(aes(x = AveLogCPM,
                               y = sqrt(sqrt(var.post)),
                               colour = type_point2),
                           alpha = 0.5) +
                geom_smooth(aes(x = AveLogCPM,
                                y = sqrt(sqrt(var.prior)),
                                colour = type_smooth),
                            linetype = 2) +
                scale_fill_manual("", values = c("black")) +
                scale_colour_manual("", values = c("indianred4", "blue")) +
                theme_bw(base_size = 14, base_family = "Arial") +
                ggtitle("Estimated QL dispersions") +
                xlab(expression(paste("Average ", log[2],"CPM"))) +
                ylab("Quarter-Root Mean Deviance")

IvsNI_Liver_Fluke_BCVQL

# Output high resolution plot
ggsave("Liver_Fluke_BCVQL-infected-vs-noinfected.pdf",
       plot = IvsNI_Liver_Fluke_BCVQL,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 14,
       units     = "in")

################################
# 27 Combine dispersions plots #
################################

# Set grids
IvsNI_Liver_Fluke_grid <- plot_grid(IvsNI_Liver_Fluke_BCV,
                                    IvsNI_Liver_Fluke_BCVQL,
                      labels = c("A)", "B)"),
                      nrow = 2)

# Check plot
IvsNI_Liver_Fluke_grid

# Export high quality image
ggsave("Liver_Fluke_grid-infected-vs-noinfected.pdf",
       plot      = IvsNI_Liver_Fluke_grid,
       device    = cairo_pdf,
       path      = imgDir,
       limitsize = FALSE,
       dpi       = 300,
       height    = 10,
       width     = 8,
       units     = "in")

#######################
# 28 Save .RData file #
#######################

save.image(file = "Liver_Fluke_BovinePBMC.RData")

##########################
# 29 Save R session info #
##########################

devtools::session_info()
