
################### Expanded Regression Discontinuity ##########################

# Clear console.
cat("\014")

# Remove Plots
dev.off(dev.list()["RStudioGD"]) # Apply dev.off() & dev.list()
dev.off()

# Remove all files from workspace - do this every time so we don't use a file archived to the workspace.
rm(list = ls())

# Change Directory
#use this for Isaac
setwd("/Users/iopper/Documents/ResearchProjects/BayesianMTE/")

#use this for aarya
#setwd("C:/Users/asuryava/Documents/Projects/Bayesian")


# Define eps for jitter
eps <- sqrt(.Machine$double.eps)

################################## Import the packages #########################
# library('ggplot2')
# library('tibble')
# library('tidyr')
# library('dplyr')
# 
# library('collapse')
# #
# library('mgcv')
#install.packages('gratia')
# library('gratia')
# library('Matrix')
# library('plgp')
#
library('haven')
# library(matrixStats)

#library("bayesianMTE")
source("bayesianMTE/BayesianMTE_project/R/main_functions.R")

library(matrixStats)
library(magrittr)
library(stats)
library(mvtnorm)
library(dplyr)
library(tibble)
library(collapse)
library(plgp)
library(ggplot2)
library(methods)


##############################################################
# Import Data
##############################################################
# Read in data
descriptive_vars <- read_dta("SampleData/OHIE_Public_Use_Files/OHIE_Data/oregonhie_descriptive_vars.dta")
stprograms_data <- read_dta("SampleData/OHIE_Public_Use_Files/OHIE_Data/oregonhie_stateprograms_vars.dta")
patterns <- read_dta("SampleData/OHIE_Public_Use_Files/OHIE_Data/oregonhie_patterns_vars.dta")

# Merge
data <- merge(descriptive_vars, stprograms_data)
data <- merge(data, patterns)

##############################################################
# Run Code
##############################################################
# Empirical Bayes Approach
estimates_eb <- bayesian_mte(data$any_visit_180p_180, data$ohp_all_ever_matchn_30sep2009, data$treatment, full_bayes = FALSE)

# Output Graphs
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/UnconditionalRandomFunctions.pdf", estimates_eb$prior_plot, width = 5, height = 6, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/ConditionalRandomFunctions.pdf", estimates_eb$posterior_plot, width = 5, height = 6, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/PosteriorMomentFunctions.pdf", estimates_eb$predictions_plot, width = 5, height = 6, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/MTEPosterior.pdf", estimates_eb$MTE_plot, width = 5, height = 6, dpi = 300)

ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/ATE_prior_posterior.pdf", estimates_eb$ATE$graph, width = 5, height = 5, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/LATE_prior_posterior.pdf", estimates_eb$LATE$graph, width = 5, height = 5, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/ATATE_prior_posterior.pdf", estimates_eb$ATATE$graph, width = 5, height = 5, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/NTATE_prior_posterior.pdf", estimates_eb$NTATE$graph, width = 5, height = 5, dpi = 300)


# Full Bayes
estimates <- bayesian_mte(data$any_visit_180p_180, data$ohp_all_ever_matchn_30sep2009, data$treatment, full_bayes = TRUE)

##############################################################
# Statistical vs Extrapolation
##############################################################
# Full Sample
estimates_freq <- bayesian_mte(data$any_visit_180p_180, data$ohp_all_ever_matchn_30sep2009, data$treatment, full_bayes = FALSE, frequentist_uncertainty = TRUE, input_hypers = estimates_eb$hypers)
estimates_extrap <- bayesian_mte(data$any_visit_180p_180, data$ohp_all_ever_matchn_30sep2009, data$treatment, full_bayes = FALSE, extrapolation_uncertainty = TRUE, input_hypers = estimates_eb$hypers)

# Output Graphs
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/PosteriorMomentFunctions_ExtrapOnly.pdf", estimates_extrap$predictions_plot, width = 5, height = 6, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/MTEPosterior_ExtrapOnly.pdf", estimates_extrap$MTE_plot, width = 5, height = 6, dpi = 300)


# N = 1,000
sampled_data <- sample_n(data, 1000, replace = TRUE)
estimates_freq <- bayesian_mte(sampled_data$any_visit_180p_180, sampled_data$ohp_all_ever_matchn_30sep2009, sampled_data$treatment, full_bayes = FALSE, frequentist_uncertainty = TRUE, input_hypers = estimates_eb$hypers)
estimates_extrap <- bayesian_mte(sampled_data$any_visit_180p_180, sampled_data$ohp_all_ever_matchn_30sep2009, sampled_data$treatment, full_bayes = FALSE, extrapolation_uncertainty = TRUE, input_hypers = estimates_eb$hypers)


