

################### Expanded Regression Discontinuity ##########################

# Clear console.
cat("\014")

# Remove Plots
#dev.off(dev.list()["RStudioGD"]) # Apply dev.off() & dev.list()
#dev.off()

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
#source("bayesianMTE/BayesianMTE_project/R/main_functions.R")

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

remove.packages("/Users/iopper/Desktop/BayesianMTEs-main/bayesianMTE")
devtools::install("/Users/iopper/Desktop/BayesianMTEs-main/bayesianMTE")
library("bayesianMTE")


library(tidyr)

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

# 10x the Sample Size
#data <- rbind(data, data, data, data, data, data, data, data, data, data)
#set.seed(1215)
#data <- sample_n(data, 1709, replace = TRUE)


##############################################################
# Bayesian MTEs 
##############################################################
# Full Bayes
library(patchwork)
estimates <- bayesianMTE::bayesian_mte(data$any_visit_180p_180, data$ohp_all_ever_matchn_30sep2009, data$treatment, full_bayes = TRUE)

##############################################################
# Bayesian MTEs vs. Linear Extrapolation vs. Heckit
##############################################################
# Helper function to output estimates and standard errors given weights & vcv & params
calculate_treatment <- function(params, vcv, weights) {
  
  # Returns point estimates and standard errors
  return(c("Estimate" = weights %*% params, "StdError" = (weights %*% vcv %*% weights)^.5))
}

# Determine Cut-Points
first_stage <- lm(formula = ohp_all_ever_matchn_30sep2009 ~ treatment, data = data)
lower_cut <- first_stage$coefficients[1]
upper_cut <- first_stage$coefficients[1] + first_stage$coefficients[2]

# Calculate Threshold Values for Linear Model 
#--> Assumes J(u) = u - .5
treat_lower <- (lower_cut - 1)/2
treat_upper <- (upper_cut - 1)/2
control_lower <- lower_cut/2
control_upper <- upper_cut/2

# Add to data
data <- data %>% mutate(linear_J = treat_lower*ohp_all_ever_matchn_30sep2009*(1 - treatment) + treat_upper*ohp_all_ever_matchn_30sep2009*treatment + control_lower*(1 - ohp_all_ever_matchn_30sep2009)*(1 - treatment) + control_upper*(1 - ohp_all_ever_matchn_30sep2009)*treatment)

# Parameter Estimation
linear_J_reg <- lm(any_visit_180p_180 ~ linear_J*(1 - ohp_all_ever_matchn_30sep2009) + ohp_all_ever_matchn_30sep2009 + linear_J*ohp_all_ever_matchn_30sep2009, data = data)
linear_J_params <- linear_J_reg$coefficients
linear_J_vcov <- vcov(linear_J_reg)

# Linear Estimates
ATE_linear <- calculate_treatment(linear_J_params, linear_J_vcov, c(0,0,1,0))
ATATE_linear <- calculate_treatment(linear_J_params, linear_J_vcov, c(0,0,1,(lower_cut -1)/2))
LATE_linear <- calculate_treatment(linear_J_params, linear_J_vcov, c(0,0,1,(upper_cut^2 - lower_cut^2 - upper_cut + lower_cut)/(2 *(upper_cut - lower_cut))))
NTATE_linear <- calculate_treatment(linear_J_params, linear_J_vcov, c(0,0,1,upper_cut/2))

# MTEs
all_MTEs <- estimates$MTE %>% mutate(tau_hat_linear_extrap = linear_J_params[3] + (eta + .005 - .5)*linear_J_params[4])
all_MTEs <- all_MTEs %>% mutate(std_tau_hat_linear_extrap = (linear_J_vcov[3,3] + (eta + .005 - .5)^2*linear_J_vcov[4,4] + 2*(eta + .005 - .5)*linear_J_vcov[3,4])^.5)

# Calculate Threshold Values for Heckit
treat_lower <- integrate(qnorm, 0, lower_cut)$value/lower_cut
treat_upper <- integrate(qnorm, 0, upper_cut)$value/upper_cut
control_lower <- integrate(qnorm, lower_cut, 1)$value/(1 - lower_cut)
control_upper <- integrate(qnorm, upper_cut, 1)$value/(1 - upper_cut)

# Add to data
data <- data %>% mutate(heckit_J = treat_lower*ohp_all_ever_matchn_30sep2009*(1 - treatment) + treat_upper*ohp_all_ever_matchn_30sep2009*treatment + control_lower*(1 - ohp_all_ever_matchn_30sep2009)*(1 - treatment) + control_upper*(1 - ohp_all_ever_matchn_30sep2009)*treatment)

# Parameter Estimation
heckit_reg <- lm(any_visit_180p_180 ~ heckit_J*(1 - ohp_all_ever_matchn_30sep2009) + ohp_all_ever_matchn_30sep2009 + heckit_J*ohp_all_ever_matchn_30sep2009, data = data)
heckit_params <- heckit_reg$coefficients
heckit_vcov <- vcov(heckit_reg)

# Heckit Estimates
ATE_heckit <- calculate_treatment(heckit_params, heckit_vcov, c(0,0,1,0))
ATATE_heckit <- calculate_treatment(heckit_params, heckit_vcov, c(0,0,1,integrate(qnorm, 0, lower_cut)$value/lower_cut))
LATE_heckit <- calculate_treatment(heckit_params, heckit_vcov, c(0,0,1,integrate(qnorm,lower_cut, upper_cut)$value/(upper_cut - lower_cut)))
NTATE_heckit <- calculate_treatment(heckit_params, heckit_vcov, c(0,0,1,integrate(qnorm, upper_cut, 1)$value/(1 - upper_cut)))

# MTEs
all_MTEs <- all_MTEs %>% mutate(inv_cdf = qnorm(eta + .005))
all_MTEs <- all_MTEs %>% mutate(tau_hat_heckit = heckit_params[3] + (inv_cdf)*heckit_params[4])
all_MTEs <- all_MTEs %>% mutate(std_tau_hat_heckit = (heckit_vcov[3,3] + (inv_cdf)^2*heckit_vcov[4,4] + 2*(inv_cdf)*heckit_vcov[3,4])^.5)

# Graph MTEs
all_MTEs <- all_MTEs %>% select(!c(q5, q95, inv_cdf)) %>% rename(std_tau_hat_bayesianMTE = std_tau_hat) %>% rename(tau_hat_bayesianMTE = tau_hat)
MTEs <- all_MTEs %>% pivot_longer(cols = -eta,   names_to = c(".value", "estimate"),
                                  names_pattern = "(.*)_(bayesianMTE|linear_extrap|heckit)")
MTEs <- MTEs %>% mutate(q5 = tau_hat - 1.96*std_tau_hat, q95 = tau_hat + 1.96*std_tau_hat)

MTE_estimates <- ggplot(MTEs, aes(x = eta, y = tau_hat, color = estimate)) +
  geom_line(size = 1.5) +
  geom_line(aes(y = q5), linetype = "dashed", size = 1.5) +
  geom_line(aes(y = q95), linetype = "dashed", size = 1.5) +
  labs(color = "Estimation Approach") + 
  ggplot2::xlab(expression(eta)) + ggplot2::ylab("Conditional average treatment effect") +
  scale_color_manual(values = c("bayesianMTE" = "blue", "heckit" = "red", "linear_extrap" = "#006400"),
                     labels = c("Bayesian MTEs", "Heckit", "Linear Extrapolation")) +
  theme_minimal()
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/MTEEstimates.pdf", MTE_estimates, width = 5, height = 6, dpi = 300)



##############################################################
# Latex Table
##############################################################
tab <- paste(
sprintf("Avg. Treat Effect (ATE) & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\", estimates$ATE$mean, estimates$ATE$variance^.5, ATE_linear[1], ATE_linear[2], ATE_heckit[1], ATE_heckit[2]),
sprintf("Always Taker ATE (ATATE) & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\", estimates$ATATE$mean, estimates$ATATE$variance^.5, ATATE_linear[1], ATATE_linear[2], ATATE_heckit[1], ATATE_heckit[2]),
sprintf("Local ATE (LATE) & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\", estimates$LATE$mean, estimates$LATE$variance^.5, LATE_linear[1], LATE_linear[2], LATE_heckit[1], LATE_heckit[2]),
sprintf("Never Taker ATE (NTATE) & %.3f & %.3f & %.3f & %.3f & %.3f & %.3f \\", estimates$NTATE$mean, estimates$NTATE$variance^.5, NTATE_linear[1], NTATE_linear[2], NTATE_heckit[1], NTATE_heckit[2])
)

c(tab)


