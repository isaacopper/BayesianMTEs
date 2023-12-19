
################### Expanded Regression Discontinuity ##########################

# Clear console.
#cat("\014")

# Remove Plots
#dev.off(dev.list()["RStudioGD"]) # Apply dev.off() & dev.list()
#dev.off()

# Remove all files from workspace - do this every time so we don't use a file archived to the workspace.
#rm(list = ls())

# Change Directory
#use this for Isaac
setwd("/Users/iopper/Documents/ResearchProjects/BayesianMTE/")

#use this for aarya
#setwd("C:/Users/asuryava/Documents/Projects/Bayesian")


# Define eps for jitter
eps <- sqrt(.Machine$double.eps)

################################## Import the packages #########################

library('haven')

#library("bayesianMTE")
source("bayesianMTE/BayesianMTE_project/R/main_functions.R")

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



# ################################## Asymptotics ##############################
ptm <- proc.time()
for (i in seq(1, 10)) {
  
  # Set Up file to save
  sim_results <- matrix(0, 9, 9)
  row_n <- 1
  
  for (sample_size in c(100, 500, 1000, 5000, 10000, 50000, 100000, 1000000, 10000000)) {
    
    # Randomly sample
    sampled_data <- sample_n(data, sample_size, replace = TRUE)
    
    # Run the code
    out <- bayesian_mte(sampled_data$any_visit_180p_180, sampled_data$ohp_all_ever_matchn_30sep2009, sampled_data$treatment, full_bayes = TRUE, hyperparameter_draws = 1000)
    
    # Add output to where we're storing it 
    sim_results[row_n, 1] <- sample_size
    sim_results[row_n, 2] <- out$ATE$mean
    sim_results[row_n, 3] <- out$ATE$variance
    sim_results[row_n, 4] <- out$LATE$mean
    sim_results[row_n, 5] <- out$LATE$variance
    sim_results[row_n, 6] <- out$ATonAT$mean
    sim_results[row_n, 7] <- out$ATonAT$variance
    sim_results[row_n, 8] <- out$ATonNT$mean
    sim_results[row_n, 9] <- out$ATonNT$variance
    row_n <- 1 + row_n
    
  }
  
  # Turn in to data frame
  sim_results <- data.frame(sim_results) %>% rename(nobs = X1, ATE_mean = X2, ATE_var = X3, LATE_mean = X4, LATE_var = X5, ATATE_mean = X6, ATATE_var = X7, NTATE_mean = X8, NTATE_var = X9)
  sim_results <- sim_results %>% mutate(iter_n = i)
  
  # Concat
  if (i == 1) {
    full_sim_results <- sim_results
  }
  if (i > 1) {
    full_sim_results <- rbind(sim_results, full_sim_results)
  }
  
  # Display
  print(i)
  print(proc.time() - ptm)
  
}

# Add MSE
full_sim_results <- full_sim_results %>% mutate(ATE_MSE = (ATE_mean)^2 + ATE_var) %>%
  mutate(LATE_MSE = (LATE_mean)^2 + LATE_var) %>% mutate(ATATE_MSE = (ATATE_mean)^2 + ATATE_var) %>%
  mutate(NTATE_MSE = (NTATE_mean)^2 + NTATE_var)

# Reshape
full_sim_results <- full_sim_results %>% pivot_longer(cols = !c(nobs, iter_n), names_to = c('estimator', 'statistic'), names_sep = "_", values_to = 'values')

# Graphs
asymptotic_variance <- ggplot(data = full_sim_results %>% filter(statistic == 'var'), aes(x = nobs, y = values, colour = estimator, linetype = estimator)) + 
  geom_smooth(se = FALSE) + 
  theme_bw() +  xlab("Sample Size") + ylab("Posterior Variance")  +
  scale_x_log10() + scale_y_log10() + theme(legend.position = c(0.2, 0.2)) + labs(colour = "Estimand:", linetype = "Estimand:")
asymptotic_variance
asymptotic_variance <- ggplot(data = full_sim_results %>% filter(statistic == 'var'), aes(x = nobs, y = values, colour = estimator, linetype = estimator)) + 
  geom_smooth(se = FALSE) + 
  theme_bw() +  xlab("Sample Size") + ylab("Posterior Variance")  +
  scale_x_log10() + theme(legend.position = c(0.8, 0.8)) + labs(colour = "Estimand:", linetype = "Estimand:")
asymptotic_variance


# 
# #MSE
# asymptotic_MSE <- ggplot(data = full_sim_results %>% filter(statistic == 'MSE'), aes(x = nobs, y = values, colour = estimator, linetype = estimator)) + 
#   geom_smooth(se = FALSE) + 
#   theme_bw() +  xlab("Sample Size") + ylab("Mean-Squared Error")  +
#   scale_x_log10() 
# asymptotic_MSE
# 
# # Bias
# asymptotic_bias <- ggplot(data = full_sim_results %>% filter(statistic == 'mean'), aes(x = nobs, y = values, colour = estimator, linetype = estimator)) + 
#   geom_smooth(se = FALSE) + 
#   theme_bw() +  xlab("Sample Size") + ylab("Bias")  +
#   scale_x_log10() 
# 
