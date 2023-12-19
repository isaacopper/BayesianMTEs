


################### Expanded Regression Discontinuity ##########################

# Clear console.
#cat("\014")

# Remove Plots
dev.off(dev.list()["RStudioGD"]) # Apply dev.off() & dev.list()
dev.off()

# Remove all files from workspace - do this every time so we don't use a file archived to the workspace.
rm(list = ls())

# Change Directory
setwd("/Users/iopper/Documents/ResearchProjects/BayesianMTE/")

# Define eps for jitter
eps <- sqrt(.Machine$double.eps)


################################## Import the packages #########################
library('ggplot2')
library('tibble')
library('tidyr')
library('dplyr')

library('collapse')

library('mgcv')
#install.packages('gratia')
library('gratia')
library('Matrix')
library('plgp')

library('haven')
library(matrixStats)
library('patchwork')

#source("bayesian_mte.R")
source("bayesianMTE/BayesianMTE_project/R/main_functions.R")


##############################################################
# Helper Functions
##############################################################

generate_functions <- function(hyperparams) {
  # Create eta Grid
  eta <- tibble(eta = seq(0, 1, by= .01), row_n = seq(1, 1/.01 + 1, 1))
  
  # From Vector of hyperparams to specific values
  if (length(hyperparams) == 4) {
    std_mu = exp(hyperparams[1])
    lengthscale_mu = exp(hyperparams[2])
    std_tau = exp(hyperparams[3])
    lengthscale_tau = exp(hyperparams[4])
  }
  if (length(hyperparams) == 2) {
    std_mu = exp(hyperparams[1])
    lengthscale_mu = exp(hyperparams[2])
    std_tau = exp(hyperparams[1])
    lengthscale_tau = exp(hyperparams[2])
  }
  
  # Create base kernels based on the hyperparams
  Sigma_mu <- (std_mu)^2*exp(-distance(eta$eta)/(2*lengthscale_mu^2))
  Sigma_tau <- (std_tau)^2*exp(-distance(eta$eta)/(2*lengthscale_tau^2))
  
  # Generate random functions
  mu <- t(rmvnorm(1, mean = matrix(0, nrow(Sigma_mu)), sigma = Sigma_mu + diag(1*eps, nrow(Sigma_mu))))
  mu = data.frame(mu)
  tau <- t(rmvnorm(1, mean = matrix(0, nrow(Sigma_tau)), sigma = Sigma_tau + diag(1*eps, nrow(Sigma_tau))))
  tau = data.frame(tau)
  
  #stack vectors together (one col eta, one col mu, one col tau)
  func = cbind(eta, cbind(mu, tau))
  
  # Return
  return(func)
  
}


create_data <- function(hypers,  treatment_cuts, iv_cuts, sample_size, error_std = .5) {
  
  # Drawn Random Functions
  true_functions <- generate_functions(hypers)
  
  # Randomly draw treatment propensities for each individual
  eta_i <- data.frame(eta = round(runif(sample_size)/.01)*.01)
  
  # Randomly Assign to Instrument
  r <- runif(sample_size)
  m <- matrix(0, nrow = sample_size, ncol = length(iv_cuts) + 1)
  m[,1] <- as.numeric(r < iv_cuts[1])
  if (length(iv_cuts) > 1) { 
    for (i in 2:length(iv_cuts)) {
      m[, i] <- as.numeric(r > iv_cuts[i-1] & r < iv_cuts[i])
    }
    m[, length(iv_cuts) + 1] <- as.numeric(r > iv_cuts[length(iv_cuts)])
  }
  
  # Determine (endogenous) treatment depending on instrument
  treat_i <- (eta_i < m %*% treatment_cuts)
  
  # Create Single categorical variable for IV
  for (i in 1:length(iv_cuts) + 1) {
    m[,i] <- m[,i]*i
  }
  treat_assign_i <- rowSums(m)
  
  # Combine variables
  d <- cbind(eta_i, treat_i, treat_assign_i)
  colnames(d) <- c("eta", "treat_i", "treat_assign_i")
  
  # Merge eta_i to true_functions dataset to determine mu_i & tau_i
  dataset = merge(x = true_functions, y = d, by = "eta", all = FALSE)
  #mu_i = true_functions_merge$mu
  #tau_i = true_functions_merge$tau
  
  # Outcome
  dataset <- dataset %>% mutate(Y = mu + tau*treat_i + error_std*rnorm(sample_size))
  
  # Return Data
  base_data <- data.frame(
    outcome = dataset$Y,
    treatment = dataset$treat_i,
    IV = dataset$treat_assign_i
  )
  #base_data <- base_data %>% rename(outcome = eta) %>% rename(treatment= eta.1)
  
  return(base_data)
  
}

transform_data <- function(base_data) {
  
  # Grid Length
  grid_length <- .01
  
  # Calculate observed moments & propensity scores
  p_score <- base_data %>% dplyr::group_by(IV) %>% fmean %>% dplyr::select(IV, treatment) %>% dplyr::rename(p_score = treatment)
  means <- base_data %>% dplyr::group_by(treatment, IV) %>% fmean
  var <- base_data %>% dplyr::group_by(treatment, IV) %>% fvar
  nobs <- base_data %>% dplyr::group_by(treatment, IV) %>% fnobs
  means <- means %>% dplyr::rename(mean = outcome)
  var <- var %>% dplyr::rename(var = outcome)
  nobs <- nobs %>% dplyr::rename(nobs = outcome)
  c_data <- merge(merge(means, nobs), var)
  c_data <- c_data %>% dplyr::mutate(var = var/nobs)
  c_data <- merge(c_data, p_score)
  
  # Determine average variance of the residuals
  residuals_std <- mean(c_data$var*c_data$nobs)^.5
  
  # Create Grid
  eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))
  
  # Round p_score to merge to eta_grid
  c_data <- c_data %>% dplyr::mutate(eta = round(p_score/grid_length)*grid_length)
  
  # Get row_number of prediction data --> Note, since treat are stacked below controls, treat row is order in the sequence + nrow(prediction_grid)
  merged_control <- merge(c_data %>% dplyr::filter(treatment == 0), eta) %>% dplyr::select(mean, var, row_n)
  merged_treat <- merge(c_data %>% dplyr::filter(treatment == 1), eta) %>% dplyr::select(mean, var, row_n) %>% dplyr::mutate(row_n = row_n + nrow(eta))
  merged <- rbind(merged_control, merged_treat)
  
  # Return 
  return(list("merged" = merged, "c_data" = c_data))
  
}


estimate_likelihoods <- function(lscale_seq, sigma_seq, analysis_data) {
  # Set up results
  ml_results <- matrix(0, length(lscale_seq)*length(sigma_seq), 5)
  
  # Loop 
  i <- 1
  for (l in lscale_seq) {
    for (sigma in sigma_seq) {
      
      # Randomly draw hyperparameters
      hypers <- c(sigma, l, sigma, l)
      
      # For each set of hypers estimate marginal likelihood
      ll <- log_m_likelihood(hypers, analysis_data, eta, ml_only = FALSE)
      
      # Results
      ml_results[i,1] <- hypers[1]
      ml_results[i,2] <- hypers[2]
      ml_results[i,3] <- ll$data_fit
      ml_results[i,4] <- ll$complexity_penalty
      ml_results[i,5] <- ll$log_m_likelihood
      i <- i + 1
      
    }
  }
  
  # Turn into dataframe
  ml_results <- data.frame(ml_results) %>% rename(sigma = X1, lengthscale = X2, datafit = X3, complexity = X4, log_ml = X5)
  
  # From log likelihood to likelihood
  ml_results <- ml_results %>% mutate(ml = exp(log_ml))
  ml_results <- ml_results %>% mutate(ml = ml/max(ml))
  
  
  # Return
  return(ml_results)
}



estimate_both <- function(lscale_seq, sigma_seq, analysis_data, collapsed_data) {
  # Set up results
  ml_results <- matrix(0, length(lscale_seq)*length(sigma_seq), 13)
  
  # Loop 
  i <- 1
  for (l in lscale_seq) {
    for (sigma in sigma_seq) {
      
      # Randomly draw hyperparameters
      hypers <- c(sigma, l, sigma, l)
      
      # For each set of hypers estimate marginal likelihood
      ll <- log_m_likelihood(hypers, analysis_data, eta, ml_only = FALSE)
      treat_effects <- bayesian_mte_sub(hypers, analysis_data, eta, collapsed_data, grid_length)
      
      # Results
      ml_results[i,1] <- hypers[1]
      ml_results[i,2] <- hypers[2]
      ml_results[i,3] <- ll$data_fit
      ml_results[i,4] <- ll$complexity_penalty
      ml_results[i,5] <- ll$log_m_likelihood
      ml_results[i,6] <- treat_effects$ATE$mean
      ml_results[i,7] <- treat_effects$ATE$variance
      ml_results[i,8] <- treat_effects$LATE$mean
      ml_results[i,9] <- treat_effects$LATE$variance
      ml_results[i,10] <- treat_effects$ATATE$mean
      ml_results[i,11] <- treat_effects$ATATE$variance
      ml_results[i,12] <- treat_effects$NTATE$mean
      ml_results[i,13] <- treat_effects$NTATE$variance
      
      i <- i + 1
      
    }
  }
  
  # Turn into dataframe
  colnames(ml_results) <- c("sigma", "lengthscale", "datafit", "complexity", "log_ml", "ATE_mean", "ATE_variance", "LATE_mean", "LATE_variance", "ATATE_mean", "ATATE_variance", "NTATE_mean", "NTATE_variance")
  ml_results <- data.frame(ml_results)
  
  #ml_results <- data.frame(ml_results) %>% rename(std = X1, lengthscale = X2, datafit = X3, complexity = X4, log_ml = X5)
  
  # From log likelihood to likelihood
  ml_results <- ml_results %>% mutate(ml = exp(log_ml))
  ml_results <- ml_results %>% mutate(ml = ml/max(ml))
  
  
  # Return
  return(ml_results)
}

##############################################################
# Continuous Instrument With Limited Support
##############################################################
# Create Data
observed_data <- create_data(c(-.2, 0), seq(.2, .390, .01), seq(.05, .95, .05), 100000, error_std = .1)

# Estimate
estimates_eb <- bayesian_mte(observed_data$outcome, observed_data$treatment, observed_data$IV, full_bayes = FALSE)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/PosteriorMomentFunctions_ContinuousSim.pdf", estimates_eb$predictions_plot, width = 5, height = 6, dpi = 300)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/MTEPosterior_ContinuousSim.pdf", estimates_eb$MTE_plot, width = 5, height = 6, dpi = 300)

estimates_fb <- bayesian_mte(observed_data$outcome, observed_data$treatment, observed_data$IV, full_bayes = TRUE)

##############################################################
# Simulate over multiple draws of data
##############################################################
# Binary
for (k in seq(1, 100)) { 
  
  # Create Data
  observed_data <- create_data(c(-.2, 0), seq(.25, .75, .5), c(.5), 100000, error_std = .01)
  analysis_data <- transform_data(observed_data)$merged
  
  # Create Grid
  grid_length <- .01
  eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))
  
  # Compute Likelihoods for Range of Hyperparameters
  ml_results <- estimate_likelihoods(seq(-5,5,.1), seq(-3,1,.1),analysis_data)
  
  # Add interation number
  ml_results <- ml_results %>% mutate(iter = k) %>% mutate(support = "binary")
  
  # Save results
  if (k == 1) {
    full_ml_results <- ml_results
  }
  else {
    full_ml_results <- rbind(full_ml_results, ml_results)
  }
  
  # print
  print(k)
}


# Continuous with Limited Support
for (k in seq(1, 100)) { 
  
  # Create Data
  observed_data <- create_data(c(-.2, 0, -.2, 0), seq(.1, .55, .05), seq(.1, .9, .1), 100000, error_std = .01)
  analysis_data <- transform_data(observed_data)$merged
  
  # Create Grid
  grid_length <- .01
  eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))
  
  # Compute Likelihoods for Range of Hyperparameters
  ml_results <- estimate_likelihoods(seq(-5,5,.1), seq(-3,1,.1),analysis_data)
  
  # Add interaction number
  ml_results <- ml_results %>% mutate(iter = k) %>% mutate(support = "continuous_partial_support")
  
  # Save results
  full_ml_results <- rbind(full_ml_results, ml_results)
  
  # print
  print(k)
}


# Continuous with Full Support
for (k in seq(1, 100)) { 
  
  # Create Data
  observed_data <- create_data(c(-.2, 0, -.2, 0), seq(.025, .975, .05), seq(.05, .95, .05), 100000, error_std = .01)
  analysis_data <- transform_data(observed_data)$merged
  
  # Create Grid
  grid_length <- .01
  eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))
  
  # Compute Likelihoods for Range of Hyperparameters
  ml_results <- estimate_likelihoods(seq(-5,5,.1), seq(-3,1,.1),analysis_data)
  
  # Add interaction number
  ml_results <- ml_results %>% mutate(iter = k) %>% mutate(support = "continuous_full_support")
  
  # Save results
  full_ml_results <- rbind(full_ml_results, ml_results)
  
  # print
  print(k)
}


# Keep highest value of ML in for each lengthscale
maximum_value <- full_ml_results %>% 
  group_by(lengthscale, iter, support) %>%
  slice_max(log_ml)

# Integrate over distributions
avg_value <- maximum_value %>%
  group_by(lengthscale, support) %>%
  fmean

# Rename support
avg_value$support <- ifelse(avg_value$support == "binary", "Binary", avg_value$support)
avg_value$support <- ifelse(avg_value$support == "continuous_full_support", "Continous w/ Full Support", avg_value$support)
avg_value$support <- ifelse(avg_value$support == "continuous_partial_support", "Continous w/ Partial Support", avg_value$support)

# Plot
#ggplot(maximum_value, aes(x = lengthscale, group = iter, linetype = support)) + geom_line(aes(y = log_ml))
ggplot(avg_value, aes(x = lengthscale, y = ml, linetype = support, colour = support)) + 
  geom_line(linewidth = 2) + theme_bw() + ylab("Normalized Marginal Likelihood") + xlab("Log Lengthscale") +
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "black") + labs(colour = "Instrument Type:", linetype = "Instrument Type:") + 
  theme(legend.position = c(0.25, 0.9)) +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 10)) 
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/IDingLengthscale.pdf", width = 5, height = 6, dpi = 300)

ggplot(avg_value %>% filter(abs(lengthscale) < 2.5), aes(x = lengthscale, y = ml, linetype = support, colour = support)) + 
  geom_line(linewidth = 2) + theme_bw() + ylab("Normalized Marginal Likelihood") + xlab("Log Lengthscale") +
  geom_vline(xintercept = 0.0, linetype = "dashed", color = "black") + labs(colour = "Instrument Type:", linetype = "Instrument Type:") + 
  theme(legend.position = c(0.25, 0.9)) +
  theme(legend.text = element_text(size = 8), legend.title = element_text(size = 10)) 
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/IDingLengthscaleZoomed.pdf", width = 5, height = 6, dpi = 300)



##############################################################
# Potential Problem
##############################################################
# Create a data frame with Y, Z, and T variables
df <- data.frame(outcome = rnorm(100),
                 IV = rep(c(0, 1), each = 50),
                 treatment = NA)

# Split T based on the specified proportions
df$treatment[df$IV == 0] <- sample(c(0, 1), 50, replace = TRUE, prob = c(0.55, 0.45))
df$treatment[df$IV == 1] <- sample(c(0, 1), 50, replace = TRUE, prob = c(0.45, 0.55))

# Subtract the mean of Y for each group from each value of Y using group_by and mutate
df <- df %>% 
  group_by(IV, treatment) %>% 
  mutate(outcome = outcome - mean(outcome) + 0.1) %>% 
  ungroup()

# Add back a small constant treatment effect
df$outcome <- ifelse(df$treatment == 1, df$outcome + 0.2, df$outcome)

# Transform data for calculating likelihood
analysis_data <- transform_data(df)$merged

# Create Grid
grid_length <- .01
eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))

# Compute Likelihoods for Range of Hyperparameters
ml_results <- estimate_likelihoods(seq(-5,5,.1), seq(-3,1,.1),analysis_data)

# Graph
ggplot(ml_results  %>% group_by(lengthscale) %>% slice_max(log_ml) %>% ungroup(), aes(x = lengthscale, y = ml)) +
  geom_line() + theme_bw()

# Graph
ggplot(ml_results  %>% group_by(lengthscale) %>% slice_max(log_ml) %>% ungroup(), aes(x = lengthscale)) +
  geom_line(aes(y = log_ml)) + geom_line(aes(y = complexity))  + geom_line(aes(y = datafit)) + theme_bw()

# Estimate EB with Optimal Hyperparameters
ml_results %>% slice_max(log_ml) %>% select(c('std','lengthscale'))
estimates_eb <- bayesian_mte(df$outcome, df$treatment, df$IV, input_hypers = c(-1.9,5), full_bayes = FALSE)
estimates_eb_extrap <- bayesian_mte(df$outcome, df$treatment, df$IV, input_hypers = c(-1.9,5), full_bayes = FALSE, extrapolation_uncertainty = TRUE)

