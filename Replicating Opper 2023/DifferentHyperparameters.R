
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
  ml_results <- data.frame(ml_results) %>% rename(std = X1, lengthscale = X2, datafit = X3, complexity = X4, log_ml = X5)
  
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

# Add Prior SD and Calculate Posterior -- both via Full Bayes and Empirical Bayes
var_estimates <- function(initial_data, prior_sd) {
  
  # Calculate priors
  initial_data <- initial_data %>% mutate(normed_lengthscale_prior = stats::dnorm((lengthscale)/prior_sd)/stats::dnorm(0))
  
  # Calculate posterior
  initial_data <- initial_data %>% mutate(posterior = normed_sigma_prior*normed_lengthscale_prior*ml)
  
  # Calculate weighted means and standard deviations
  ATE_variance <- sum(initial_data$ATE_variance * initial_data$posterior) / sum(initial_data$posterior) + sum(initial_data$posterior * (initial_data$ATE_mean - weighted.mean(initial_data$ATE_mean, initial_data$posterior))^2) / (sum(initial_data$posterior))
  LATE_variance <- sum(initial_data$LATE_variance * initial_data$posterior) / sum(initial_data$posterior) + sum(initial_data$posterior * (initial_data$LATE_mean - weighted.mean(initial_data$LATE_mean, initial_data$posterior))^2) / (sum(initial_data$posterior))
  ATATE_variance <- sum(initial_data$ATATE_variance * initial_data$posterior) / sum(initial_data$posterior) + sum(initial_data$posterior * (initial_data$ATATE_mean - weighted.mean(initial_data$ATATE_mean, initial_data$posterior))^2) / (sum(initial_data$posterior))
  NTATE_variance <- sum(initial_data$NTATE_variance * initial_data$posterior) / sum(initial_data$posterior) + sum(initial_data$posterior * (initial_data$NTATE_mean - weighted.mean(initial_data$NTATE_mean, initial_data$posterior))^2) / (sum(initial_data$posterior))
  
  # Empirical Bayes
  eb <- initial_data %>%
    slice_max(posterior)
  
  # Return
  return(list("ATE" = ATE_variance, "LATE" = LATE_variance, "ATATE" = ATATE_variance, "NTATE" = NTATE_variance, "ATE_EB" = eb$ATE_variance, "LATE_EB" = eb$LATE_variance, "ATATE_EB" = eb$ATATE_variance, "NTATE_EB" = eb$NTATE_variance))
}

##############################################################
# Likelihood for OHIE data
##############################################################
# Data
descriptive_vars <- read_dta("SampleData/OHIE_Public_Use_Files/OHIE_Data/oregonhie_descriptive_vars.dta")
stprograms_data <- read_dta("SampleData/OHIE_Public_Use_Files/OHIE_Data/oregonhie_stateprograms_vars.dta")
patterns <- read_dta("SampleData/OHIE_Public_Use_Files/OHIE_Data/oregonhie_patterns_vars.dta")

# Merge
data <- merge(descriptive_vars, stprograms_data)
data <- merge(data, patterns)

base_data <- data.frame(
  outcome = data$any_visit_180p_180,
  treatment = data$ohp_all_ever_matchn_30sep2009,
  IV = data$treatment
)

# Turn OHIE to format we need
analysis_data <- transform_data(base_data)$merged
collapsed_data <- transform_data(base_data)$c_data


# Create Grid
grid_length <- .01
eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))

# Compute Likelihoods for Range of Hyperparameters
ml_results <- estimate_likelihoods(seq(-5,5,.1), seq(-4,4,.1),analysis_data)

# Add prior and posterior
ml_results <- ml_results %>% mutate(prior = stats::dnorm(lengthscale/.5)*stats::dnorm((std+0.5)/1.5))
ml_results <- ml_results %>% mutate(prior = prior/max(prior))

# Add posterior
ml_results <- ml_results %>% mutate(posterior = prior*ml)

# Keep highest value of ML in for each lengthscale
maximum_value <- ml_results %>%
  group_by(lengthscale) %>%
  slice_max(log_ml)

# Add prior
maximum_value <- maximum_value %>% mutate(normed_prior = stats::dnorm(lengthscale/.5)/stats::dnorm(0/.5))

# Plot
ggplot(maximum_value %>% filter(abs(lengthscale) < 3), aes(x = lengthscale)) + geom_line(aes(y = ml)) + 
  geom_line(aes(y = normed_prior),linetype="twodash") + xlab("Log-Lengthscale") + ylab("Normalized likelihood") +
  theme_bw() + theme(legend.position = "right") +
  scale_linetype_manual(name = "Likelihood", values = c("twodash", "solid"), 
                        labels = c("Prior", "Marginal likelihood"))
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/Lengthscale_MLandPrior.pdf", width = 5, height = 5, dpi = 300)

# Keep highest value of posterior in for each lengthscale
maximum_value <- ml_results %>%
  group_by(lengthscale) %>%
  slice_max(posterior)

# Graph
ggplot(maximum_value %>% filter(abs(lengthscale) < 3), aes(x = lengthscale, y = posterior)) + 
  geom_line() + xlab("Log-Lengthscale") + ylab("Normalized likelihood") + theme_bw()
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/Lengthscale_Posterior.pdf", width = 5, height = 5, dpi = 300)

# Relationship
ggplot(maximum_value, aes(x = lengthscale)) + geom_line(aes(y = std))

# Keep highest value of ML in for each std
maximum_value <- ml_results %>%
  group_by(std) %>%
  slice_max(log_ml)

# Add prior
maximum_value <- maximum_value %>% mutate(normed_prior = stats::dnorm((std+0.5)/1.5)/stats::dnorm(0))

# Plot
ggplot(maximum_value %>% filter(abs(std) < 3), aes(x = std)) + geom_line(aes(y = ml, linetype = "Marginal likelihood"), linetype = "solid", show.legend = TRUE) + 
  geom_line(aes(y = normed_prior, linetype = "Prior"),linetype="twodash",show.legend = TRUE) + xlab("Log-Sigma") + ylab("Normalized likelihood") +
  theme_bw() + theme(legend.position = "right") + 
  scale_linetype_manual(name = "Legend", values = c("twodash", "solid"), 
                                    labels = c("Prior", "Marginal likelihood"))
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/Sigma_MLandPrior.pdf", width = 5, height = 5, dpi = 300)

# Keep highest value of ML in for each std
maximum_value <- ml_results %>%
  group_by(std) %>%
  slice_max(posterior)

# Graph posterior
ggplot(maximum_value  %>% filter(abs(std) < 3), aes(x = std, y = posterior)) + 
  geom_line() + xlab("Log-Sigma") + ylab("Normalized likelihood") + theme_bw()
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/Sigma_Posterior.pdf", width = 5, height = 5, dpi = 300)


##############################################################
# Estimate for Various Parameters
##############################################################
# Keep highest value of posterior in for each lengthscale
maximum_value <- ml_results %>%
  group_by(lengthscale) %>%
  slice_max(posterior)

# Smooth the assumed value of sigma
maximum_value <- as_tibble(maximum_value)
gam_est <- gam(std ~ s(lengthscale), data = maximum_value)
maximum_value <- maximum_value %>% mutate(sigma_pred = predict(gam_est, newdata = maximum_value))

# From tibble to data frame b/c R is weird and annoying
maximum_value <- data.frame(maximum_value)

# Estimate for every lengthscale value
param_results <- matrix(0, 101, 11)

for (k in seq(1, 101)) {
  treat_effects <- bayesian_mte_sub(c(maximum_value[k,9], maximum_value[k,2]), analysis_data, eta, collapsed_data, grid_length)

  # Results
  param_results[k,1] <- maximum_value[k,9]
  param_results[k,2] <- maximum_value[k,2]
  param_results[k,3] <- treat_effects$ATE$mean
  param_results[k,4] <- treat_effects$ATE$variance
  param_results[k,5] <- treat_effects$LATE$mean
  param_results[k,6] <- treat_effects$LATE$variance
  param_results[k,7] <- treat_effects$ATATE$mean
  param_results[k,8] <- treat_effects$ATATE$variance
  param_results[k,9] <- treat_effects$NTATE$mean
  param_results[k,10] <- treat_effects$NTATE$variance
  param_results[k,11] <- maximum_value[k,6]
}

# Name columns
colnames(param_results) <- c("sigma", "lengthscale", "ATE_mean", "ATE_variance", "LATE_mean", "LATE_variance", "ATATE_mean", "ATATE_variance", "NTATE_mean", "NTATE_variance", "ML")
param_results <- data.frame(param_results)

# Pivot Longer
est_results <- param_results %>% pivot_longer(cols = !c("lengthscale", "sigma", "ML"), names_to = c("estimand", "mean_or_var"), 
                                              names_sep = "_", values_to = "values")

# Graph
ggplot(est_results %>% filter(mean_or_var == "variance"), aes(x = lengthscale, y = values, linetype = estimand, colour = estimand)) + 
  geom_line(size = 2) + theme_bw() + ylab("Posterior Variance") + xlab("Log-Lengthscale") + 
  labs(colour = "Estimand:", linetype = "Estimand:") + theme(legend.position = c(0.8, 0.8)) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/VarianceAndLengthscaleFull.pdf", width = 5, height = 6, dpi = 300)

ggplot(est_results %>% filter(abs(lengthscale -.25) < 1.75) %>% filter(mean_or_var == "variance"), aes(x = lengthscale, y = values, linetype = estimand, colour = estimand)) + 
  geom_line(size = 2) + theme_bw() + ylab("Posterior Variance") + xlab("Log-Lengthscale") + 
  labs(colour = "Estimand:", linetype = "Estimand:") + theme(legend.position = c(0.8, 0.8)) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/VarianceAndLengthscale.pdf", width = 5, height = 6, dpi = 300)

##############################################################
# How does the prior variance
##############################################################
# Estimate both variance and marginal likelihood
full_estimates <- estimate_both(seq(-5,5,.1), seq(-4,4,.1),analysis_data, collapsed_data)

# Add prior for sigma
full_estimates <- full_estimates %>% mutate(normed_sigma_prior = stats::dnorm((sigma+0.5)/1.5)/stats::dnorm(0))

# Estimate for every SD value
results <- matrix(0, 41, 9)
i <- 1
for (v in seq(0,2,.05)) {
  if (v == 0) {
   v <- .001 
  }
  
  # Estimate
  est <- var_estimates(full_estimates, v)
  
  # store
  results[i,1] = v
  results[i,2] = est$ATE
  results[i,3] = est$LATE
  results[i,4] = est$ATATE
  results[i,5] = est$NTATE
  results[i,6] = est$ATE_EB
  results[i,7] = est$LATE_EB
  results[i,8] = est$ATATE_EB
  results[i,9] = est$NTATE_EB
  
  i <- i + 1
  
}

# Save
colnames(results) <- c("variance", "ATE", "LATE", "ATATE", "NTATE", "ATE_EB", "LATE_EB", "ATATE_EB", "NTATE_EB")
results <- data.frame(results)

# Graph
ggplot(results, aes(x = variance, y = LATE)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of LATE") + xlab("Prior SD of Lengthscale")
ggplot(results, aes(x = variance, y = ATE)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of ATE") + xlab("Prior SD of Lengthscale")
ggplot(results, aes(x = variance, y = NTATE)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of ATATE") + xlab("Prior SD of Lengthscale")
ggplot(results, aes(x = variance, y = ATATE)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of NTATE") + xlab("Prior SD of Lengthscale")

# Graph Empirical Bayes
ggplot(results, aes(x = variance, y = LATE_EB)) + geom_smooth(se = FALSE, colour = 'black') + theme_bw() +
  ylab("Posterior Variance of LATE") + xlab("Prior SD of Lengthscale")
ggplot(results, aes(x = variance, y = ATE_EB)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of ATE") + xlab("Prior SD of Lengthscale")
ggplot(results, aes(x = variance, y = NTATE_EB)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of ATATE") + xlab("Prior SD of Lengthscale")
ggplot(results, aes(x = variance, y = ATATE_EB)) + geom_line() + theme_bw() +
  ylab("Posterior Variance of NTATE") + xlab("Prior SD of Lengthscale")

# Normalize
results <- results %>% mutate(LATE = LATE/max(LATE))
results <- results %>% mutate(LATE_EB = LATE_EB/max(LATE_EB))
results <- results %>% mutate(ATE = ATE/max(ATE))
results <- results %>% mutate(ATE_EB = ATE_EB/max(ATE_EB))
results <- results %>% mutate(ATATE = ATATE/max(ATATE))
results <- results %>% mutate(ATATE_EB = ATATE_EB/max(ATATE_EB))
results <- results %>% mutate(NTATE = NTATE/max(NTATE))
results <- results %>% mutate(NTATE_EB = NTATE_EB/max(NTATE_EB))

# Reshape
results <- results %>% pivot_longer(cols = !c(variance), names_to = "estimand", values_to = "estimate")

# Graph
ggplot(results %>% filter(!stringr::str_ends(estimand, "_EB")), aes(x = variance, y = estimate, colour = estimand, linetype = estimand)) + geom_line(size = 2) +
  theme_bw() + ylab("Normalized Posterior Variance") + xlab("SD of Lengthscale Hyperprior") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") + labs(colour = "Estimand:", linetype = "Estimand:") + theme(legend.position = c(0.8, 0.2)) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/VarianceAndPriorSD.pdf", width = 5, height = 6, dpi = 300)

ggplot(results %>% filter(stringr::str_ends(estimand, "_EB")), aes(x = variance, y = estimate, colour = estimand, linetype = estimand)) + geom_smooth(se = FALSE) +
  theme_bw() + ylab("Normalized Posterior Variance") + xlab("SD of Lengthscale Hyperprior")

eb_results <- results %>% filter(stringr::str_ends(estimand, "_EB")) %>%
  mutate(estimand = gsub("_EB", "", estimand))
ggplot(eb_results, aes(x = variance, y = estimate, colour = estimand, linetype = estimand)) + geom_line(size = 2) +
  theme_bw() + ylab("Normalized Posterior Variance") + xlab("SD of Lengthscale Hyperprior") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "black") + labs(colour = "Estimand:", linetype = "Estimand:") + theme(legend.position = c(0.15, 0.15)) +
  theme(legend.text = element_text(size = 12), legend.title = element_text(size = 14))
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/EBVarianceAndPriorSD.pdf", width = 5, height = 6, dpi = 300)

# Add prior for lengthscale, posterior, and then normalize
full_estimates <- full_estimates %>% mutate(prior_l = stats::dnorm(lengthscale/.5)) %>%
  mutate(posterior = normed_sigma_prior*prior_l*ml) %>% mutate(posterior = posterior/max(posterior))

# Plot posterior
ggplot(full_estimates %>% filter(abs(lengthscale) < 1.5) %>% filter(abs(sigma + 2) < 1.5), aes(x = sigma, y = lengthscale, z = posterior)) + geom_contour_filled() 

# Accept/reject
full_estimates <- full_estimates %>% mutate(r = runif(nrow(full_estimates), min = 0, max = 1))
reshaped <- full_estimates %>% filter(r < posterior) %>% select(!c("datafit", "complexity","log_ml","ml","normed_sigma_prior","prior_l", "posterior","r")) %>%
  pivot_longer(cols = !c(sigma, lengthscale), names_to = c('estimator', 'statistic'), names_sep = "_", values_to = 'values')
reshaped <- reshaped %>% pivot_wider(names_from = statistic, values_from = values)

variation_by_hypers <- ggplot(data = reshaped %>% filter(estimator != 'ATE'), aes(x = mean, y = variance, color = estimator, shape = estimator)) + geom_point(alpha = .7) +
  scale_y_continuous(trans='log10') + theme_bw() +
  xlab("Posterior Mean") + ylab("Posterior Variance") + labs(colour = "Estimand:", shape = "Estimand:") + theme(legend.position = c(0.8, 0.8))
variation_by_hypers
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/VariationOverHypers.pdf", variation_by_hypers, width = 5, height = 6, dpi = 300)







