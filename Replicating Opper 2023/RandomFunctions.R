

################### Lengthscale Examples ##########################

# Clear console.
#cat("\014")

# Remove Plots
#dev.off(dev.list()["RStudioGD"]) # Apply dev.off() & dev.list()
#dev.off()

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

# used to calculate the moving average
library(purrr)


######################## Define the Function ###################################
# Create Grid
grid_length = .01
eta <- tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))

# Define funciton that plots potential functions
plot_potential_base_functions <- function(hypers, n = 3) {
  # Turn into parameters
  std_mu = exp(hypers[1])
  lengthscale_mu = exp(hypers[2])
  std_tau = exp(hypers[3])
  lengthscale_tau = exp(hypers[4])
  
  # Create base kernels based on the hyperparams
  Sigma <- (std_mu)^2*exp(-distance(eta$eta)/(2*lengthscale_mu^2))

  # Draw 
  Y <- rmvnorm(n, mean = matrix(0, nrow = nrow(Sigma)), sigma = Sigma + diag(1*eps, nrow(Sigma)))
  
  # Save results
  random_functions <- data.frame(
    x = rep(1:nrow(Sigma), times = 1),
    f1 = Y[1,],
    f2 = Y[2,],
    f3 = Y[3,]
  )
  
  # Calculate Moving Averages
  random_functions <- random_functions %>% mutate(f1_tilde = map_dbl(row_number(), ~mean(f1[.x:101]))) %>%
    mutate(f2_tilde = map_dbl(row_number(), ~mean(f2[.x:101]))) %>%
    mutate(f3_tilde = map_dbl(row_number(), ~mean(f3[.x:101])))
  
  
  #random_functions <- random_functions %>% mutate(treatment = as.numeric(rownames(random_functions)) > max(random_functions$x))
  p <- ggplot(data = random_functions, aes(x = x/nrow(Sigma))) + theme_bw() + 
    geom_line(aes(y = f1, colour = 'red')) +
    geom_line(aes(y = f2, colour = 'blue', linetype = 'dashed')) +
    geom_line(aes(y = f3, colour = 'green', linetype = 'dotdash')) +
    scale_x_continuous(breaks=c(0, 0.5, 1.0)) + scale_y_continuous(n.breaks = 3) +
    labs(title = bquote('l=' ~ e^.(hypers[2]))) + theme(plot.title = element_text(size=11)) +
    theme(legend.position="none") + xlab(expression(eta)) + ylab(expression(mu(eta)))
  
  tilde_p <- ggplot(data = random_functions, aes(x = x/nrow(Sigma))) + theme_bw() + 
    geom_line(aes(y = f1_tilde, colour = 'red')) +
    geom_line(aes(y = f2_tilde, colour = 'blue', linetype = 'dashed')) +
    geom_line(aes(y = f3_tilde, colour = 'green', linetype = 'dotdash')) +
    xlab(expression(eta)) + ylab(bquote(y[0](eta))) + theme(legend.position="none") +
    scale_x_continuous(breaks=c(0, 0.5, 1.0)) + scale_y_continuous(n.breaks = 3) +
    labs(title = bquote('l=' ~ e^.(hypers[2]))) + theme(plot.title = element_text(size=11))

  return(list("p" = p, "tilde_p" = tilde_p))
}

######################## Define the Function ###################################
# Estimate three random functions for a range of lengthscales
p1 <- plot_potential_base_functions(c(-3,-3,1,1))
p2 <- plot_potential_base_functions(c(-3,-2,1,1))
p3 <- plot_potential_base_functions(c(-3,-1,1,1))
p4 <- plot_potential_base_functions(c(-3,-0.5,1,1))
p5 <- plot_potential_base_functions(c(-3,0,1,1))
p6 <- plot_potential_base_functions(c(-3,0.5,1,1))
p7 <- plot_potential_base_functions(c(-3,1,1,1))
p8 <- plot_potential_base_functions(c(-3,2,1,1))
p9 <- plot_potential_base_functions(c(-3,3,1,1))

# Combine the 9 subplots into one plot, separately for mu and y_0
p <- (p1$p + p2$p + p3$p) / (p4$p + p5$p + p6$p) /(p7$p + p8$p + p9$p)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/VariousLengthscales.pdf", p, width = 6, height = 8, dpi = 300)

tilde_p <- (p1$tilde_p + p2$tilde_p + p3$tilde_p) / (p4$tilde_p + p5$tilde_p + p6$tilde_p) /(p7$tilde_p + p8$tilde_p + p9$tilde_p)
ggsave("/Users/iopper/Dropbox/Research Papers/FromLATEtoATE/Figures/VariousLengthscalesY0.pdf", p, width = 6, height = 8, dpi = 300)

############### Look at a very small lengthscale ###############################
# Extremly small lengthscale
ptiny <- plot_potential_base_functions(c(-3,-10,1,1))


