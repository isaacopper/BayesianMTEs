#this imports dependencies using roxygen
#' @import matrixStats
#' @import magrittr
#' @import stats
#' @import mvtnorm
#' @import dplyr
#' @import tibble
#' @import collapse
#' @import plgp
#' @import ggplot2
#' @import methods
#' @import patchwork



##############################################################
##############################################################
# Support Functions
##############################################################
##############################################################
# Function to transform stationary kernel to non-stationary kernel

transform_kernel <- function(Sigma_mu, Sigma_tau) {
  # Create matrix to divide by
  i <- matrixStats::colCumsums(matrix(1, nrow = nrow(Sigma_mu)))
  ij <- i %*% t(i)

  # Create m_1 based only on Sigma_mu
  # --> Do this first to create m_0 and m_01
  m_1 <- matrixStats::colCumsums(matrixStats::rowCumsums(Sigma_mu))
  m_1 <- m_1/ij

  # Reverse to create m_0
  m_0 <- apply(apply(m_1, 1, rev), 1, rev)

  # Compute covariance
  m_01 <- apply(matrixStats::colCumsums(matrixStats::rowCumsums(apply(Sigma_mu, 1, rev)))/ij, 1, rev)

  # Add Sigma_tau to m_1
  m_1_tau <- matrixStats::colCumsums(matrixStats::rowCumsums(Sigma_tau))
  m_1_tau <- m_1_tau/ij
  m_1 <- m_1 + m_1_tau

  # Return
  return_list <- list("m_0" = m_0, "m_1" = m_1, "m_01" = m_01)    #these are the moments that are used to find the posterior
  return(return_list)

}

##############################################################
# Function to estimate GP
##############################################################
posterior_predictions <- function(outcomes, error_var, cov, ix, prediction_grid, frequentist_uncertainty_only = FALSE) {

  # Re partition Sigma matrix to separate points we observe and points we don't
  S11 <- cov[ix, ix] + diag(error_var)
  S21 <- cov[, ix]
  S22 <- cov

  # Save inverse since using it twice
  S11_inv <- solve(S11)

  # Calculate posterior mean
  mu_hat <- S21 %*% S11_inv %*% outcomes

  # If focusing only on statistical uncertainty
  if (frequentist_uncertainty_only == TRUE) {
    Sigma_hat <- S21 %*% S11_inv %*% diag(error_var) %*% S11_inv %*% t(S21)
  }
  else {
    Sigma_hat <- S22 - S21 %*% S11_inv %*% t(S21)
  }

  # Force the Sigma_hat to be symmetric
  #--> Theoretical it is, but rounding errors in the inverse mess it up and it needs to be symmetric later
  Sigma_hat = (Sigma_hat + t(Sigma_hat))/2

  # Return
  return_list <- list("mu_hat" = mu_hat, "Sigma_hat" = Sigma_hat)
  return(return_list)

}

##############################################################
# Moment Transformation
##############################################################
moment_transformation <- function(p_lower, p_upper, mu_hat, Sigma_hat, eta) {

  # Get ix's
  lower_index <- eta %>% dplyr::filter(abs(eta - p_lower) < .0001)
  lower_index <- lower_index$row_n
  upper_index <- eta %>% dplyr::filter(abs(eta - p_upper) < .0001)
  upper_index <- upper_index$row_n
  ix <- matrix(c(lower_index, upper_index, lower_index + nrow(eta), upper_index + nrow(eta)))

  # Get weights
  weights <- matrix(c(-1*(1 - p_lower), (1 - p_upper), -1*p_lower, p_upper) )

  # Estimate
  MTE_hat <- (t(weights) %*% mu_hat[ix])/(p_upper - p_lower)
  MTE_hat_var <- (t(weights) %*% Sigma_hat[ix, ix] %*% weights)/((p_upper - p_lower)^2)

  # Return
  return_list <- list("mean" = MTE_hat, "variance" = MTE_hat_var)
  return(return_list)
}

##############################################################
# Marginal likelihood --> Used to estimate hyperparameters
##############################################################
log_m_likelihood <- function(hyperparams, data, eta, ml_only = TRUE, log_prior = NULL) {

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
  Sigma_mu <- (std_mu)^2*exp(-plgp::distance(eta$eta)/(2*lengthscale_mu^2))
  Sigma_tau <- (std_tau)^2*exp(-plgp::distance(eta$eta)/(2*lengthscale_tau^2))

  # Transform kernels
  m_covs <- transform_kernel(Sigma_mu, Sigma_tau)
  m_0 <- m_covs$m_0 + diag(eps, ncol(m_covs$m_0))
  m_1 <- m_covs$m_1 + diag(eps, ncol(m_covs$m_1))
  m_01 <- m_covs$m_01

  # Create big matrix
  Sigma <- cbind(rbind(m_0, m_01), rbind(t(m_01), m_1))

  # Relevant matrix
  K_y <- Sigma[matrix(data$row_n), matrix(data$row_n)] + diag(data$var)

  # Log Marginal Likelihood
  data_fit <- -1e8
  try(data_fit <- -.5*(t(matrix(data$mean)) %*% solve(K_y) %*% matrix(data$mean)))

  complexity_penalty <- -.5*log(det(K_y), base = exp(1))
  log_m_likelihood <- data_fit + complexity_penalty

  # If log prior
  if (!missing(log_prior)) {
    # Transform hyperparameters
    l_prior <- log_prior(hyperparams)
    log_m_likelihood <- log_m_likelihood + l_prior
  }

  # Return
  if (ml_only == TRUE) {
    return(log_m_likelihood)
  }
  else {
    if (!missing(log_prior)) {
      return_list <- list("data_fit" = data_fit, "complexity_penalty" = complexity_penalty, "log_m_likelihood" = log_m_likelihood, "log_prior" = l_prior)
      return(return_list)
    }
    else {
      return_list <- list("data_fit" = data_fit, "complexity_penalty" = complexity_penalty, "log_m_likelihood" = log_m_likelihood)
      return(return_list)
    }
  }
}


##############################################################
# Plot some stuff
##############################################################
plot_potential_functions <- function(mu, sigma, n = 10, poly_kernel = FALSE, poly_degree = 3) {

  # If poly_kernel == FALSE, can just simulate the data via the sigma matrix
  if (poly_kernel == FALSE) {
    Y <- mvtnorm::rmvnorm(n, mean = mu, sigma = sigma + diag(1*eps, nrow(mu)))
  }
  
  # If poly_kernel == TRUE, then sigma is low rank and so need to 
  # account for that
  if (poly_kernel == TRUE) {
    # Do an eigendecomposition and split out eigen vectors vs values
    max_component <- poly_degree*2 + 2
    eigen_decomp <- eigen(sigma)
    eigenvalues <- eigen_decomp$values[1:max_component]
    eigenvectors <- eigen_decomp$vectors[,1:max_component]
    
    # Draw random variables
    Y <- matrix(rnorm(n * dim(eigenvectors)[2]), ncol = dim(eigenvectors)[2]) %*% diag(sqrt(eigenvalues)) %*% t(eigenvectors)
  
    # Make mean matrix to same dimensions as Y
    adj_mu <- mu
    for (i in seq(1, n-1)) {
      adj_mu <- cbind(mu, adj_mu)
    }
    
    # Add mean
    Y <- Y + t(adj_mu)
  }
  
  # Save results
  top_n = ncol(Y)/2
  random_functions <- data.frame(
    x = rep(1:top_n, times = 2),
    f1 = Y[1,],
    f2 = Y[2,],
    f3 = Y[3,],
    f4 = Y[4,],
    f5 = Y[5,],
    f6 = Y[6,],
    f7 = Y[7,],
    f8 = Y[8,],
    f9 = Y[9,],
    f10 = Y[10,]
  )

  random_functions <- random_functions %>% dplyr::mutate(treatment = as.numeric(rownames(random_functions)) > max(random_functions$x))
  p <- ggplot2::ggplot(data = random_functions, aes(x = x/top_n)) + theme_bw() +
    ggplot2::geom_line(aes(y = f1, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f2, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f3, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f4, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f5, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f6, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f7, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f8, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f9, colour = treatment, linetype = treatment)) +
    ggplot2::geom_line(aes(y = f10, colour = treatment, linetype = treatment)) +
    ggplot2::xlab("Treatment Thresholds") + ggplot2::ylab("Conditional Average of Outcome")
  return(p)
}

##############################################################
# Calculate a Polynomial Kernel
##############################################################
poly_kernel <- function(hyperparams, degree, eta) {
  # Create Polynomial Matrices (with Intercept)
  a <- poly(eta$eta - .5, degree = degree, raw = TRUE)
  a <- cbind(1, a)
  
  # Penalty Matrix
  pen <- diag(hyperparams)
  
  # Kernel
  kernel <- a %*% pen %*% t(a)
  
  # Return
  return(kernel)
  
}


# ##############################################################
# # From Hypers to Squared Exponential Covariate
# ##############################################################
# squared_exponential <- function(hyperparams, eta) {
#   # From Vector of hyperparams to specific values
#   if (length(hyperparams) == 4) {
#     std_mu = exp(hyperparams[1])
#     lengthscale_mu = exp(hyperparams[2])
#     std_tau = exp(hyperparams[3])
#     lengthscale_tau = exp(hyperparams[4])
#   }
#   if (length(hyperparams) == 2) {
#     std_mu = exp(hyperparams[1])
#     lengthscale_mu = exp(hyperparams[2])
#     std_tau = exp(hyperparams[1])
#     lengthscale_tau = exp(hyperparams[2])
#   }
#   
#   # Create base kernels based on the hyperparams
#   Sigma_mu <- (std_mu)^2*exp(-distance(eta$eta)/(2*lengthscale_mu^2))
#   Sigma_tau <- (std_tau)^2*exp(-distance(eta$eta)/(2*lengthscale_tau^2))
#   
#   # Return
#   return(c("Sigma_mu" = Sigma_mu, "Sigma_tau" = Sigma_tau))
# }

##############################################################
# Bayesian MTE -- Subroutine
##############################################################
bayesian_mte_sub <- function(Sigma_mu, Sigma_tau, merged, eta, c_data, grid_length, frequentist_uncertainty = FALSE)  {

  # Transform kernels
  m_covs <- transform_kernel(Sigma_mu, Sigma_tau)
  m_0 <- m_covs$m_0 + diag(eps, ncol(m_covs$m_0))
  m_1 <- m_covs$m_1 + diag(eps, ncol(m_covs$m_1))
  m_01 <- m_covs$m_01

  # Create big matrix
  Sigma <- cbind(rbind(m_0, m_01), rbind(t(m_01), m_1)) #other code plots potential prior functions after this

  # Posterior predictions about m
  posterior_m <- posterior_predictions(merged$mean, merged$var, Sigma, matrix(merged$row_n), eta, frequentist_uncertainty_only = frequentist_uncertainty)
  mu_hat <- posterior_m$mu_hat
  Sigma_hat <- posterior_m$Sigma_hat

  #other code plots potential posterior functions as well as plot observed moments and predictions

  # Calculate Weights
  eta_length <- nrow(mu_hat)/2
  w <- matrix(0, nrow = nrow(mu_hat)/2 - 1, ncol = nrow(mu_hat))
  for (i in seq(1,eta_length - 1,1)) {
    w[i,i] <- -(eta_length-i)
    w[i,i+1] <- (eta_length-i-1)
    w[i,eta_length+i] <- -(i-1)
    w[i,eta_length+i+1] <- i
  }
  Sigma_tau_hat <- w %*% Sigma_hat %*% t(w)
  tau_hat <- w %*% mu_hat

  # ATE
  ATE_transform <- moment_transformation(0, 1, mu_hat, Sigma_hat, eta)

  # LATE
  p_lower <- min(c_data$eta)
  p_upper <- max(c_data$eta)
  LATE_transform <- moment_transformation(p_lower, p_upper, mu_hat, Sigma_hat, eta)

  # Average on Always Takers
  ATATE_transform <- moment_transformation(0, p_lower, mu_hat, Sigma_hat, eta)

  # Average on Never Takers
  NTATE_transform <- moment_transformation(p_upper, 1, mu_hat, Sigma_hat, eta)

  return_list <- list("MTE_hat" = tau_hat, "var_MTE_hat" = Sigma_tau_hat, "ATE" = ATE_transform, "LATE" = LATE_transform, "ATATE" = ATATE_transform, "NTATE" = NTATE_transform)
  return(return_list)

}


##############################################################
##############################################################
# Main Function
##############################################################
##############################################################
#' implements the Bayesian marginal treatment effect approach described in Opper (2022).
#'
#'
#' @param outcome is the outcome; input should be a vector.
#' @param treatment_choice is the treatment status; input should be a vector.
#' @param treatment_assignment is the treatment assignment; input should be a vector.
#' @param grid_length is the grid length of eta; default is 0.01.
#' @param demean determines whether to remove the mean from the outcome before the analysis; default is FALSE.
#' @param full_bayes determines to integrate over the prior distribution of hyperparameters (if TRUE) or to estimate the hyperparameters (if FALSE); default is true.
#' @param input_hypers hyperparameters; default is NULL in which they are estimated.
#' @param log_prior is a function that specifies prior -- see documentation for description of function; default is NULL.
#' @param hyperparameter_draws is number of draws from the hyper-prior distribution; default is 10,000.
#' @param extrapolation_uncertainty determines whether the observed moments are assumed to be estimated without error (if TRUE) or with error (if FALSE);  default is FALSE.
#' @param frequentist_uncertainty determines whether the reported standard deviations are the standard errors of the MAP (if TRUE) or the posterior standard deviation (if FALSE);  default is FALSE.
#' @param poly_kernel determines whether to use a polynomial kernel (if TRUE) as opposed a squared-exponential kernel (if FALSE). Default is FALSE, i.e., to use a squared expoential. If using a polynomial kernel, you need to specify an input_hypers of size degree + 1, which specify the prior variance of each coefficient on the polynomials, starting with the constant term and cannot specify full_bayes == TRUE.
#' @param poly_degree determines the degree of the polynomial to use if specifying a polynomial kernel. Default is 3.
#' @return ATE is the posterior mean (ATE$mean) and variance (ATE$variance) of the average treatment effect
#' @return LATE, ATATE, NTATE are identical to ATE, but for the local average treatment effect (LATE), always taker average treatment effect (ATATE), and never taker average treatment effect (NTATE)
#' @return prior_plot, posterior_plot are plots showing 10 random functions generated according to the estimated hyperparameter. Prior plot does not condition on the observed moments, while posterior plot does. Only returned when full_bayes == FALSE.
#' @return MTE, MTE_plot are the posterior mean/variance of the marginal treatment effect function at each point on the eta-grid. MTE is a matrix with the information and MTE_plot returns a graph with the posterior mean and 95% CI. Only returned when full_bayes == FALSE.
#' @return predictions_plot is a plot showing the posterior mean and 95% CI of the m_1 and m_0 functions, as defined in Opper (2022), along with the observed moments. Only returned when full_bayes == FALSE.
#' @return hypers is a list of the hyperparameters that were used. Only returned when full_bayes == FALSE.
#' @return posteriors_and_priors is a graph showing the prior and posterior distribution of the hyperparameters. Only returned when full_bayes == TRUE.
#' @return number_of_hypers_used is the number of hyperparameters that are "accepted" by the accept/reject algorithm. Only returned when full_bayes == TRUE.
#' @export
bayesian_mte <- function(outcome, treatment_choice, treatment_assignment, grid_length = .01, demean = FALSE, full_bayes = TRUE, input_hypers = NULL, log_prior = NULL, hyperparameter_draws = 10000, extrapolation_uncertainty = FALSE, frequentist_uncertainty = FALSE, poly_kernel = FALSE, poly_degree = 3) {

  ##############################################################
  # Error Checks
  ##############################################################
  # outcome, treatment_choice, and treatment_assignment are all inputted as vectors - list of numbers
  stopifnot("outcome, treatment_choice, and treatment_assignment must all be inputted as vectors" = (length(treatment_choice) > 1) & (length(treatment_assignment) > 1) & (length(outcome) > 1 & all(sapply(outcome, is.numeric))))

  # log_prior is a function that both returns the log_prior if given a set of hyperparameters and randomly draws from the prior when specified that draw = TRUE
  if (!is.null(log_prior)) {
    #ensure log_prior has draws input and hypers input (with default null)
    stopifnot("log_prior must have two inputs: a variable named \"hypers\" with default NULL, and a Boolean variable named \"draw\"", "draws" %in% methods::formalArgs(log_prior) & "hypers" %in% methods::formalArgs(log_prior))
  }

  # Can't specify that extrapolation_uncertainty and frequentist_uncertainty are both true
  stopifnot("extrapolation_uncertainty and frequentist_uncertainty cannot both be TRUE" = !(extrapolation_uncertainty == TRUE & frequentist_uncertainty == TRUE))


  ##############################################################
  # Collapse data to moments
  ##############################################################
  # Combine into new dataframe
  base_data <- data.frame(
    outcome = outcome,
    treatment = treatment_choice,
    IV = treatment_assignment
  )
  #specify column names
  colnames(base_data) = c('outcome', 'treatment', 'IV')
  # Drop any observations that are missing outcomes, treatment, or IV
  base_data = base_data[stats::complete.cases(base_data), ]

  # Demean Outcomes if asked to
  if (demean == TRUE) {
    base_data$outcome <- base_data$outcome - mean(base_data$outcome)
  }

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

  # Zero out variance if want to focus only on extrapolation uncertainty
  if (extrapolation_uncertainty == TRUE) {
    c_data <- c_data %>% dplyr::mutate(var = 1e-8)
  }

  # Create Grid
  eta <- tibble::tibble(eta = seq(0, 1, by= grid_length), row_n = seq(1, 1/grid_length + 1, 1))

  # Round p_score to merge to eta_grid
  c_data <- c_data %>% dplyr::mutate(eta = round(p_score/grid_length)*grid_length)

  # Get row_number of prediction data --> Note, since treat are stacked below controls, treat row is order in the sequence + nrow(prediction_grid)
  merged_control <- merge(c_data %>% dplyr::filter(treatment == 0), eta) %>% dplyr::select(mean, var, row_n)
  merged_treat <- merge(c_data %>% dplyr::filter(treatment == 1), eta) %>% dplyr::select(mean, var, row_n) %>% dplyr::mutate(row_n = row_n + nrow(eta))
  merged <- rbind(merged_control, merged_treat)

  ##############################################################
  # If using the default prior, define the log_prior function
  ##############################################################
  if (is.null(log_prior)) {
    # Mean/Variance Parameters
    #corr_mean <- .5
    #l_mean <- (-1/log(corr_mean))^.5
    l_mean <- 1
    l_std <- .5
    std_mean <- .5*residuals_std
    std_std <- 1.25

    # Define the log_prior using that info
    log_prior <- function(hypers = NULL, draws = FALSE) {
      if (draws == FALSE) {
        #
        adj_hypers <- hypers

        # Adjust first and third by simply adjusting mean & standard deviation
        adj_hypers[1] <- (hypers[1] - log(std_mean))/(std_std)
        adj_hypers[3] <- (hypers[1] - log(std_mean))/(std_std)


        # Adjust second and fourth
        adj_hypers[2] <- (hypers[2] - log(l_mean))/(l_std)
        adj_hypers[4] <- (hypers[2] - log(l_mean))/(l_std)

        # Return log_prior
        l_prior <- log(stats::dnorm(adj_hypers[1])) + log(stats::dnorm(adj_hypers[2])) + log(stats::dnorm(adj_hypers[3])) + log(stats::dnorm(adj_hypers[4]))
        return(l_prior)
      }
      else {
        # Start with draw from standard normal
        #hypers <- stats::rnorm(4)
        hyper_draws <- stats::rnorm(2)
        hypers <- c(hyper_draws[1], hyper_draws[2], hyper_draws[1], hyper_draws[2])

        # Adjust the first and third, so the standard deviation parameters, i.e., exp(hypers[1]), is distributed log-normal with mean std_mean and std of std_std
        hypers[1] <- (hypers[1]*(std_std) + log(std_mean))
        hypers[3] <- (hypers[3]*(std_std) + log(std_mean))

        # Turn from correlation parameters to log-lengthscale parameters
        hypers[2] <- (hypers[2]*(l_std) + log(l_mean))
        hypers[4] <- (hypers[4]*(l_std) + log(l_mean))

        return(hypers)
      }
    }
  }

  ##############################################################
  # If going the empirical Bayes route, as opposed to drawing from the hyperparameters
  ##############################################################
  if (full_bayes == FALSE) {
    ##############################################################
    # If needed, estimated hyperparameters
    ##############################################################
    if (missing(input_hypers)) {
      ml_hyper <- stats::optim(c(0,0), log_m_likelihood, gr = NULL, merged, eta, log_prior = log_prior, control = list(fnscale = -1))
      ml_hyper <- ml_hyper$par
    }
    else {
      ml_hyper = input_hypers
    }
    
    ##############################################################
    # Create Kernels
    ##############################################################
    if (poly_kernel == FALSE) {
      # From Vector of hyperparams to specific values
      if (length(ml_hyper) == 4) {
        std_mu = exp(ml_hyper[1])
        lengthscale_mu = exp(ml_hyper[2])
        std_tau = exp(ml_hyper[3])
        lengthscale_tau = exp(ml_hyper[4])
      }
      if (length(ml_hyper) == 2) {
        std_mu = exp(ml_hyper[1])
        lengthscale_mu = exp(ml_hyper[2])
        std_tau = exp(ml_hyper[1])
        lengthscale_tau = exp(ml_hyper[2])
      }

      # Create base kernels based on the hyperparams
      Sigma_mu <- (std_mu)^2*exp(-distance(eta$eta)/(2*lengthscale_mu^2))
      Sigma_tau <- (std_tau)^2*exp(-distance(eta$eta)/(2*lengthscale_tau^2))
    }
    if (poly_kernel == TRUE) {
      Sigmas <- poly_kernel(ml_hyper, degree = poly_degree, eta)
      Sigma_mu <- Sigmas
      Sigma_tau <- Sigmas
    }

    ##############################################################
    # Bayesian MTE -- Subroutine
    ##############################################################
    # Transform kernels
    m_covs <- transform_kernel(Sigma_mu, Sigma_tau)
    m_0 <- m_covs$m_0 + diag(eps, ncol(m_covs$m_0))
    m_1 <- m_covs$m_1 + diag(eps, ncol(m_covs$m_1))
    m_01 <- m_covs$m_01

    # Create big matrix
    Sigma <- cbind(rbind(m_0, m_01), rbind(t(m_01), m_1))

    # Plot potential prior functions
    prior_plot <- plot_potential_functions(matrix(0, nrow = nrow(Sigma)), Sigma, poly_kernel = poly_kernel, poly_degree = poly_degree)
    
    # Posterior predictions about m
    posterior_m <- posterior_predictions(merged$mean, merged$var, Sigma, matrix(merged$row_n), eta, frequentist_uncertainty_only = frequentist_uncertainty)
    mu_hat <- posterior_m$mu_hat
    Sigma_hat <- posterior_m$Sigma_hat

    # Plot potential posterior functions
    posterior_plot <- plot_potential_functions(mu_hat, Sigma_hat, poly_kernel = poly_kernel, poly_degree = poly_degree)
    
    # Plot Observed Moments & Predictions
    # Make predictions dataframe
    predictions <- data.frame(
      mu_hat = mu_hat,
      var = diag(Sigma_hat)
    )

    # Add pscore values
    predictions <- cbind(predictions, rbind(eta, eta))

    # Split treatment/control
    predictions <- predictions %>% dplyr::mutate(treatment = as.numeric(rownames(predictions)) > max(predictions$row_n))

    # Add CIs
    predictions <- predictions %>% dplyr::mutate(p5 = mu_hat - 1.96*var^.5, p95 = mu_hat + 1.96*var^.5)

    # Plot
    predictions_plot <- ggplot2::ggplot() + ggplot2::theme_bw() +
      ggplot2::geom_line(data = predictions, mapping = aes(x = eta, y = mu_hat, colour = treatment, linetype = treatment)) +
      ggplot2::geom_line(data = predictions, mapping = aes(x = eta, y = p5, colour = treatment, linetype = treatment)) +
      ggplot2::geom_line(data = predictions, mapping = aes(x = eta, y = p95, colour = treatment, linetype = treatment)) +
      ggplot2::geom_point(data = c_data, mapping = aes(x = p_score, y = mean)) +
      ggplot2::xlab("Treatment thresholds") + ggplot2::ylab("Conditional average outcome")

    # From m to tau
    eta_length <- nrow(mu_hat)/2
    w <- matrix(0, nrow = nrow(mu_hat)/2 - 1, ncol = nrow(mu_hat))
    for (i in seq(1,eta_length - 1,1)) {
      w[i,i] <- -(eta_length-i)
      w[i,i+1] <- (eta_length-i-1)
      w[i,eta_length+i] <- -(i-1)
      w[i,eta_length+i+1] <- i
    }
    Sigma_tau_hat <- w %*% Sigma_hat %*% t(w)
    std_tau_hat <- diag(Sigma_tau_hat)^.5
    tau_hat <- w %*% mu_hat
    MTEs <- tibble::tibble(tau_hat, std_tau_hat)

    # ATE
    ATE_transform <- moment_transformation(0, 1, mu_hat, Sigma_hat, eta)

    # LATE
    p_lower <- min(c_data$eta)
    p_upper <- max(c_data$eta)
    LATE_transform <- moment_transformation(p_lower, p_upper, mu_hat, Sigma_hat, eta)

    # Average on Always Takers
    ATATE_transform <- moment_transformation(0, p_lower, mu_hat, Sigma_hat, eta)

    # Average on Never Takers
    NTATE_transform <- moment_transformation(p_upper, 1, mu_hat, Sigma_hat, eta)

    ##############################################################
    # Plot MTEs
    ##############################################################
    MTEs <- MTEs %>% dplyr::mutate(q5 = tau_hat - 1.96*std_tau_hat, q95 = tau_hat + 1.96*std_tau_hat) %>%
      dplyr::mutate(eta = seq(0, 1 - .01, .01))

    MTE_plot <- ggplot2::ggplot(MTEs) + ggplot2::theme_bw() +
      ggplot2::geom_line(aes(eta, tau_hat)) +
      ggplot2::geom_line(aes(eta, q95), linetype = 'dashed') +
      ggplot2::geom_line(aes(eta, q5), linetype = 'dashed') +
      ggplot2::xlab(expression(eta)) + ggplot2::ylab("Conditional average treatment effect")

    ##############################################################
    # Prior vs. Posterior ATE/LATE/ATATE/NTATE Distributions
    ##############################################################
    # Calculate prior variances based on Sigma_tau
    # Get ix's
    lower_index <- eta %>% dplyr::filter(abs(eta - p_lower) < .0001)
    lower_index <- lower_index$row_n
    upper_index <- eta %>% dplyr::filter(abs(eta - p_upper) < .0001)
    upper_index <- upper_index$row_n

    ATE_weight <- matrix(1/101, nrow = 101, ncol = 1)
    LATE_weight <- matrix(0, nrow = 101, ncol = 1)
    LATE_weight[lower_index:upper_index, 1] <- 1/(upper_index - lower_index + 1)
    ATATE_weight <- matrix(0, nrow = 101, ncol = 1)
    ATATE_weight[1:lower_index, 1] <- 1/(lower_index)
    NTATE_weight <- matrix(0, nrow = 101, ncol = 1)
    NTATE_weight[upper_index:101, 1] <- 1/(101-upper_index + 1)

    ATE_prior_var <- t(ATE_weight) %*% Sigma_tau %*% ATE_weight
    LATE_prior_var <- t(LATE_weight) %*% Sigma_tau %*% LATE_weight
    ATATE_prior_var <- t(ATATE_weight) %*% Sigma_tau %*% ATATE_weight
    NTATE_prior_var <- t(NTATE_weight) %*% Sigma_tau %*% NTATE_weight

    # Plot prior and posteriors
    ATE_graph <- ggplot(data = data.frame(prior = rnorm(100000, 0, ATE_prior_var^.5),posterior = rnorm(100000,ATE_transform$mean, ATE_transform$variance^.5))) +
      geom_density(aes(x = prior), fill = 'blue', alpha = 0.3) +
      geom_density(aes(x = posterior), fill = 'red', alpha = 0.3) +
      labs(x = "Treatment Effect", y = "Density") + theme_bw()
    LATE_graph <- ggplot(data = data.frame(prior = rnorm(100000, 0, LATE_prior_var^.5),posterior = rnorm(100000,LATE_transform$mean, LATE_transform$variance^.5))) +
      geom_density(aes(x = prior), fill = 'blue', alpha = 0.3) +
      geom_density(aes(x = posterior), fill = 'red', alpha = 0.3) +
      labs(x = "Treatment Effect", y = "Density") + theme_bw()
    NTATE_graph <- ggplot(data = data.frame(prior = rnorm(100000, 0, NTATE_prior_var^.5),posterior = rnorm(100000,NTATE_transform$mean, NTATE_transform$variance^.5))) +
      geom_density(aes(x = prior), fill = 'blue', alpha = 0.3) +
      geom_density(aes(x = posterior), fill = 'red', alpha = 0.3) +
      labs(x = "Treatment Effect", y = "Density") + theme_bw()
    ATATE_graph <- ggplot(data = data.frame(prior = rnorm(100000, 0, ATATE_prior_var^.5),posterior = rnorm(100000,ATATE_transform$mean, ATATE_transform$variance^.5))) +
      geom_density(aes(x = prior), fill = 'blue', alpha = 0.3) +
      geom_density(aes(x = posterior), fill = 'red', alpha = 0.3) +
      labs(x = "Treatment Effect", y = "Density") + theme_bw()

    # Add prior vars to ATE transform
    ATE_transform$prior_var <- ATE_prior_var
    ATE_transform$graph <- ATE_graph
    LATE_transform$prior_var <- LATE_prior_var
    LATE_transform$graph <- LATE_graph
    ATATE_transform$prior_var <- ATATE_prior_var
    ATATE_transform$graph <- ATATE_graph
    NTATE_transform$prior_var <- NTATE_prior_var
    NTATE_transform$graph <- NTATE_graph

    ##############################################################
    # Return
    ##############################################################
    return_list <- list("MTE" = MTEs, "full_post_var" = Sigma_tau_hat, "ATE" = ATE_transform, "LATE" = LATE_transform, "ATATE" = ATATE_transform, "NTATE" = NTATE_transform,  "prior_plot" = prior_plot, "posterior_plot" = posterior_plot, "predictions_plot" = predictions_plot, "MTE_plot" = MTE_plot, "hypers" = ml_hyper)
    return(return_list)
  }

  ##############################################################
  # If going the fully Bayesian route
  ##############################################################
  if (full_bayes == TRUE) {

    # Set up results matrix
    ml_results <- matrix(0, hyperparameter_draws, 5)

    # Loop
    for (k in seq(1, hyperparameter_draws)) {

      # Randomly draw hyperparameters
      hypers <- log_prior(draw = TRUE)

      # For each set of hypers estimate marginal likelihood
      ll <- log_m_likelihood(hypers, merged, eta, ml_only = TRUE)

      # Results
      ml_results[k,1] <- hypers[1]
      ml_results[k,2] <- hypers[2]
      ml_results[k,3] <- hypers[1]
      ml_results[k,4] <- hypers[2]
      ml_results[k,5] <- ll[1]
    }

    # Turn into dataframe
    ml_results <- data.frame(ml_results) %>% dplyr::rename(std_mu = X1, lengthscale_mu = X2, std_tau = X3, lengthscale_tau = X4, log_ml = X5)

    # From log likelihood to likelihood
    ml_results <- ml_results %>% dplyr::mutate(ml = exp(log_ml))
    ml_results <- ml_results %>% dplyr::mutate(ml = ml/max(ml))

    # Plot Prior vs Posterior Distributions
    length_mu <- ggplot2::ggplot(ml_results, aes(lengthscale_mu)) + ggplot2::theme_bw() +
      ggplot2::geom_histogram(color="black", alpha = .4, bins = 30) + ggplot2::geom_histogram(aes(weight = ml/mean(ml)), bins = 30, color = 'blue', alpha=.8, fill = "lightblue") +
      ggplot2::xlab("Lengthscale mu") + ggplot2::ylab(" ")
    length_tau <- ggplot2::ggplot(ml_results, aes(lengthscale_tau)) + ggplot2::theme_bw() +
      ggplot2::geom_histogram(color="black", alpha = .4, bins = 30) + ggplot2::geom_histogram(aes(weight = ml/mean(ml)), bins = 30, color = 'blue', alpha=.8, fill = "lightblue") +
      ggplot2::xlab("Lengthscale tau") + ggplot2::ylab(" ")
    std_tau <- ggplot2::ggplot(ml_results, aes(std_tau)) + ggplot2::theme_bw() +
      ggplot2::geom_histogram(color="black", alpha = .4, bins = 30) + ggplot2::geom_histogram(aes(weight = ml/mean(ml)), bins = 30, color = 'blue', alpha=.8, fill = "lightblue") +
      ggplot2::xlab("Standard Deviation tau") + ggplot2::ylab(" ")
    std_mu <- ggplot2::ggplot(ml_results, aes(std_mu)) + ggplot2::theme_bw() +
      ggplot2::geom_histogram(color="black", alpha = .4, bins = 30) + ggplot2::geom_histogram(aes(weight = ml/mean(ml)), bins = 30, color = 'blue', alpha=.8, fill = "lightblue") +
      ggplot2::xlab("Standard Deviation mu") + ggplot2::ylab(" ")
    posteriors_and_priors <- (length_mu + length_tau ) / (std_mu + std_tau)
    posteriors_and_priors + patchwork::plot_annotation(
      title = "Prior and posterior distributions of hyperparameters",
      caption = 'Blue = posterior distribution; black = prior distribution'
    )

    # Accept/Reject Sampling
    ml_results <- ml_results %>% dplyr::mutate(runif = stats::runif(hyperparameter_draws))  %>% dplyr::mutate(keep = (ml > runif))

    # Turn hyperparameters into results
    ml_results <- ml_results %>% dplyr::filter(keep == 1)
    num = nrow(ml_results)
    param_results <- matrix(0, num, 12)
    MTE_results <- data.frame(matrix(0,num,nrow(eta)-1))
    for (k in seq(1, num)) {
      # Adjust Hyperparams
      std_mu = exp(ml_results[k,1])
      lengthscale_mu = exp(ml_results[k,2])
      std_tau = exp(ml_results[k,3])
      lengthscale_tau = exp(ml_results[k,4])
      
      # Calculate Kernels
      Sigma_mu <- (std_mu)^2*exp(-distance(eta$eta)/(2*lengthscale_mu^2))
      Sigma_tau <- (std_tau)^2*exp(-distance(eta$eta)/(2*lengthscale_tau^2))
      
      # Calculate Treatment
      treat_effects <- bayesian_mte_sub(Sigma_mu, Sigma_tau, merged, eta, c_data, grid_length, frequentist_uncertainty = frequentist_uncertainty)

      # Results
      param_results[k,1] <- ml_results[k,1]
      param_results[k,2] <- ml_results[k,2]
      param_results[k,3] <- ml_results[k,3]
      param_results[k,4] <- ml_results[k,4]
      param_results[k,5] <- treat_effects$ATE$mean
      param_results[k,6] <- treat_effects$ATE$variance
      param_results[k,7] <- treat_effects$LATE$mean
      param_results[k,8] <- treat_effects$LATE$variance
      param_results[k,9] <- treat_effects$ATATE$mean
      param_results[k,10] <- treat_effects$ATATE$variance
      param_results[k,11] <- treat_effects$NTATE$mean
      param_results[k,12] <- treat_effects$NTATE$variance

      # MTE Results
      MTE_results[k,] <- treat_effects$MTE_hat
      if (k == 1) {
        var_MTE_hat <- treat_effects$var_MTE_hat
      }
      else {
        var_MTE_hat <- ((k-1)*var_MTE_hat + treat_effects$var_MTE_hat)/k
      }
    }

    # Combine to get MTE results
    Sigma_tau_hat <- var_MTE_hat + cov(MTE_results)
    std_tau_hat <- diag(Sigma_tau_hat)^.5
    tau_hat <- colMeans(MTE_results)
    MTEs <- tibble::tibble(tau_hat, std_tau_hat)

    # Save MTEs in a tibble
    MTEs <- MTEs %>% dplyr::mutate(q5 = tau_hat - 1.96*std_tau_hat, q95 = tau_hat + 1.96*std_tau_hat) %>%
      dplyr::mutate(eta = seq(0, 1 - .01, .01))

    # Save MTE_plot
    MTE_plot <- ggplot2::ggplot(MTEs) + ggplot2::theme_bw() +
      ggplot2::geom_line(aes(eta, tau_hat)) +
      ggplot2::geom_line(aes(eta, q95), linetype = 'dashed') +
      ggplot2::geom_line(aes(eta, q5), linetype = 'dashed') +
      ggplot2::xlab(expression(eta)) + ggplot2::ylab("Conditional average treatment effect")

    # Make param_results into a data.frame with correct variable names
    full_results <- data.frame(param_results) %>% dplyr::rename(sigma_mu = X1, lengthscale_mu = X2, sigma_tau = X3, lengthscale_tau = X4, ATE_mean = X5, ATE_variance = X6,
                                                                LATE_mean = X7, LATE_variance = X8, ATATE_mean = X9, ATATE_variance = X10, NTATE_mean = X11, NTATE_variance = X12)

    # Turn mixtures of Gaussians into mean/std of effects
    return_list <- list("MTE" = MTEs, "full_post_var" = Sigma_tau_hat, "ATE" = list("mean" = mean(param_results[1:num,5]), "variance" = sd(param_results[1:num,5])^2 + mean(param_results[1:num,6])), "LATE" = list("mean" = mean(param_results[1:num,7]), "variance" = stats::sd(param_results[1:num,7])^2 + mean(param_results[1:num,8])), "ATATE" = list("mean" = mean(param_results[1:num,9]), "variance" = stats::sd(param_results[1:num,9])^2 + mean(param_results[1:num,10])), "NTATE" = list("mean" = mean(param_results[1:num,11]), "variance" = stats::sd(param_results[1:num,11])^2 + mean(param_results[1:num,12])), "posteriors_and_priors" = posteriors_and_priors, "number_of_hypers_used" = num, 'full_results' = full_results, "MTE_plot" = MTE_plot)
    return(return_list)

  }

}
