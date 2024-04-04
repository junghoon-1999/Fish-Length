### Main script of MT4113 Assignment 2
# Contiains teamEM and all other relevant functions

# Load required libraries
library(dplyr)

# Function to be called by the user. Implements and returns the result of 
# the EM algorithm
#
# INPUTS:
#  data: df with columns FishID, Length, Age
#  epsilon: tolerance value [default 1e-08]
#  maxit: maximum number of iterations [default 1000]
#
# OUTPUTS:
#  estimates: estimates df; rows Age1, Age2, Age3; columns mu, sigma, and lambda
#  inits: initial value df; rows Age1, Age2, Age3; columns mu, sigma, and lambda
#  converged: Boolean value showing if algorithm converged
#  posterior: df of posterior probabilities; N rows; k columns
#  likelihood: log-likelihood value vector; length is no. of iterations
teamEM <- function(data, epsilon = 1e-08, maxit = 1000) {
    
  # Check user input for validity
  if(!is.data.frame(data) | 
      !all(sapply(colnames(data), 
                  function(x) x %in% c("FishID", "Length", "Age"))) |
      maxit < 0 |
      maxit != round(maxit) |
      !is.numeric(maxit) |
      epsilon < 0 |
      !is.numeric(epsilon)) {
      stop("invalid argument")
  }
  
  # Initialise variables
  inits  <-  initialisation(data)
  estimates <- inits$estimates
  likelihood  <-  c(inits$likelihood)
  converged  <-  FALSE
  iteration_no  <-  1

  # Run the algorithm until it has converged or maxit is reached
  while (!converged) {
      posterior <- expectation(data, estimates)
      estimates <- maximisation(data, posterior, estimates)
      
      currLikelihood <- logLikelihood(estimates, data)
      
      # Check if algorithm has converged
      if (abs(currLikelihood - likelihood[length(likelihood)]) < epsilon){
        converged <- TRUE
      }
      
      likelihood <- c(likelihood, currLikelihood)

      # append prevLikelihood to likelihood
      # likelihood <- c(likelihood, convergenceResult[["prevLikelihood"]])

      # Check if maxit has been reached
      if (iteration_no >= maxit) {
          break
      } else {
          iteration_no = iteration_no + 1
      }
  }

  # Create output list
  output <- list(estimates = estimates,
                 inits = inits,
                 converged = converged,
                 posterior = posterior,
                 likelihood = likelihood)

    return(output)
}


# Initialises values for the algorithm
#
# INPUTS
#  data: df with columns FishID, Length, Age
#
# OUPUTS:
#  estimates: initial estimates df; rows Age1, Age2, Age3; 
#     columns mu, sigma, lambda
#  likelihood: initial likelihood
initialisation <- function(data) {
    estimates  <-  NULL
    known_age_data  <-  data %>% filter(!is.na(Age))
    unknown_age_data  <-  data %>% filter(is.na(Age))

    # determine group mean of know age data
    mu1 <- mean((known_age_data %>% filter(Age == 1))$Length)
    mu2 <- mean((known_age_data %>% filter(Age == 2))$Length)
    mu3 <- mean((known_age_data %>% filter(Age == 3))$Length)
    
    # Make initial group guess for each fish with unknown age
    unknown_age_data$Age <- sapply(unknown_age_data$Length, 
                                   FUN = function(x) which.min(c(abs(x - mu1), 
                                                                 abs(x - mu2), 
                                                                 abs(x - mu3))))

    # Reombine the data frames
    data <- rbind(known_age_data, unknown_age_data)
    
    # Create empty data frame for estimates
    estimates <- data.frame(matrix(ncol = 3, nrow = 3))
    colnames(estimates) <- c("mu", "sigma", "lambda")

    # Make estimates for each age group
    for (k in c(1, 2, 3)) {
        groupData <- data %>% filter(Age == k)
        mu <-  mean(groupData$Length)
        sigma <- sd(groupData$Length)
        lambda <- nrow(groupData) / nrow(data)
        estimates[k,] <- c(mu, sigma, lambda)
    }

    # Format data frame
    rownames(estimates) <-  c("Age1", "Age2", "Age3")
    colnames(estimates) <- c("mu", "sigma", "lambda")

    # Calculate current expectation
    posterior  <-  expectation(data, estimates)
    
    # Perform maximisation based on current expectation
    estimates  <-  maximisation(data, posterior, estimates)

    # Calculate current likelihood
    likelihood <- logLikelihood(estimates, data)
    
    # Create output
    output <- list(estimates = estimates, 
                   likelihood = likelihood)

    return(output)
}


# Calculate the log likelihood for data in the standard project format
#
# INPUTS:
#  data: df with columns FishID, Length, Age
#  estimates: estimates df; rows Age1, Age2, Age3; columns mu, sigma, and lambda
#
# OUTPUTS:
#  logL: Log Likelihood (logarithmised)
logLikelihood <- function(estimates, data) {
  
  # Setup initial values
  lambdas <- estimates$lambda
  sigmas <- estimates$sigma
  mus <- estimates$mu
  K <- length(lambdas)

  # Setup empty matrix of likelihoods
  Lik <- matrix(NA, nrow = length(data$Length), ncol = K)
  
  # Iterate over all groups
  for (k in 1:K) {
    Lik[, k] <- dnorm(data$Length, mean = mus[k], sd = sigmas[k]) * lambdas[k]  
  }
  
  # calculate log likelihood as the sum of logarithmised row sums of Lik
  logL <- sum(log(rowSums(Lik)))
  
  # Note that this is the lograrithm of the variable we want for accuracy
  return(logL)
}


# Calculate the probability of each observation belonging to each group
#
# INPUTS:
#  data: dataframe with columns FishID, Length, Age
#  estimates: dataframe of estimates with rows Age1, Age2,
#    Age3 and columns mu, sigma, and lambda
#
# OUTPUTS:
#   posterior: dataframe of posterior probabilities with N
#    (= number of observations) rows and k (= number of components) columns

expectation <- function(data,estimates) {
  
  # Compute the Gauss densities
  gauss_dens1 <- dnorm(data$Length, estimates["Age1", "mu"], 
                       estimates["Age1","sigma"])
  gauss_dens2 <- dnorm(data$Length, estimates["Age2","mu"], 
                       estimates["Age2","sigma"])
  gauss_dens3 <- dnorm(data$Length, estimates["Age3","mu"], 
                       estimates["Age3","sigma"])
  
  # Calculate P(y_i)
  pyi <- gauss_dens1 * estimates["Age1","lambda"] + 
    gauss_dens2 * estimates["Age2","lambda"] + 
    gauss_dens3 * estimates["Age3","lambda"]
  
  # Create data frame of expectations per Age group
  posterior = data.frame(Age1 = gauss_dens1 * estimates["Age1","lambda"] / pyi, 
                         Age2 = gauss_dens2 * estimates["Age2","lambda"] / pyi, 
                         Age3 = gauss_dens3 * estimates["Age3","lambda"] / pyi)
  return(posterior)
}


# Update best estimates of the means and variances given the current grouping
#
# INPUTS:
#   data: dataframe with columns FishID, Length, Age
#   posterior: dataframe of posterior probabilities with N
#   (= number of observations) rows and k (= number of components) columns
#   estimates: dataframe of estimates with rows Age1, Age2,
#   Age3 and columns mu, sigma, and lambda
#
# OUTPUTS:
#   estimates: dataframe of estimates with rows Age1, Age2,
#   Age3 and columns mu, sigma, and lambda

maximisation <- function(data, posterior, estimates) {
  
  # Iterate over all groups
  for (k in 1:3) {
    
    numerator_mu <- sum(posterior[, k] * data$Length)
    
    sum_posterior_values <- sum(posterior[, k])
    
    mu <- numerator_mu / sum_posterior_values
    
    numerator_sigma <- sum(posterior[, k] * (data$Length - mu) ^ 2)
    
    sigma <- sqrt(numerator_sigma / sum_posterior_values)
    
    lambda <- sum_posterior_values / nrow(posterior)
    
    # Update estimates
    estimates[k,] <- c(mu, sigma, lambda)
  
  }
  return(estimates)
}
