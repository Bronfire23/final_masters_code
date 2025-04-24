#####################File info#####################
# Code that stores all functions used by all models.
# The functions created is:
#     1. calculate_data_belongings
#     2. calculate_node_local_sufficient_statistics
#     3. estimate_parameters_with_sufficient_statistics
#     4. calculate_log_likelihood_with_full_data
#     5. calculate_data_belongings_multivariate
#     6. calculate_node_local_sufficient_statistics_multivariate
#     7. estimate_parameters_with_sufficient_statistics_multivariate
#     8. calculate_log_likelihood_with_full_data_multivariate

# Libraries
library(mvtnorm)

calculate_data_belongings <- function(params, data){
  ### Extracting parameters and checking if they are they same.
  means <- params$means
  stds <- params$stds
  probs <- params$mixing_probs
  if (length(means) != length(stds) || length(means) != length(probs)){
    stop("Number of means, standard deviations, and weights must be equal.")
  }
  
  ### Get the data sample and number of components
  n_components <- length(means) 
  n <- length(data)
  
  ### Calculating the belongings 
  gamma_num <- matrix(nrow = n, ncol = n_components)
  for(i in 1:n_components){
    gamma_num[,i] <- probs[i] * dnorm(data, mean = means[i], sd = stds[i])
  }
  gamma_den <- rowSums(gamma_num)
  
  # Handle cases where gamma_den might be zero by replacing zero values with a small number
  gamma_den[gamma_den == 0] <- .Machine$double.eps
  
  # Calculate gamma, handling division carefully to avoid NA values
  gamma <- gamma_num / gamma_den
  
  ### Check if data belongings sum to 1 (allowing for rounding errors).
  if (any(is.na(gamma)) || any(abs(rowSums(gamma) - 1) > .Machine$double.eps * n)){
    stop("Calculated belongings for each data point do not sum to 1.")
  }
  
  ### Check if any data belonging is negative
  if (any(gamma < 0)){
    stop("Calculated belongings contain negative values.")
  }
  
  return(gamma)
}

calculate_node_local_sufficient_statistics <- function(data, params){
    ### Calculating the data belongings
    belongings <- calculate_data_belongings(params = params, data = data)
    
    ### Calculating the local sufficient statistics for the node
    #### Number 1: sum of data belongings
    sufficient_stat_one <- colSums(belongings)
    
    #### Number 2: sum of product (data and data belongings)
    sufficient_stat_two <- colSums(belongings * matrix(data, nrow = length(data), ncol = ncol(belongings), byrow = FALSE))
    
    #### Number 3: sum of product (data squared and data belongings)
    sufficient_stat_three <- colSums(belongings * matrix(data^2, nrow = length(data), ncol = ncol(belongings), byrow = FALSE))
    
    ## Check if sufficient statistic is valid
    if(any(sufficient_stat_three < 0)){
        stop("Third sufficient statistic is negative")
    }
    return(list("suff_stat_one" = sufficient_stat_one, 
                "suff_stat_two" = sufficient_stat_two,
                "suff_stat_three" = sufficient_stat_three)
    )
}

estimate_parameters_with_sufficient_statistics <- function(global_sufficient_statistics){
    
    # Initialize vectors for means, standard deviations, and probabilities
    means <- numeric(length(global_sufficient_statistics$global_suff_stat_one))
    stds <- numeric(length(global_sufficient_statistics$global_suff_stat_one))
    probs <- numeric(length(global_sufficient_statistics$global_suff_stat_one))
    
    for (k in 1:length(global_sufficient_statistics$global_suff_stat_one)) {
        mean <- global_sufficient_statistics$global_suff_stat_two[k] / global_sufficient_statistics$global_suff_stat_one[k]
        var <- (global_sufficient_statistics$global_suff_stat_three[k] / global_sufficient_statistics$global_suff_stat_one[k]) - (mean^2)
        std <- sqrt(var)
        prob <- global_sufficient_statistics$global_suff_stat_one[k] / sum(global_sufficient_statistics$global_suff_stat_one)
        
        # Assign calculated values to the respective vectors
        means[k] <- mean
        stds[k] <- std
        probs[k] <- prob
    }
    
    return(list("mixing_probs" = probs,
                "means" = means,
                "stds" = stds))
}


calculate_log_likelihood_with_full_data <- function(data, params){
    # Initialize vectors for means, standard deviations, and probabilities
    mu <- params$means
    sigma <- params$stds
    weight <- params$mixing_probs
    k <- length(mu)
    n <- length(data)
    likelihood <- numeric(n)
    for (i in 1:n) {
        point_likelihood <- 0
        for (j in 1:k) {
            component_density <- weight[j] * dnorm(data[i], mean = mu[j], sd = sigma[j])
            point_likelihood <- point_likelihood + component_density
        }
        likelihood[i] <- point_likelihood
    }
    return(sum(log(likelihood)))
}

########### Multivariate code ###########
calculate_data_belongings_multivariate <- function(params, data) {
  
  # Extracting parameters
  means <- params$means
  covs <- params$covs
  probs <- params$mixing_probs
  
  # Ensure the number of means, covariances, and probabilities are the same
  if (length(means) != length(covs) || length(means) != length(probs)) {
    stop("Number of means, covariance matrices, and component probabilities must be equal.")
  }
  
  # Number of components and data points
  n_components <- length(means)
  n <- nrow(data)
  
  # Initialize a matrix to store the numerator for gamma
  gamma_num <- matrix(nrow = n, ncol = n_components)
  
  # Calculate the numerator for each component
  for (i in 1:n_components) {
    gamma_num[, i] <- probs[i] * mvtnorm::dmvnorm(data, mean = means[[i]], 
                                                  sigma = covs[[i]])  
  }
  
  # Calculate the denominator (sum across all components for each data point)
  gamma_den <- rowSums(gamma_num)
  
  # Calculate gamma (posterior probabilities for each component)
  gamma <- gamma_num / gamma_den
  
  # Check if gamma rows sum to 1 (allowing small rounding errors)
  if (any(abs(rowSums(gamma) - 1) > .Machine$double.eps * n)) {
    stop("Calculated belongings for each data point do not sum to 1.")
  }
  
  # Check for negative values in gamma
  if (any(gamma < 0)) {
    stop("Calculated belongings contain negative values.")
  }
  
  return(gamma)
}

calculate_node_local_sufficient_statistics_multivariate <- function(data, 
                                                                    params){
  ### Initiate needed values
  sufficient_stat_two <- list()
  sufficient_stat_three <- list()
  means <- params$means
  
  ### Calculating the data belongings
  belongings <- calculate_data_belongings_multivariate(params = params, 
                                                       data = data)
  k <- ncol(belongings) 
  
  ### Calculating the local sufficient statistics for the node
  #### Number 1: sum of data belongings
  sufficient_stat_one <- colSums(belongings)
  
  for(i in 1:k){
    #### Number 2: sum of product (data and data belongings)
    sufficient_stat_two[[i]] <- colSums(belongings[, i] * data)
    #### Number 3: sum of product (data squared and data belongings)
    centered_data <- t(apply(data, 1, function(row) row - means[[i]]))
    sufficient_stat_three[[i]] <-(t(centered_data) %*% (belongings[,i] * centered_data))
  }
  
  return(list("suff_stat_one" = sufficient_stat_one, 
              "suff_stat_two" = sufficient_stat_two,
              "suff_stat_three" = sufficient_stat_three)
  )
}

estimate_parameters_with_sufficient_statistics_multivariate <- function(global_sufficient_statistics){
  
  # Initialize vectors for means, standard deviations, and probabilities
  means <- list()
  covs <- list()
  probs <- c()
  k <- length(global_sufficient_statistics$global_suff_stat_one)
  
  for (i in 1:k) {
    means[[i]] <- global_sufficient_statistics$global_suff_stat_two[[i]] / global_sufficient_statistics$global_suff_stat_one[i]
    covs[[i]] <- global_sufficient_statistics$global_suff_stat_three[[i]] / global_sufficient_statistics$global_suff_stat_one[i]
    probs <- cbind(probs,
                   global_sufficient_statistics$global_suff_stat_one[i] / sum(global_sufficient_statistics$global_suff_stat_one))
  }
  
  return(list("mixing_probs" = probs,
              "means" = means,
              "covs" = covs))
}

calculate_log_likelihood_with_full_data_multivariate <- function(data, params) {
  
  # Initialize parameters for means, covariances, and mixing probabilities
  mu <- params$means
  sigma <- params$covs
  weight <- params$mixing_probs
  k <- length(mu)  # Number of components in the mixture
  n <- nrow(data)  # Number of data points
  likelihood <- numeric(n)  # To store the likelihood for each data point
  
  for (i in 1:n) {
    point_likelihood <- 0
    for (j in 1:k) {
      # Compute the component density using dmvnorm for multivariate normal
      component_density <- weight[j] * dmvnorm(data[i, ], mean = mu[[j]], sigma = sigma[[j]])
      point_likelihood <- point_likelihood + component_density
    }
    likelihood[i] <- point_likelihood
  }
  
  # Return the total log-likelihood by summing over the log of individual likelihoods
  return(sum(log(likelihood)))
}

