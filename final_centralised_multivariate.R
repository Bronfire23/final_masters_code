#####################File info#####################
# Code that stores logic for Centralised expectation maximisation model (multivariate).

# The functions created is:
#     1. estimate_parameters_with_full_data_multivariate
#     2. MGMM_using_centralised

file_path <- '>>files location on your computer<<'
source(paste(file_path, 'final_data_generation.R', sep = ""))
source(paste(file_path, 'final_all_models_use.R', sep = ""))
library(mvtnorm)


estimate_parameters_with_full_data_multivariate <- function(data, data_belongings){
  n_k <- colSums(data_belongings)
  n <- nrow(data)
  k <- ncol(data_belongings)
  
  probs <- numeric(k)
  means <- list()
  covs <- list()
  
  for(i in 1:k){
    # Mixing probabilities
    probs[i] <- n_k[i]/n
    
    # Means
    means[[i]] <- colSums(data_belongings[, i] * data) / n_k[i]
    
    # Covariance matrices
    centered_data <- t(apply(data, 1, function(row) row - means[[i]]))
    covs[[i]] <- (t(centered_data) %*% (data_belongings[,i] * centered_data))/n_k[i]
  }
  
  return(list('mixing_probs' = probs,
              'means' = means,
              'covs' = covs))
}

MGMM_using_centralised <- function(data, params_init,
                                tol = 0.05,
                                max_iter = 1000,
                                worker_nodes = length(data)){
  params_old <- params_init
  
  # Determine the number of components (k)
  k <- length(params_old$means)
  
  iteration <- 0

  param_diff <- 1000
  while(param_diff > tol){
    iteration <- iteration + 1
    if (iteration > max_iter){
      break
    }
    gamma <- calculate_data_belongings_multivariate(params_old, data)
    params_new <- estimate_parameters_with_full_data_multivariate(data,
                                                                  gamma)
    
    param_diff <- max(abs(unlist(params_new) - unlist(params_old)))
    params_old <- params_new
  }
  return(list("params" = params_old,
              "performance" = iteration))
}



