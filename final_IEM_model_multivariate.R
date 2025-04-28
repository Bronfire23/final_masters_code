#####################File info#####################
# Code that stores logic for Incremental expectation maximisation model (multivariate).

# The functions created is:
#     1. update_multivariate_global_statistics
#     2. MGMM_using_IEM
file_path <- '>>files location on your computer<<'
source(paste(file_path, 'final_all_models_use.R', sep = ""))

update_multivariate_global_statistics <- function(global_statistics_old,
                                                local_statistics_old, 
                                                local_statistics_new){
  # Update global stats
  global_suff_stat_one <- global_statistics_old[[1]] + local_statistics_new[[1]] - local_statistics_old[[1]]
  global_suff_stat_two <- Map(function(global_old, local_old, local_new){
    global_old - local_old + local_new
  }, global_statistics_old[[2]], local_statistics_old[[2]], local_statistics_new[[2]])
  global_suff_stat_three <- Map(function(global_old, local_old, local_new){
    global_old - local_old + local_new
  }, global_statistics_old[[3]], local_statistics_old[[3]], local_statistics_new[[3]])
  
  # Check validity
  if(sum(global_suff_stat_one) <= 0){
    stop("The sum of the first global sufficient statistics must be positive")
  }
  
  return(list("global_suff_stat_one" = global_suff_stat_one,
              "global_suff_stat_two" = global_suff_stat_two,
              "global_suff_stat_three" = global_suff_stat_three))
}

MGMM_using_IEM <- function(data,
                                  initial_local_statistics,
                                  initial_global_statistics,
                                  initial_params,
                                  max_iter = 100,
                                  tol = 0.05){
  n_nodes <- length(data)
  params_old <- initial_params
  global_statistics_old <- initial_global_statistics
  local_statistics_set <- initial_local_statistics
  
  # Determine the number of components (k)
  k <- length(params_old$means)
  
  # Create column names dynamically
  column_names <- c("iteration_outer", "iteration_inner")
  iteration_outer <- 0
  iteration_inner <- 0
  
  # Initialize results matrix with the correct number of columns
  results <- matrix(ncol = length(column_names), nrow = 0)
  colnames(results) <- column_names
  params_diff <- 10000
  
  # Update parameters
  while (params_diff > tol){
    iteration_outer <- iteration_outer + 1
    iteration_inner <- 0
    for (i in 1:n_nodes){
      iteration_inner <- iteration_inner + 1
      local_statistics_old <- local_statistics_set[[i]]
      local_statistics_new <- calculate_node_local_sufficient_statistics_multivariate(
        data = data[[i]],
        params = params_old)
      global_statistics_new <- update_multivariate_global_statistics(
        global_statistics_old = global_statistics_old,
        local_statistics_old = local_statistics_old,
        local_statistics_new = local_statistics_new)
      params_new <- estimate_parameters_with_sufficient_statistics_multivariate(
        global_sufficient_statistics = global_statistics_new)
      # Combine new results into a row
      new_row <- c(iteration_outer, iteration_inner)
      
      # Append the new row to results
      results <- rbind(results, new_row)
      params_diff <- max(abs(unlist(params_new) - unlist(params_old)))
      if (params_diff < 0.05){ #break to check at each update
        break
      }
      local_statistics_set[[i]] <- local_statistics_new
      global_statistics_old <- global_statistics_new
      params_old <- params_new
    }
  }
  return(list('performance' = results, 'params' = params_new))
}


