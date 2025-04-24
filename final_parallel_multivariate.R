#####################File info#####################
# Code that stores logic for Parallel expectation maximisation model (multivariate).

# The functions created is:
#     1. calculate_global_statistics_multivariate
#     2. MGMM_using_parallel

# libraries
file_path <- '/Users/bronwynmccall/Documents/GitHub/final_masters_code/'
source(paste(file_path, 'final_all_models_use.R', sep = ""))
library(parallel)

## Function to calculate global statistics for parallel approach
calculate_global_statistics_multivariate <- function(node_local_statistics) {
  k <- length(node_local_statistics[[1]]$suff_stat_two)  # Number of components
  
  # Summing the first sufficient statistic across nodes
  global_suff_stat_one <- Reduce(`+`, lapply(node_local_statistics, `[[`, "suff_stat_one"))
  
  # Summing the second and third sufficient statistics across nodes
  global_suff_stat_two <- lapply(1:k, function(i) {
    Reduce(`+`, lapply(node_local_statistics, function(node) node$suff_stat_two[[i]]))
  })
  
  global_suff_stat_three <- lapply(1:k, function(i) {
    Reduce(`+`, lapply(node_local_statistics, function(node) node$suff_stat_three[[i]]))
  })
  
  return(list(
    "global_suff_stat_one" = global_suff_stat_one,
    "global_suff_stat_two" = global_suff_stat_two,
    "global_suff_stat_three" = global_suff_stat_three
  ))
}

MGMM_using_parallel <- function(data, params_init,
                               tol = 0.05,
                               max_iter = 1000,
                               worker_nodes = length(data)){
  params_old <- params_init
  
  # Determine the number of components (k)
  k <- length(params_old$means)
  
  iteration <- 0
  
  # Set up cluster
  cl <- makeCluster(worker_nodes, type = "PSOCK")
  clusterExport(cl, list("calculate_node_local_sufficient_statistics_multivariate", 
                         "calculate_data_belongings_multivariate",
                         "dmvnorm"
  ))
  param_diff <- 1000
  # Update parameters
  while(param_diff > tol){
    iteration <- iteration + 1
    if (iteration > max_iter){
      break
    }
    
    # Export current parameters to the workers
    clusterExport(cl, varlist = "params_old", envir = environment())
    
    # E-Step (done on workers)
    local_sufficient_statistics <- parLapply(cl, data, function(d){
      calculate_node_local_sufficient_statistics_multivariate(data = d, 
                                                 params = params_old)
    })
    
    # Get global sufficient statistics (manager node)
    global_sufficient_statistics <- calculate_global_statistics_multivariate(
      local_sufficient_statistics)
    
    
    # Update the parameters
    params_new <- estimate_parameters_with_sufficient_statistics_multivariate(
      global_sufficient_statistics = global_sufficient_statistics)

    param_diff <- max(abs(unlist(params_new) - unlist(params_old)))
    params_old <- params_new
  }
  stopCluster(cl)
  
  return(list("params" = params_old,
              "performance" = iteration))
}


