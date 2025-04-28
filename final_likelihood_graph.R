#####################File info#####################
# Code that stores logic for likelihood graph in dissertation (univariate simulation).

# The functions created is:
#     1. GMM_using_IEM
#     2. GMM_using_parallel
#     3. centralised_model
#     4. run_and_measure_IEM
#     5. run_and_measure_parallel
#     6. run_and_measure_centralised

# Compare models 
file_path <- '>>files location on your computer<<'
source(paste(file_path, 'final_IEM_model_univariate.R', sep = ""))
source(paste(file_path, 'final_parallel_model_univariate.R', sep = ""))
source(paste(file_path, 'final_centralised_model_univariate.R', sep = ""))
library(ggplot2)
library(dplyr)
set.seed(123)


# parameters
# SET MODEL HYPERPARAMETERS AND SETTINGS
K <- 3 # number of components
n <- 1000# sample size
m <- 4 # number of nodes
n_m <- n/m # sample size per node

# TRUE PARAMETERS
pi_1 <- 0.3
pi_2 <- 0.2
pi_3 <- 0.5
mu_1 <- 0
mu_2 <- 5
mu_3 <- 2
si_1 <- 1
si_2 <- 0.3
si_3 <- 3
probs <- c(pi_1, pi_2, pi_3)
means <- c(mu_1, mu_2, mu_3)
stds <- c(si_1, si_2, si_3)
params <- list("mixing_probs" = probs, "means" = means, "stds" = stds)

# INITIAL ESTIMATES
pi_1_init <- 0.4
pi_2_init <- 0.2
pi_3_init <- 0.4
mu_1_init <- 0.2
mu_2_init <- 4
mu_3_init <- 1
si_1_init <- 1
si_2_init <- 0.5
si_3_init <- 2
probs_init <- c(pi_1_init, pi_2_init, pi_3_init)
means_init <- c(mu_1_init, mu_2_init, mu_3_init)
stds_init <- c(si_1_init, si_2_init, si_3_init)
params_init <- list("mixing_probs" = probs_init,
                    "means" = means_init,
                    "stds" = stds_init)

# DATA GENERATION
data_chunks <- lapply(1:m, function(x) {
  generate_univariate_mixture_gaussian_data(means = means,
                                            stds = stds,
                                            probs = probs,
                                            n = n_m)
})

GMM_using_IEM <- function(data,
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
  column_names <- c("iteration", "likelihood")
  
  iteration <- 0
  
  # Initialize results matrix with the correct number of columns
  results <- matrix(ncol = length(column_names), nrow = 0)
  colnames(results) <- column_names
  params_diff <- 10000
  while (params_diff > tol){
    iteration <- iteration + 1
    for (i in 1:n_nodes){
      iteration_inner <- i/n_nodes
      local_statistics_old <- local_statistics_set[[i]]
      local_statistics_new <- calculate_node_local_sufficient_statistics(
        data = data[[i]],
        params = params_old)
      global_statistics_new <- update_univariate_global_statistics(
        global_statistics_old = global_statistics_old,
        local_statistics_old = local_statistics_old,
        local_statistics_new = local_statistics_new)
      params_new <- estimate_parameters_with_sufficient_statistics(
        global_sufficient_statistics = global_statistics_new)
      # Combine new results into a row
      likelihood <- calculate_log_likelihood_with_full_data(unlist(data),
                                                            params_new)
      new_row <- c(iteration + iteration_inner, likelihood)
      
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
  return(list('results' = results, 'params' = params_new))
}

GMM_using_parallel <- function(data, params_init,
                               tol = 0.05,
                               max_iter = 1000,
                               worker_nodes = length(data)){
  params_old <- params_init
  
  # Determine the number of components (k)
  k <- length(params_old$means)
  
  # Create column names dynamically
  column_names <- c("iteration", "likelihood")
  
  iteration <- 0.00
  
  # Initialize results matrix with the correct number of columns
  results <- matrix(ncol = length(column_names), nrow = 0)
  colnames(results) <- column_names
  
  # Set up cluster
  cl <- makeCluster(worker_nodes, type = "PSOCK")
  clusterExport(cl, list("calculate_node_local_sufficient_statistics", 
                         "calculate_data_belongings"
  ))
  param_diff <- 1000
  while(param_diff > tol){
    iteration <- iteration + 1
    if (iteration > max_iter){
      break
    }
    
    # Export current parameters to the workers
    clusterExport(cl, varlist = "params_old", envir = environment())
    
    # E-Step (done on workers)
    local_sufficient_statistics <- parLapply(cl, data, function(d){
      calculate_node_local_sufficient_statistics(data = d, 
                                                 params = params_old)
    })
    
    # Get global sufficient statistics (manager node)
    global_sufficient_statistics <- do.call(calculate_global_statistics, 
                                            local_sufficient_statistics)
    
    # Update the parameters
    params_new <- estimate_parameters_with_sufficient_statistics(
      global_sufficient_statistics = global_sufficient_statistics)
    log_likelihood_new <- calculate_log_likelihood_with_full_data(
        data = unlist(data), params = params_new)
    
    param_diff <- max(abs(unlist(params_new) - unlist(params_old)))
    new_row <- c(iteration, log_likelihood_new)
    
    # Append the new row to results
    results <- rbind(results, new_row)
    params_old <- params_new
  }
  stopCluster(cl)
  
  return(list("params" = params_old,
              "performance" = results))
}

centralised_model <- function(data,
                              initial_params,
                              tol = 0.05,
                              max_iter = 1000){
  params_old <- initial_params
  # Determine the number of components (k)
  k <- length(params_old$means)
  
  # Create column names dynamically
  column_names <- c("iteration", "likelihood")
  
  iteration <- 0.00
  
  # Initialize results matrix with the correct number of columns
  results <- matrix(ncol = length(column_names), nrow = 0)
  colnames(results) <- column_names
  params_diff <- 10000
  while(params_diff > tol){
    iteration <- iteration + 1
    if (iteration > max_iter){
      break
    }
    gamma <- calculate_data_belongings(params_old, data)
    params_new <- estimate_parameters_with_full_data(data,
                                                     gamma)
    
    params_diff <- max(abs(unlist(params_new) - unlist(params_old)))
    likelihood <- calculate_log_likelihood_with_full_data(data, params_new)
    new_row <- c(iteration, likelihood)
    
    # Append the new row to results
    results <- rbind(results, new_row)
    params_old <- params_new
  }
  return(list("params" = params_old,
              "performance" = results))
}



run_and_measure_IEM <- function(data_chunks, params_init){
  # Get initial global and local statistics
  initial_local_statistics <- lapply(1:m, function(x) calculate_node_local_sufficient_statistics(
    data_chunks[[x]], params = params_init
  ))
  initial_global_statistics <- list('global_suff_stat_one' = 0,
                                    'global_suff_stat_two' = 0,
                                    'global_suff_stat_three' = 0)
  
  for(local_sufficient_stat in initial_local_statistics){
    initial_global_statistics$global_suff_stat_one <- initial_global_statistics$global_suff_stat_one + local_sufficient_stat$suff_stat_one
    initial_global_statistics$global_suff_stat_two <- initial_global_statistics$global_suff_stat_two + local_sufficient_stat$suff_stat_two
    initial_global_statistics$global_suff_stat_three <- initial_global_statistics$global_suff_stat_three + local_sufficient_stat$suff_stat_three
  }
  start_time <- Sys.time()
  iem_distributed_model_results <- GMM_using_IEM(data_chunks,
                                                 initial_local_statistics,
                                                 initial_global_statistics,
                                                 params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = iem_distributed_model_results,
              'time_taken' = total_time))
}
incremental_test_results <- run_and_measure_IEM(data_chunks,
                                                params_init)

run_and_measure_parallel <- function(data_chunks, 
                                     params_init){
  start_time <- Sys.time()
  parallel_distributed_model_results <- GMM_using_parallel(data_chunks,
                                                           params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = parallel_distributed_model_results,
              'time_taken' = total_time))
}
parallel_test_results <- run_and_measure_parallel(data_chunks,
                                                  params_init)
run_and_measure_centralised <- function(data_chunks,
                                        params_init){
  start_time <- Sys.time()
  data <- NULL
  for (i in 1:length(data_chunks)){
    data <- rbind(data, data_chunks[[i]])
  }
  centralised_model_results <- centralised_model(data,
                                                 params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = centralised_model_results,
              'time_taken' = total_time))
}
centralised_test_results <- run_and_measure_centralised(data_chunks,
                                                        params_init)


parallel_results <- as.data.frame(parallel_test_results$model_results$performance)
incremental_results <- as.data.frame(incremental_test_results$model_results$results)
centralised_results <- as.data.frame(centralised_test_results$model_results$performance)
parallel_results$type <- "Parallel"
incremental_results$type <- "Incremental"
centralised_results$type <- "Centralised"

combined_data <- bind_rows(parallel_results, incremental_results, centralised_results)
ggplot(combined_data, aes(x = iteration, y = likelihood, color = type, group = type)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  labs(title = "Log-Likelihoods Across Iterations",
       x = "Iterations",
       y = "Log-Likelihood") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        legend.title = element_blank())

