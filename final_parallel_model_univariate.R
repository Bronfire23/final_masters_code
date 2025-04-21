#####################File info#####################
# Code that stores logic for Parallel expectation maximisation model (univariate).

# The functions stored is:
#     1. calculate_global_statistics 
#     2. GMM_using_parallel

# libraries
file_path <- '/Users/bronwynmccall/Documents/GitHub/final_masters_code/'
source(paste(file_path, 'final_all_models_use.R', sep = ""))
source(paste(file_path,'final_data_generation.R', sep = ""))
library(parallel)

calculate_global_statistics <- function(...){
    local_statistics <- list(...)
    
    # Set initial global ss to nothing
    global_suff_stat_one <- NULL
    global_suff_stat_two <- NULL
    global_suff_stat_three <- NULL
    
    # Update global ss
    for(local_stat in local_statistics){
        if(is.null(global_suff_stat_one)){
            global_suff_stat_one <- local_stat$suff_stat_one
            global_suff_stat_two <- local_stat$suff_stat_two
            global_suff_stat_three <- local_stat$suff_stat_three
        } else{
            global_suff_stat_one <- global_suff_stat_one + local_stat$suff_stat_one
            global_suff_stat_two <- global_suff_stat_two + local_stat$suff_stat_two
            global_suff_stat_three <- global_suff_stat_three + local_stat$suff_stat_three
        }
    }
    
    # checks for validity
    if(sum(global_suff_stat_one) <= 0){
        stop("The sum of the first global sufficient statistics must be positi`ve")
    }
    
    if(any(global_suff_stat_three < 0)){
        stop("Third global sufficient statistics is negative")
    }
    return(list("global_suff_stat_one" = global_suff_stat_one,
                "global_suff_stat_two" = global_suff_stat_two,
                "global_suff_stat_three" = global_suff_stat_three))
}

GMM_using_parallel <- function(data, params_init,
                               tol = 0.05,
                               max_iter = 1000,
                               worker_nodes = length(data)){
    params_old <- params_init
    
    # Determine the number of components (k)
    k <- length(params_old$means)
    
    # Create column names dynamically
    column_names <- c("iteration")
    column_names <- c(column_names,
                      paste0("component_", 1:k, "_mean"),
                      paste0("component_", 1:k, "_stds"),
                      paste0("component_", 1:k, "_mixing_probs"))
    
    iteration <- 0
    
    # Initialize results matrix with the correct number of columns
    results <- matrix(ncol = length(column_names), nrow = 0)
    colnames(results) <- column_names
    
    # Set up cluster
    cl <- makeCluster(worker_nodes, type = "PSOCK")
    clusterExport(cl, list("calculate_node_local_sufficient_statistics", 
                           "calculate_data_belongings"
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
            calculate_node_local_sufficient_statistics(data = d, 
                                                       params = params_old)
        })
        
        # Get global sufficient statistics (manager node)
        global_sufficient_statistics <- do.call(calculate_global_statistics, 
                                                local_sufficient_statistics)
        
        # Update the parameters
        params_new <- estimate_parameters_with_sufficient_statistics(
            global_sufficient_statistics = global_sufficient_statistics)
        # log_likelihood_new <- calculate_log_likelihood_with_full_data(
        #     data = full_data, params = params_new)

        param_diff <- max(abs(unlist(params_new) - unlist(params_old)))
        new_row <- c(iteration, 
                     params_new$means, params_new$stds, params_new$mixing_probs)
        
        # Append the new row to results
        results <- rbind(results, new_row)
        params_old <- params_new
    }
    stopCluster(cl)
    
    return(list("params" = params_old,
                "performance" = results))
}

