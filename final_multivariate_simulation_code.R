#####################File info#####################
# Code that stores logic for multivariate simulation.

# The functions created is:
#     1. run_and_measure_IEM_multivariate 
#     2. run_and_measure_parallel_multivariate
#     3. run_and_measure_centralised_multivariate

# libraries
file_path <- '>>files location on your computer<<'
source(paste(file_path, 'final_data_generation.R', sep = ""))
source(paste(file_path, 'final_all_models_use.R', sep = ""))
source(paste(file_path, 'final_IEM_model_multivariate.R', sep = ""))
source(paste(file_path, 'final_parallel_multivariate.R', sep = ""))
source(paste(file_path, 'final_centralised_multivariate.R', sep = ""))
library(psych)
library(openxlsx)
library(patchwork)
library(ggplot2)

## parameters
n <- 100 # data points 100, 1000, 10000 -> set variable
m <- 4 # number of nodes
n_m <- n/m
p <- 2 # number of components
## parameters
set.seed(123) # For reproducibility

# True parameters
mu1 <- c(-2, 2)
sigma1 <- matrix(c(1.2, 0.5, 0.5, 0.8), nrow = 2)
mu2 <- c(1, -1)
sigma2 <- matrix(c(1.5, 0.9, 0.9, 2.0), nrow = 2)
probs <- c(0.4, 0.6)
means <- list(mu1, mu2)
covs <- list(sigma1, sigma2)
params <- list('means' = means,
               'covs' = covs,
               'mixing_probs' = probs)
params_init <- list()
params_init$mixing_probs <- c(0.4018663, 0.5981337)
params_init$means[[1]] <- c(-3.120951, 1.539645)
params_init$means[[2]] <- c(4.1174166, -0.8589832)
params_init$covs[[1]] <- matrix(c(1.2129288, 0.6087991,
                                  0.6087991, 0.6734939), ncol = 2, byrow = T)
params_init$covs[[2]] <- matrix(c(1.431315, 0.938921,
                                  0.938921, 2.035981), ncol = 2, byrow = T)

# # Generating less biased initial values - uncomment out if you want to use
# generate_noisy_params <- function(true_means, true_covs, true_probs, noise_factor_mean, noise_factor_cov, noise_factor_probs) {
#   means_init <- lapply(true_means, function(mu) {
#     mu + rnorm(length(mu), mean = 0, sd = noise_factor_mean)
#   })
# 
#   covs_init <- lapply(true_covs, function(sigma) {
#     sigma + matrix(rnorm(length(sigma), mean = 0, sd = noise_factor_cov), nrow = nrow(sigma), ncol = ncol(sigma))
#   })
# 
#   probs_init <- true_probs + runif(length(true_probs),
#                                    min = -noise_factor_probs,
#                                    max = noise_factor_probs)
#   probs_init <- probs_init / sum(probs_init)
# 
#   list("mixing_probs" = probs_init,
#        "means" = means_init,
#        "covs" = covs_init)
# }
# 
# # Adjust noise levels as needed
# noise_factor_mean <- 2
# noise_factor_cov <- 0.1
# noise_factor_probs <- 0.1
# 
# params_init <- generate_noisy_params(means, covs, probs, noise_factor_mean,
#                                      noise_factor_cov, noise_factor_probs)
# 
# # Function to make a matrix symmetric
# make_symmetric <- function(mat) {
#   (mat + t(mat)) / 2
# }
# 
# # Apply to initial covariance matrices
# params_init$covs[[1]] <- make_symmetric(params_init$covs[[1]])
# params_init$covs[[2]] <- make_symmetric(params_init$covs[[2]])
# print(params_init)

## Generate data for nodes
data_chunks <- lapply(1:m, 
                      function(x) generate_multivariate_mixture_gaussian_data(
                        means, 
                        covs,
                        probs,
                        n/m))

data_full <- NULL
for (i in 1:length(data_chunks)){
  data_full <- rbind(data_full, data_chunks[[i]])
} 

data_df <- as.data.frame(data_full)
colnames(data_df) <- c('x', 'y')

# Plot data
p_main <- ggplot(data_df, aes(x, y)) +
  geom_point(alpha = 0.3) +
  stat_density_2d(aes(color = after_stat(level)), geom = "contour") +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(x = "Var 1", y = "Var 2") +
  theme(legend.position = "none")

# Create density plot for x with corrected legends
p_x <- ggplot(data_df, aes(x)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  stat_function(aes(color = "Component 1"), 
                fun = function(x) dnorm(x, mean = mu1[1], 
                                        sd = sqrt(sigma1[1,1]))) +
  stat_function(aes(color = "Component 2"), 
                fun = function(x) dnorm(x, mean = mu2[1], 
                                        sd = sqrt(sigma2[1,1]))) +
  theme_void() +
  scale_color_manual(name = "Legend", 
                     values = c("Component 1" = "red", 
                                "Component 2" = "green"))

# Create density plot for y with corrected legends
p_y <- ggplot(data_df, aes(y)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  stat_function(color = "red", 
                fun = function(x) dnorm(x, mean = mu1[2], 
                                        sd = sqrt(sigma1[2,2]))) +
  stat_function(color = "green", 
                fun = function(x) dnorm(x, mean = mu2[2], 
                                        sd = sqrt(sigma2[2,2]))) +
  theme_void() +
  coord_flip() +
  scale_x_reverse() +  
  scale_y_reverse() 


layout <- "
AB
AC
"

combined_plot <- p_y + p_main + p_x + plot_layout(design = layout, 
                                                  widths = c(1, 4), 
                                                  heights = c(4, 1))

combined_plot

run_and_measure_IEM_multivariate <- function(data_chunks, params_init){
  # Get initial global and local statistics
  p <- length(params_init$means)
  initial_local_statistics <- lapply(1:m, function(x) 
    calculate_node_local_sufficient_statistics_multivariate(data_chunks[[x]], 
                                                            params = params_init
    ))
  initial_global_statistics <- list('global_suff_stat_one' = 0,
                                    'global_suff_stat_two' = lapply(
                                      1:p, function(x) 0),
                                    'global_suff_stat_three' = lapply(
                                      1:p, function(x) 0))
  
  for(local_sufficient_stat in initial_local_statistics){
    initial_global_statistics$global_suff_stat_one <- 
      initial_global_statistics$global_suff_stat_one + 
      local_sufficient_stat$suff_stat_one
    initial_global_statistics$global_suff_stat_two <- 
      Map(`+`, 
          initial_global_statistics$global_suff_stat_two,
          local_sufficient_stat$suff_stat_two)
    initial_global_statistics$global_suff_stat_three <- 
      Map(`+`,
          initial_global_statistics$global_suff_stat_three,
          local_sufficient_stat$suff_stat_three)
  }
  start_time <- Sys.time()
  iem_distributed_model_results <- MGMM_using_IEM(
    data_chunks,
    initial_local_statistics,
    initial_global_statistics,
    params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = iem_distributed_model_results,
              'time_taken' = total_time))
}
incremental_test_results <- run_and_measure_IEM_multivariate(data_chunks, 
                                                             params_init)

run_and_measure_parallel_multivariate <- function(data_chunks, 
                                                  params_init){
  start_time <- Sys.time()
  parallel_distributed_model_results <- MGMM_using_parallel(data_chunks,
                                                            params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = parallel_distributed_model_results,
              'time_taken' = total_time))
}
parallel_test_results <- run_and_measure_parallel_multivariate(data_chunks,
                                                               params_init)

run_and_measure_centralised_multivariate <- function(data_chunks,
                                                     params_init){
  start_time <- Sys.time()
  data <- NULL
  for (i in 1:length(data_chunks)){
    data <- rbind(data, data_chunks[[i]])
  }
  centralised_model_results <- MGMM_using_centralised(data,
                                                      params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = centralised_model_results,
              'time_taken' = total_time))
}
centralised_test_results <- run_and_measure_centralised_multivariate(
  data_chunks,
  params_init)

# Initialize storage for the results of each model
incremental_means_1 <- matrix(nrow = 1000, ncol = p)
incremental_means_2 <- matrix(nrow = 1000, ncol = p)
incremental_cov <- matrix(nrow = 1000, ncol = p)
incremental_var_1 <- matrix(nrow = 1000, ncol = p)
incremental_var_2 <- matrix(nrow = 1000, ncol = p)
incremental_mixing_probs <- matrix(nrow = 1000, ncol = p)
incremental_time <- c()

parallel_means_1 <- matrix(nrow = 1000, ncol = p)
parallel_means_2 <- matrix(nrow = 1000, ncol = p)
parallel_cov <- matrix(nrow = 1000, ncol = p)
parallel_var_1 <- matrix(nrow = 1000, ncol = p)
parallel_var_2 <- matrix(nrow = 1000, ncol = p)
parallel_mixing_probs <- matrix(nrow = 1000, ncol = p)
parallel_time <- c()

centralised_means_1 <- matrix(nrow = 1000, ncol = p)
centralised_means_2 <- matrix(nrow = 1000, ncol = p)
centralised_cov <- matrix(nrow = 1000, ncol = p)
centralised_var_1 <- matrix(nrow = 1000, ncol = p)
centralised_var_2 <- matrix(nrow = 1000, ncol = p)
centralised_mixing_probs <- matrix(nrow = 1000, ncol = p)
centralised_time <- c()

# Number of iterations
iterations <- 1000

# Main loop for 1000 iterations
for (i in 1:iterations) {
  # Generate data chunks for each iteration
  data_chunks <- lapply(1:m, function(x) {
    generate_multivariate_mixture_gaussian_data(means = params$means,
                                                covariances = params$covs,
                                                probs = params$mixing_probs,
                                                n = n_m)
  })
  full_data <- NULL # for calculating likelihood
  for (j in 1:length(data_chunks)){
    full_data <- rbind(full_data, data_chunks[[j]])
  }
  # Incremental EM
  IEM_result <- run_and_measure_IEM_multivariate(data_chunks, params_init)
  incremental_mixing_probs[i,] <- IEM_result$model_results$params$mixing_probs
  for(it in 1:p){
    incremental_means_1[i, it] <- IEM_result$model_results$params$means[[it]][1]
    incremental_means_2[i, it] <- IEM_result$model_results$params$means[[it]][2]
    incremental_cov[i, it] <- IEM_result$model_results$params$covs[[it]][1,2]
    incremental_var_1[i, it] <- IEM_result$model_results$params$covs[[it]][1,1]
    incremental_var_2[i, it] <- IEM_result$model_results$params$covs[[it]][2,2]
  }
  incremental_time <- rbind(incremental_time, as.numeric(IEM_result$time_taken))
  
  # Parallel EM
  parallel_result <- run_and_measure_parallel_multivariate(data_chunks, 
                                                           params_init)
  parallel_mixing_probs[i,] <- parallel_result$model_results$params$mixing_probs
  for(it in 1:p){
    parallel_means_1[i, it] <- parallel_result$
      model_results$params$means[[it]][1]
    parallel_means_2[i, it] <- parallel_result$
      model_results$params$means[[it]][2]
    parallel_cov[i, it] <- parallel_result$model_results$
      params$covs[[it]][1,2]
    parallel_var_1[i, it] <- parallel_result$model_results$
      params$covs[[it]][1,1]
    parallel_var_2[i, it] <- parallel_result$model_results$
      params$covs[[it]][2,2]
  }
  parallel_time <- rbind(parallel_time, as.numeric(parallel_result$time_taken))
  
  # Centralised EM
  centralised_result <- run_and_measure_centralised_multivariate(
    data_chunks, params_init)
  centralised_mixing_probs[i,] <- centralised_result$model_results$
    params$mixing_probs
  for(it in 1:p){
    centralised_means_1[i, it] <- centralised_result$model_results$
      params$means[[it]][1]
    centralised_means_2[i, it] <- centralised_result$model_results$
      params$means[[it]][2]
    centralised_cov[i, it] <- centralised_result$model_results$
      params$covs[[it]][1,2]
    centralised_var_1[i, it] <- centralised_result$model_results$
      params$covs[[it]][1,1]
    centralised_var_2[i, it] <- centralised_result$model_results$
      params$covs[[it]][2,2]
  }
  centralised_time <- rbind(centralised_time, 
                            as.numeric(centralised_result$time_taken))
}

# # Set up a 3x1 grid for subplots (3 rows, 1 column for each component)
lapply(1:p, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns

  # Plot incremental means
  hist(incremental_means_1[, i],
       main = paste("Incremental (Component", i, ")"),
       xlab = "Mean 1 Estimates",
       col = "lightblue")

  # Plot parallel means
  hist(parallel_means_1[, i],
       main = paste("Parallel (Component", i, ")"),
       xlab = "Mean 1 Estimates",
       col = "lightgreen")

  # Plot centralised means
  hist(centralised_means_1[, i],
       main = paste("Centralised (Component", i, ")"),
       xlab = "Mean 1 Estimates",
       col = "lightpink")
})

lapply(1:p, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns

  # Plot incremental means
  hist(incremental_means_2[, i],
       main = paste("Incremental (Component", i, ")"),
       xlab = "Mean 2 Estimates",
       col = "lightblue")

  # Plot parallel means
  hist(parallel_means_2[, i],
       main = paste("Parallel (Component", i, ")"),
       xlab = "Mean 2 Estimates",
       col = "lightgreen")

  # Plot centralised means
  hist(centralised_means_2[, i],
       main = paste("Centralised (Component", i, ")"),
       xlab = "Mean 2 Estimates",
       col = "lightpink")
})

lapply(1:p, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns

  # Plot incremental means
  hist(incremental_cov[, i],
       main = paste("Incremental (Component", i, ")"),
       xlab = "Covariance Estimates",
       col = "lightblue")

  # Plot parallel means
  hist(parallel_cov[, i],
       main = paste("Parallel (Component", i, ")"),
       xlab = "Covariance Estimates",
       col = "lightgreen")

  # Plot centralised means
  hist(centralised_cov[, i],
       main = paste("Centralised (Component", i, ")"),
       xlab = "Covariance Estimates",
       col = "lightpink")
})

lapply(1:p, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns

  # Plot incremental means
  hist(incremental_var_1[, i],
       main = paste("Incremental (Component", i, ")"),
       xlab = "Variance 1 Estimates",
       col = "lightblue")

  # Plot parallel means
  hist(parallel_var_1[, i],
       main = paste("Parallel (Component", i, ")"),
       xlab = "Variance 1 Estimates",
       col = "lightgreen")

  # Plot centralised means
  hist(centralised_var_1[, i],
       main = paste("Centralised (Component", i, ")"),
       xlab = "Variance 1 Estimates",
       col = "lightpink")
})

lapply(1:p, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns

  # Plot incremental means
  hist(incremental_var_2[, i],
       main = paste("Incremental (Component", i, ")"),
       xlab = "Variance 2 Estimates",
       col = "lightblue")

  # Plot parallel means
  hist(parallel_var_2[, i],
       main = paste("Parallel (Component", i, ")"),
       xlab = "Variance 2 Estimates",
       col = "lightgreen")

  # Plot centralised means
  hist(centralised_var_2[, i],
       main = paste("Centralised (Component", i, ")"),
       xlab = "Variance 2 Estimates",
       col = "lightpink")
})

# Set up a 3x1 grid for subplots (3 rows, 1 column for each component)
lapply(1:p, function(i) {
  # Create a new plot window for each component
  par(mfrow = c(1, 3)) # 1 row, 3 columns

  # Plot incremental means
  hist(incremental_mixing_probs[, i],
       main = paste("Incremental (Component", i, ")"),
       xlab = "Mixing Probability Estimates",
       col = "lightblue")

  # Plot parallel means
  hist(parallel_mixing_probs[, i],
       main = paste("Parallel (Component", i, ")"),
       xlab = "Mixing Probability Estimates",
       col = "lightgreen")

  # Plot centralised means
  hist(centralised_mixing_probs[, i],
       main = paste("Centralised (Component", i, ")"),
       xlab = "Mixing Probability Estimates",
       col = "lightpink")
})

# Display time taken for each model
par(mfrow = c(1,3))
hist(incremental_time, main = "Incremental time", xlab ='Time', col = "lightblue")
hist(parallel_time, main = "Parallel time", xlab = 'Time', col = 'lightgreen')
hist(centralised_time, main = "Centralised time", xlab = 'Time', col = 'lightpink')


desc_iem_means_1 <- describe(incremental_means_1)
desc_iem_means_2 <- describe(incremental_means_2)
desc_parallel_means_1 <- describe(parallel_means_1)
desc_parallel_means_2 <- describe(parallel_means_2)
decs_centralised_means_1 <- describe(centralised_means_1)
decs_centralised_means_2 <- describe(centralised_means_2)

desc_iem_covs <- describe(incremental_cov)
desc_iem_var_1 <- describe(incremental_var_1)
desc_iem_var_2 <- describe(incremental_var_2)
desc_parallel_covs <- describe(parallel_cov)
desc_parallel_var_1 <- describe(parallel_var_1)
desc_parallel_var_2 <- describe(parallel_var_2)
desc_centralised_covs <- describe(centralised_cov)
desc_centralised_var_1 <- describe(centralised_var_1)
desc_centralised_var_2 <- describe(centralised_var_2)

desc_iem_mixing_probs <- describe(incremental_mixing_probs)
desc_parallel_mixing_probs <- describe(parallel_mixing_probs)
desc_centralised_mixing_probs <- describe(centralised_mixing_probs)

# Create a new workbook
wb <- createWorkbook()

# Add worksheets and write data
addWorksheet(wb, "Incremental Means1")
writeData(wb, sheet = "Incremental Means1", desc_iem_means_1)

addWorksheet(wb, "Incremental Means2")
writeData(wb, sheet = "Incremental Means2", desc_iem_means_2)

addWorksheet(wb, "Parallel Means1")
writeData(wb, sheet = "Parallel Means1", desc_parallel_means_1)

addWorksheet(wb, "Parallel Means2")
writeData(wb, sheet = "Parallel Means2", desc_parallel_means_2)

addWorksheet(wb, "Centralised Means1")
writeData(wb, sheet = "Centralised Means1", decs_centralised_means_1)

addWorksheet(wb, "Centralised Means2")
writeData(wb, sheet = "Centralised Means2", decs_centralised_means_2)

addWorksheet(wb, "Incremental Covs")
writeData(wb, sheet = "Incremental Covs", desc_iem_covs)


addWorksheet(wb, "Incremental Var1")
writeData(wb, sheet = "Incremental Var1", desc_iem_var_1)

addWorksheet(wb, "Incremental Var2")
writeData(wb, sheet = "Incremental Var2", desc_iem_var_2)

addWorksheet(wb, "Parallel Covs")
writeData(wb, sheet = "Parallel Covs", desc_parallel_covs)

addWorksheet(wb, "Parallel Var1")
writeData(wb, sheet = "Parallel Var1", desc_parallel_var_1)

addWorksheet(wb, "Parallel Var2")
writeData(wb, sheet = "Parallel Var2", desc_parallel_var_2)

addWorksheet(wb, "Centralised Covs")
writeData(wb, sheet = "Centralised Covs", desc_centralised_covs)

addWorksheet(wb, "Centralised Var1")
writeData(wb, sheet = "Centralised Var1", desc_centralised_var_1)

addWorksheet(wb, "Centralised Var2")
writeData(wb, sheet = "Centralised Var2", desc_centralised_var_2)

addWorksheet(wb, "Incremental Mixing Probs")
writeData(wb, sheet = "Incremental Mixing Probs", desc_iem_mixing_probs)

addWorksheet(wb, "Parallel Mixing Probs")
writeData(wb, sheet = "Parallel Mixing Probs", desc_parallel_mixing_probs)

addWorksheet(wb, "Centralised Mixing Probs")
writeData(wb, sheet = "Centralised Mixing Probs", desc_centralised_mixing_probs)

# Save the workbook
saveWorkbook(wb, file = "medium_cov_new_descriptive_statistics_MV_2_100.xlsx", 
             overwrite = TRUE)

cat("Results saved to descriptive_statistics.xlsx")

