#####################File info#####################
# This files contents are not used in paper but is nice to have to check model
# results of the multivariate model.

# libraries
file_path <- '/Users/bronwynmccall/Documents/GitHub/final_masters_code/'
source(paste(file_path, 'final_data_generation.R', sep = ""))
source(paste(file_path, 'final_all_models_use.R', sep = ""))
source(paste(file_path, 'final_IEM_model_multivariate.R', sep = ""))
source(paste(file_path, 'final_parallel_multivariate.R', sep = ""))
source(paste(file_path, 'final_centralised_multivariate.R', sep = ""))
library(psych)
library(openxlsx)
library(patchwork)
library(mclust)
library(ggplot2)


## parameters
n <- 1000 # data points 100, 1000, 10000 -> set variable
m <- 4 # number of nodes
n_m <- n/m
p <- 2 # number of components

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
                fun = function(x) probs[1]*dnorm(x, mean = mu1[1], 
                                        sd = sqrt(sigma1[1,1]))) +
  stat_function(aes(color = "Component 2"), 
                fun = function(x) probs[2]*dnorm(x, mean = mu2[1], 
                                        sd = sqrt(sigma2[1,1]))) +
  theme_void() +
  scale_color_manual(name = "Legend", 
                     values = c("Component 1" = "red", 
                                "Component 2" = "green"))

# Create density plot for y with corrected legends
p_y <- ggplot(data_df, aes(y)) +
  geom_density(fill = "lightblue", alpha = 0.5) +
  stat_function(color = "red", 
                fun = function(x) probs[1]*dnorm(x, mean = mu1[2], 
                                        sd = sqrt(sigma1[2,2]))) +
  stat_function(color = "green", 
                fun = function(x) probs[2]*dnorm(x, mean = mu2[2], 
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

#######################################################################
# Models
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



############################################################################
# Plot IEM results
calculate_ellipse <- function(mu, sigma, prob = 0.95) {
  # Get eigenvalues and eigenvectors of covariance matrix
  eigen_decomp <- eigen(sigma)
  # Calculate scaling factor for desired probability
  chi_sq <- qchisq(prob, df = 2)
  # Generate circle points
  theta <- seq(0, 2*pi, length.out = 100)
  circle <- cbind(cos(theta), sin(theta))
  # Transform circle to ellipse
  ellipse <- t(mu + t(circle %*% diag(sqrt(chi_sq * eigen_decomp$values)) %*% t(eigen_decomp$vectors)))
  return(as.data.frame(ellipse))
}

# Calculate ellipse points for each component
ellipse1 <- calculate_ellipse(
  mu = unlist(incremental_test_results$model_results$params$means[[1]]),
  sigma = incremental_test_results$model_results$params$covs[[1]]
)
ellipse2 <- calculate_ellipse(
  mu = unlist(incremental_test_results$model_results$params$means[[2]]),
  sigma = incremental_test_results$model_results$params$covs[[2]]
)

p_main <- ggplot(data_df, aes(x, y)) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  labs(x = "Var 1", y = "Var 2") +
  theme(legend.position = "none")

p_main + geom_point(data = data.frame(
  x = c(centralised_test_results$model_results$params$means[[1]][1],
        centralised_test_results$model_results$params$means[[2]][1]),
  y = c(centralised_test_results$model_results$params$means[[1]][2],
        centralised_test_results$model_results$params$means[[2]][2])), 
  size = 5, col = c('red', 'blue')) +   
  geom_path(data = ellipse1, aes(V1, V2), color = "red", linetype = "dashed") +
  geom_path(data = ellipse2, aes(V1, V2), color = "blue", linetype = "dashed")

p_main + geom_point(data = data.frame(
  x = c(parallel_test_results$model_results$params$means[[1]][1],
        parallel_test_results$model_results$params$means[[2]][1]),
  y = c(parallel_test_results$model_results$params$means[[1]][2],
        parallel_test_results$model_results$params$means[[2]][2])), 
  size = 5, col = c('green', 'orange')) +   
  geom_path(data = ellipse1, aes(V1, V2), color = "green", linetype = "dashed") +
  geom_path(data = ellipse2, aes(V1, V2), color = "orange", linetype = "dashed")

p_main + geom_point(data = data.frame(
  x = c(incremental_test_results$model_results$params$means[[1]][1],
        incremental_test_results$model_results$params$means[[2]][1]),
  y = c(incremental_test_results$model_results$params$means[[1]][2],
        incremental_test_results$model_results$params$means[[2]][2])), 
  size = 5, col = c('purple', 'yellow')) +   
  geom_path(data = ellipse1, aes(V1, V2), color = "purple", linetype = "dashed") +
  geom_path(data = ellipse2, aes(V1, V2), color = "yellow", linetype = "dashed")

