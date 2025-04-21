#####################File info#####################
# Code that stores the logic for data generation (both univariate and multivariate).
# Code for both generating data and plotting data

# The functions stored is:
#     1. generate_univariate_mixture_gaussian_data
#     2. plot_univariate_gaussian_data
#     3. plot_gmm_predictions
#     4. generate_multivariate_mixture_gaussian_data
#     5. plot_mixture_histograms_with_density

#libraries
library(MASS)
library(reshape2)
library(ggplot2)

generate_univariate_mixture_gaussian_data <- function(means, stds, probs, n){
    ## Ensure different parameter lengths match
    if(length(means) != length(stds) || length(means) != length(probs)){
        stop("Number of means, standard deviations and components probabilities 
             must be equal.")
    }
    ## Generate data
    data <- matrix(nrow = n)
    for(i in 1:n){
        component <- which(rmultinom(1, 1, probs) == 1)
        data[i] <- rnorm(1, mean = means[component], sd = stds[component])
    }
    return(data)
}

plot_univariate_gaussian_data <- function(data, n_components, params, 
                                          overlay = TRUE, freq = FALSE, 
                                          main = "Histogram of data") {
  ## Extracting parameters
  means <- params$means
  stds <- params$stds
  probs <- params$mixing_probs
  
  ## Histogram
  hist(data, breaks = 50, col = "lightblue", xlab = "Value", ylim = c(0, 1.3),
       prob = TRUE, main = main)
  
  ## Overlayed line plot
  if (overlay) {
    colors <- rainbow(n_components)
    for (i in 1:n_components) {
      curve(dnorm(x, mean = means[i], sd = stds[i]), add = TRUE, 
            col = colors[i], lwd = 2)
    }
    # Adding the legend
    legend("topright", legend = paste("Component", 1:n_components), 
           col = colors, lwd = 2, bty = "n")  # bty = "n" removes the box around the legend
  }
}

plot_gmm_predictions <- function(data_nodes, params, bins = 50, 
                                 hist_fill = "lightblue", line_color = "red",
                                 main = "Histogram and Density Estimation of GMM") 
  {
  
  # Function to calculate the predicted density
  calculate_prediction <- function(data, params) {
    n_components <- length(params$means)
    prediction <- rep(0, length(data))
    
    for (i in 1:n_components) {
      prediction <- prediction + (
        params$mixing_probs[i] * dnorm(data,
                                       mean = params$means[i],
                                       sd = params$stds[i]))
    }
    return(prediction)
  }
  
  # Calculate the predicted densities
  densities <- calculate_prediction(unlist(data_nodes), params = params)
  
  # Create a data frame for ggplot
  data_frame <- data.frame(
    x = unlist(data_nodes),
    density = densities
  )
  
  # Plot histogram with density estimation line
  density_plot <- ggplot(data_frame, aes(x = x)) +
    geom_histogram(aes(y = ..density..), bins = bins, fill = hist_fill, 
                   color = "black", alpha = 0.7) +
    geom_line(aes(y = density), color = line_color, size = 1) +
    labs(
      title = main,
      x = "Value",
      y = "Density"
    ) +
    theme_minimal()
  density_plot
  return(density_plot)
}


####################################################################
#####Multivariate#####
generate_multivariate_mixture_gaussian_data <- function(means, covariances, probs, n) {
  # Ensure lengths match
  if (length(means) != length(covariances) || length(means) != length(probs)) {
    stop("Number of means, covariance matrices, and component probabilities must be equal.")
  }
  
  data <- matrix(nrow = n, ncol = length(means[[1]]))
  components <- numeric(n)  # Initialize components vector to store component indices
  
  for (i in 1:n) {
    component <- which.max(rmultinom(1, 1, probs))  # Get the component index
    components[i] <- component                      # Store component number
    data[i, ] <- mvrnorm(1, 
                         mu = means[[component]], 
                         Sigma = covariances[[component]])
  }
  
  return(data)  # Return both data and components
}

plot_mixture_histograms_with_density <- function(data, means, covariances, probs) {
  num_dims <- ncol(data)
  num_components <- length(means)
  n <- nrow(data)
  
  # Calculate component assignment based on probabilities
  components <- numeric(n)
  
  for (i in 1:n) {
    components[i] <- which.max(rmultinom(1, 1, probs))  # Assign component based on probabilities
  }
  
  components <- as.factor(components)  # Convert components to factor for discrete coloring
  
  # Create a data frame for plotting
  df <- as.data.frame(data)
  df$components <- components  # Add the calculated components
  
  # Create a long-form data for easier plotting
  df_long <- melt(df, id.vars = "components", variable.name = "dimension")
  
  # Initialize plot
  plot <- ggplot(df_long, aes(x = value, fill = components)) + 
    geom_histogram(aes(y = after_stat(density)), bins = 30, alpha = 0.6, position = "identity") +
    facet_wrap(~dimension, scales = "free", ncol = 2) +  # Create a grid for each dimension
    scale_fill_manual(values = c("#F8766D", "#00BA38", "#619CFF")) +
    theme_minimal() +
    labs(title = "Histogram of Mixture Components with Density Overlays", x = "Value", y = "Density")
  
  # Display the plot
  print(plot)
}

