#####################File info#####################
# Code that stores the logic for the application section in dissertation
# The functions stored is:
#     1. run_and_measure_IEM
#     2. run_and_measure_parallel
#     3. run_and_measure_centralised

library(readr)
library(wordcloud)
library(tm)
file_path <- '/Users/bronwynmccall/Documents/GitHub/final_masters_code/'
source(paste(file_path, 'final_IEM_model_univariate.R', sep = ""))
source(paste(file_path, 'final_parallel_model_univariate.R', sep = ""))
source(paste(file_path, 'final_centralised_model_univariate.R', sep = ""))

# Real life application 
m <- 4
data <- read_csv('/Users/bronwynmccall/Downloads/google_play_store/googleplaystore_user_reviews.csv')

# data visualisation
sentiment_polarity <- data$Sentiment_Polarity
hist(sentiment_polarity, 
     main = "Distribution for the sentiment polarity of user reviews
     of Google store applications", 
     xlab = "Sentiment Polarity", 
     col = "lightblue",
     prob = T)
boxplot(sentiment_polarity, 
        main = "Boxplot for the sentiment polarity of user reviews
     of Google store applications",
        col = "lightblue")
null_polarity <- sum(is.na(data$Sentiment_Polarity))
print(null_polarity)
print(null_polarity/length(sentiment_polarity)*100)

# data cleaning 
sentiment_polarity <- as.matrix(na.omit(sentiment_polarity))
length(sentiment_polarity)

# try apply models to data
## Randomly shuffle the rows of the cleaned data
set.seed(123)
shuffled_data <- sentiment_polarity[sample(nrow(sentiment_polarity)), ]
chunk_size <- ceiling(nrow(shuffled_data) / m)
data_chunks <- split(shuffled_data, rep(1:m, each = chunk_size, 
                                        length.out = nrow(shuffled_data)))
chunk_sizes <- sapply(data_chunks, length)
print(paste("Sizes of the chunks:", paste(chunk_sizes, collapse = ", ")))

## Set initial values
means_init <- c(0, -0.5, 0.4)
probs_init <- c(0.4, 0.1, 0.5)
stds_init <- c(0.1, 0.4, 0.4)
params_init <- list("mixing_probs" = probs_init,
                    "means" = means_init,
                    "stds" = stds_init)

## Run the approaches on the use case
### Incremental
run_and_measure_IEM <- function(data_chunks, params_init){
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
incremental_results <- run_and_measure_IEM(data_chunks,
                                                params_init)

### Parallel
run_and_measure_parallel <- function(data_chunks, 
                                     params_init){
  start_time <- Sys.time()
  parallel_distributed_model_results <- GMM_using_parallel(data_chunks,
                                                           params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = parallel_distributed_model_results,
              'time_taken' = total_time))
}
parallel_results <- run_and_measure_parallel(data_chunks,
                                                  params_init)

### Centralised
run_and_measure_centralised <- function(data_chunks,
                                        params_init){
  start_time <- Sys.time()
  data <- c()
  for (i in 1:length(data_chunks)){
    chunk <- data_chunks[[i]]
    data <- c(data, chunk)
  }
  centralised_model_results <- centralised_model(data,
                                                 params_init)
  total_time <- Sys.time() - start_time
  return(list('model_results' = centralised_model_results,
              'time_taken' = total_time))
}
centralised_results <- run_and_measure_centralised(data_chunks,
                                                        params_init)

# Times
print(incremental_results$time_taken)
print(parallel_results$time_taken)
print(centralised_results$time_taken)

# Evaluation
iem_params <- incremental_results$model_results$params
parallel_params <- parallel_results$model_results$params
cent_params <- centralised_results$model_results$params
par(mfrow = c(1, 2))
plot_univariate_gaussian_data(sentiment_polarity, 3, iem_params,
                              main = "GMM components IEM")
plot_univariate_gaussian_data(sentiment_polarity, 3, parallel_params,
                              main = "GMM components Parallel")
par(mfrow = c(1, 1))
plot_univariate_gaussian_data(sentiment_polarity, 3, cent_params,
                              main = "GMM components")
iem_plot <- plot_gmm_predictions(data_chunks, iem_params, 
                     main = "Histogram and Density Estimation of GMM using IEM")
parallel_plot <- plot_gmm_predictions(
  data_chunks, parallel_params,
  main = "Histogram and Density Estimation of GMM using Parallel")
plot_gmm_predictions(data_chunks, cent_params,
                     main = "Histogram and Density Estimation of GMM")
iem_plot + parallel_plot

# Further Applications
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
iem_soft_clust <- calculate_data_belongings(iem_params, sentiment_polarity)
cent_soft_clust <- calculate_data_belongings(cent_params, sentiment_polarity)
parallel_soft_clust <- calculate_data_belongings(parallel_params, 
                                                 sentiment_polarity)
iem_max_col <- max.col(iem_soft_clust, ties.method = "random")
parallel_max_col <- max.col(parallel_soft_clust, ties.method = "random")
data <- data[!is.na(data$Sentiment_Polarity), ]
data$IEM_cluster <- iem_max_col
data$parallel_cluster <- parallel_max_col
### Word clouds
create_wordcloud <- function(data_subset, cluster_number) {
  text <- paste(data_subset$Translated_Review, collapse = " ")
  
  corpus <- Corpus(VectorSource(text))
  
  corpus <- tm_map(corpus, content_transformer(tolower))       
  corpus <- tm_map(corpus, removePunctuation)                 
  corpus <- tm_map(corpus, removeNumbers)                    
  corpus <- tm_map(corpus, removeWords, stopwords("english")) 
  
  tdm <- TermDocumentMatrix(corpus)
  matrix <- as.matrix(tdm)
  word_freq <- sort(rowSums(matrix), decreasing = TRUE)
  
  wordcloud(
    names(word_freq), word_freq,
    max.words = 100,    
    random.order = FALSE,
    colors = brewer.pal(8, "Dark2"),
    main = paste("Word Cloud for Cluster", cluster_number)
  )
}


# Generate word clouds for each cluster
for (cluster in 1:3) {
  cluster_data <- subset(data, IEM_cluster == cluster)
  create_wordcloud(cluster_data, cluster)
}

for (cluster in 1:3) {
  cluster_data <- subset(data, parallel_cluster == cluster)
  create_wordcloud(cluster_data, cluster)
}

