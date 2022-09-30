#helper functions for covariate simulation and performance evaluation
#calculate statistics for comparison
get_statistics <- function(data_set, columns = NULL) {
  if(is.null(columns)) {
    columns <- colnames(data_set)
  }
  stats <- function(x) {
    x <- x[!is.na(x)]
    c(mean = mean(x), median = median(x), sd = sd(x), min = min(x), max = max(x))
  }
  covariance_mat <- cov(data_set[, columns], use = "pairwise.complete.obs")
  cor_mat <- cor(data_set[, columns], use = "pairwise.complete.obs")
  univariate <- as.data.frame(apply(data_set[, columns], 2, stats))
  
  multivariate <- reshape2::melt(replace(covariance_mat, lower.tri(covariance_mat, TRUE), NA), na.rm = TRUE) %>% 
    mutate(statistic = "covariance", covariate = paste(Var1, Var2, sep = "_")) %>% 
    dplyr::select(statistic, covariate, value) %>% 
    bind_rows(reshape2::melt(replace(cor_mat, lower.tri(cor_mat, TRUE), NA), na.rm = TRUE) %>% 
                mutate(statistic = "correlation", covariate = paste(Var1, Var2, sep = "_")))
  
  long_format <- univariate %>% 
    rownames_to_column("statistic") %>% 
    pivot_longer(-statistic, names_to = "covariate") %>% 
    bind_rows(multivariate) %>% 
    mutate(Var1 = as.character(Var1), Var2 = as.character(Var2))
    
  
  return(long_format)
}

#calculate statistics for comparison for each set of simulations
get_statistics_multiple_sims <- function(data_set, m, columns = NULL, type = NULL) {
  data_set <- data_set[, colSums(is.na(data_set)) != nrow(data_set)]
  
  if (is.null(columns)) {
    columns <- setdiff(colnames(data_set), c("simulation_nr", "simulation"))
  }
  n_statistics <- length(columns)*5  + 2*choose(length(columns), 2)
  full_results <- as.data.frame(matrix(nrow = m*n_statistics, ncol = 5))
  names(full_results) <- c("statistic", "covariate", "value", "cov1", "cov2")
  full_results$simulation_nr <- rep(1:m, each = n_statistics)
  
  
  
  for(i in unique(data_set$simulation_nr)) {
    sim_results <- get_statistics(data_set[data_set$simulation_nr == i, ], columns = columns)
    full_results[full_results$simulation_nr == i, c("statistic", "covariate", "value",  "cov1", "cov2")] <- sim_results
  }
  if (!is.null(type)) {
    full_results$type <- type
  }
  
  return(full_results)
}

#estimate splines using fitdistrplus package and actuar package for extra distribution options
estimate_spline_marginal <- function(covariate, xmin = NaN) {
  covariate <- covariate[!is.na(covariate)]
  param <- kde1d(covariate, xmin = xmin)
  marg <- list(pdf = function(u) qkde1d(u, param),
               pit = function(x) pkde1d(x, param),
               rdist = function(n) rkde1d(n, param),
               density = function(x) dkde1d(x, param),
               dist = param)
  return(marg)
}

#transform a variable to a uniform distribution using it's marginal distribution
tranform_to_uniform <- function(covariate, marg_dist) {
  if (any(is.na(covariate))) {
    ind_not_na <- which(!is.na(covariate))
    unif_vector <- numeric(length = length(covariate))
    unif_vector[] <- NA
    unif_vector[ind_not_na] <- marg_dist$pit(covariate[ind_not_na])
  } else {
    unif_vector <- marg_dist$pit(covariate)
  }
  return(unif_vector)
}


#create color data frame
create_colors <- function(labels = NULL, selected = 1:length(labels), demo = FALSE) {
  colors <- c("#F3D4DF", "#E7AAC0", "#DB7FA1", "#CF5581", "#C32A62", "#B70043", "#1784E4", 
              "#1784E4", "#E4AB01", "#04715F", "#B52807", "#C3C3C3", "#3ABAC1", "#C37121", 
              "#26B72C", "#8DD2FF", "#9DF7A1", "#F7F18B", "#F3A492", "#083A9C", "#DB7FA1",
              "#001158", "#f46e32", "#969696", "#F46E32", "black")
  names(colors) <- c("grey pink", "light pink", "midlight pink", "mid pink", "middark pink", 
                     "dark pink", "midlight blue", "blue", "dark yellow", "dark green", "red", 
                     "light grey", "turquoise", "brown", "green", "light blue", "light green", 
                     "yellow", "peach", "dark blue", "pink", "leiden blue", "leiden orange", "grey", "orange", "black")
  color_palette <- colors[selected]
  if (!is.null(labels)) {
    names(color_palette) <- labels[1:length(selected)]
  }
  
  if (!demo) {
    return(color_palette)
  }
  
  demo_colors <- function(color_df, blocks) {
    color_plot <- color_df %>% 
      mutate(x_demo = rep(1:blocks[1], each = blocks[2]),
             y_demo = rep(1:blocks[2], times = blocks[1])) %>% 
      ggplot(aes(x = x_demo, y = y_demo)) +
      geom_tile(aes(fill = color)) +
      scale_fill_identity() +
      geom_text(aes(label = color, color = font.color)) +
      scale_color_identity() +
      theme_void()
    return(color_plot)
  }
  
  return(demo_colors(data.frame(color = color_palette, font.color = "black"), c(13, 2)))
}

