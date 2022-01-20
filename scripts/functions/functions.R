

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
  univariate <- as.data.frame(apply(data_set[, columns], 2, stats))
  
  multivariate <- reshape2::melt(replace(covariance_mat, lower.tri(covariance_mat, TRUE), NA), na.rm = TRUE) %>% 
    mutate(statistic = paste("cov", Var1, Var2, sep = "_"), covariate = "all") %>% 
    select(statistic, covariate, value)
  
  long_format <- univariate %>% 
    rownames_to_column("statistic") %>% 
    pivot_longer(-statistic, names_to = "covariate") %>% 
    bind_rows(multivariate)
  
  return(long_format)
}

#calculate statistics for comparison for each set of simulations
get_statistics_multiple_sims <- function(data_set, m, n_statistics, columns = NULL, type = NULL) {
  full_results <- as.data.frame(matrix(nrow = m*n_statistics, ncol = 3))
  names(full_results) <- c("statistic", "covariate", "value")
  full_results$simulation_nr <- rep(1:m, each = n_statistics)
  for(i in unique(data_set$simulation_nr)) {
    sim_results <- get_statistics(data_set[data_set$simulation_nr == i, ], columns = c("age", "BW", "CREA"))
    full_results[full_results$simulation_nr == i, c("statistic", "covariate", "value")] <- sim_results
  }
  if (!is.null(type)) {
    full_results$type <- type
  }
  
  return(full_results)
}

#estimate splines using fitdistrplus package and actuar package for extra distribution options
estimate_spline_marginal <- function(covariate, xmin = NaN) {
  param <- kde1d(covariate, xmin = xmin)
  marg <- list(pdf = function(u) qkde1d(u, param),
               pit = function(x) pkde1d(x, param),
               rdist = function(n) rkde1d(n, param),
               density = function(x) dkde1d(x, param),
               dist = param)
  return(marg)
}

#estimate parametric distribution using kde1d package
estimate_parametric_marginal <- function(covariate, distribution, param = NULL) {
  if (is.null(param)) {
    param <- fitdist(covariate, distribution)
  }
  marg <- list(pdf = function(u) {
    do.call(eval(paste0("q", distribution)), c(list(p = u), param$estimate))},
    pit = function(x) {
      do.call(eval(paste0("p", distribution)), c(list(q = x), param$estimate))}, 
    rdist = function(n) {
      do.call(eval(paste0("r", distribution)), c(list(n = n), param$estimate))},
    density = function(x) {
      do.call(eval(paste0("d", distribution)), c(list(x = x), param$estimate))},
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


#check a marginal plot fit (base R)
check_fit_plot <- function(x, marg_density) {
  x_dens <- seq(min(x, na.rm = T), max(x, na.rm = TRUE), by = 1)
  hist(x, probability = TRUE, breaks = 30)
  lines(x_dens, marg_density(x_dens), col = "red", lwd = 2)
}
