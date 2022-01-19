#Comparison of different simulation methods
##### - Packages - #####
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(mvtnorm)
library(kde1d)

##### - Data - #####
# Clean data
age_no_crea <-data_t0 %>% group_by(age) %>% 
  summarize(crea_imputed = sum(CREA == 72)/n() == 1) %>% 
  ungroup() %>% 
  filter(crea_imputed) %>% unlist() %>% max()

data_clean <- data_t0 %>% 
  filter(!ID %in% c(2004, 2175, 2141)) %>% 
  filter(GRP %in% c(1, 2)) %>%
  filter(age > age_no_crea) %>%
  select(age, BW, CREA)


##### - Functions - #####
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
               dist = param)
  return(marg)
}
#estimate parametric distribution using kde1d package
estimate_parametric_marginal <- function(covariate, distribution) {
  param <- fitdist(covariate, distribution)
  marg <- list(pdf = function(u) {
    do.call(eval(paste0("q", distribution)), c(list(p = u), param$estimate))},
    pit = function(x) {
      do.call(eval(paste0("p", distribution)), c(list(q = x), param$estimate))}, 
    rdist = function(n) {
      do.call(eval(paste0("r", distribution)), c(list(n = n), param$estimate))},
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


##### - Simulations - #####
# Simulation settings
n_sim <- nrow(data_clean)
m <- 30
n_statistics <- ncol(data_clean)*5  + choose(ncol(data_clean), 2)
seed_nr <- 794594
# Statistics for observed distribution (truth)
obs_results <- get_statistics(data_clean, columns = c("age", "BW", "CREA"))

#### - Copulas - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_parametric_marginal(data_clean$age, "nbinom")
marg_BW <- estimate_spline_marginal(data_clean$BW)

data_unif <- data.frame(CREA = tranform_to_uniform(data_clean$CREA, marg_crea),
                        age = tranform_to_uniform(data_clean$age, marg_age),
                        BW = tranform_to_uniform(data_clean$BW, marg_BW),
                        age_d = tranform_to_uniform(data_clean$age - 1, marg_age))

vine_fit <- vinecop(data_unif, family_set = "parametric", var_types = c("c", "d", "c"))
vine_distribution <- vinecop_dist(vine_fit$pair_copulas, vine_fit$structure, var_types =  c("c", "d", "c"))
# Draw sample
set.seed(seed_nr)
sim_unif <- as.data.frame(rvinecop(n_sim*m, vine_distribution))
names(sim_unif) <- names(data_unif)[1:3]
sim_original <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                           age = marg_age$pdf(sim_unif$age),
                           BW = marg_BW$pdf(sim_unif$BW),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim))
# Get result statistics
copula_results <- get_statistics_multiple_sims(sim_original, m, n_statistics, type = "copula")

#### - Marginals - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_parametric_marginal(data_clean$age, "nbinom")
marg_BW <- estimate_spline_marginal(data_clean$BW)
# Draw sample
set.seed(seed_nr)
sim_marginal <- data.frame(CREA = marg_crea$rdist(n_sim*m),
                           age = marg_age$rdist(n_sim*m),
                           BW = marg_BW$rdist(n_sim*m),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim))
# Get result statistics
marginal_results <- get_statistics_multiple_sims(sim_marginal, m, n_statistics, type = "marginal")

#### - Multivariate Normal Distribution - ####
# Fit distribution
data_clean_trans <- data_clean %>% 
  mutate(log_BW = log(BW),
         log_age = log(age)) %>% 
  select(log_age, log_BW, CREA)
mvnorm_param <- Rfast::mvnorm.mle(as.matrix(data_clean_trans))
mvnorm_fit <- list(pdf = function(u) qmvnorm(u, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   pit = function(x) pmvnorm(x, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   rdist = function(n) rmvnorm(n, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma))
# Draw sample
set.seed(seed_nr)
sim_mvnorm_raw <- as.data.frame(mvnorm_fit$rdist(n_sim*m))
names(sim_mvnorm_raw) <- names(data_clean_trans)
sim_mvnorm <- data.frame(CREA = sim_mvnorm_raw$CREA,
                         age = exp(sim_mvnorm_raw$log_age),
                         BW = exp(sim_mvnorm_raw$log_BW),
                         type = "simulated", 
                         simulation_nr = rep(1:m, each = n_sim))
# Get result statistics
mvnorm_results <- get_statistics_multiple_sims(sim_mvnorm, m, n_statistics, type = "multivariate normal")

#### - Conditional Distribution - ####
# Import MICE simulation function simCovMICE
source("scripts/Smania_Jonsson_MICE_simulation.R")
# Draw sample
set.seed(seed_nr)
sim_cd <- simCovMICE(m = m, data_clean, catCovs = NULL)
names(sim_cd)[names(sim_cd) == "NSIM"] <- "simulation_nr"
# Get result statistics
cd_results <- get_statistics_multiple_sims(sim_cd, m, n_statistics, type = "conditional distribution")

##### - Plot Results - #####
results_plot <- copula_results %>% 
  bind_rows(marginal_results) %>% 
  bind_rows(mvnorm_results) %>% 
  bind_rows(cd_results) %>% 
  filter(statistic %in% c("mean", "median", "sd") | covariate == "all") %>% 
  ggplot(aes(y = value, x = type, color = covariate)) +
  geom_boxplot() +
  geom_hline(data = obs_results %>% mutate(type = "observed") %>% 
               filter(statistic %in% c("mean", "median", "sd") | covariate == "all") , 
             aes(yintercept = value, linetype = type)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(statistic ~covariate , scales = "free_y", nrow = 3, dir = "v") +
  scale_linetype_manual(values = 2, name = NULL) +
  theme_bw()

pdf("results/figures/simulation_statitsics.pdf", height = 7, width = 10)
print(results_plot)
dev.off()