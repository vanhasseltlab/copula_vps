#Comparison of different simulation methods
##### - Packages - #####
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(EnvStats)
library(mvtnorm)
library(kde1d)

##### - Data - #####
# Load data
pediatric_data_raw <- read.csv("data/GTV_final_roos_final2.csv", na.strings = ".", colClasses = "numeric")

# Clean data
data_t0 <- pediatric_data_raw %>% 
  filter(TIME == 0) %>% 
  select(-c(TIME, DV, MDV)) %>% 
  mutate(AGED = round(AGEW*7)) %>% 
  mutate(age = ifelse(AGED == 0, AGED, AGED - 1))

age_no_crea <- data_t0 %>% group_by(age) %>% 
  summarize(crea_imputed = sum(CREA == 72)/n() == 1) %>% 
  ungroup() %>% 
  filter(crea_imputed) %>% unlist() %>% max()

data_total <- data_t0 %>% 
  filter(!ID %in% c(2004, 2175, 2141)) %>% 
  filter(GRP %in% c(1, 2)) %>%
  filter(age > age_no_crea) %>%
  select(age, BW, CREA)


##### - Functions - #####
source("scripts/functions/functions.R")
source("scripts/functions/Smania_Jonsson_MICE_simulation.R")

##### - Simulations - #####
# Simulation settings
#set type of simulation study: "full_data", "random_split" or "range_split" set range for range split
simulation_type <- "range_split"
range <- c(median(data_total$age) + 1, max(data_total$age))
range_variable <- "age"


if (simulation_type == "full_data") {
  data_clean <- data_total
  data_test <- data_total
} else {
  if (simulation_type == "random_split") {
    set.seed(7910462)
    ind_train <- sample(1:nrow(data_total), round(nrow(data_total)/2), replace = F)
    }
  if (simulation_type == "range_split") {
    ind_train <- which(data_total[, range_variable] >= range[1] & data_total[, range_variable] <= range[2])
    simulation_type <- paste(simulation_type, range_variable,  range[1],  range[2], sep = "_")
    }
  data_clean <- data_total[ind_train, ]
  data_test <- data_total[-ind_train, ]
}

#other options
n_sim <- nrow(data_test)
m <- 30
n_statistics <- ncol(data_clean)*5  + choose(ncol(data_clean), 2)
seed_nr <- 794594



# Statistics for observed distribution (truth)
obs_results <- get_statistics(data_test, columns = c("age", "BW", "CREA"))

#### - Copulas I: parametric - ####
# Fit distribution
param_crea_raw <- mixtools::normalmixEM(data_clean$CREA)
param_crea <- list(estimate = c(mean1 = param_crea_raw$mu[1], sd1 = param_crea_raw$sigma[1], 
                                mean2 = param_crea_raw$mu[2], sd2 = param_crea_raw$sigma[2], 
                                p.mix = param_crea_raw$lambda[2]))
marg_crea <- estimate_parametric_marginal(data_clean$CREA, "normMix", param = param_crea)
marg_age <- estimate_parametric_marginal(data_clean$age, "nbinom")
marg_BW <- estimate_parametric_marginal(data_clean$BW, "lnorm")

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


range_set <- sim_original[, range_variable] < range[1]
sim_train <- sim_original[range_set, ]
hist(sim_train$age, breaks = 30)
range_results <- get_statistics_multiple_sims(sim_train, m, n_statistics, type = "copula_I")

# Get result statistics
copula_I_results <- get_statistics_multiple_sims(sim_original, m, n_statistics, type = "copula_I")

#### - Copulas II: non-parametric - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_spline_marginal(data_clean$age)
marg_BW <- estimate_spline_marginal(data_clean$BW)

data_unif <- data.frame(CREA = tranform_to_uniform(data_clean$CREA, marg_crea),
                        age = tranform_to_uniform(data_clean$age, marg_age),
                        BW = tranform_to_uniform(data_clean$BW, marg_BW))

vine_fit <- vinecop(data_unif, family_set = "all", var_types = c("c", "c", "c"))
vine_distribution <- vinecop_dist(vine_fit$pair_copulas, vine_fit$structure, var_types =  c("c", "c", "c"))
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
copula_II_results <- get_statistics_multiple_sims(sim_original, m, n_statistics, type = "copula_II")


#### - Marginals I: parametric - ####
# Fit distribution
param_crea_raw <- mixtools::normalmixEM(data_clean$CREA)
param_crea <- list(estimate = c(mean1 = param_crea_raw$mu[1], sd1 = param_crea_raw$sigma[1], 
                                mean2 = param_crea_raw$mu[2], sd2 = param_crea_raw$sigma[2], 
                                p.mix = param_crea_raw$lambda[2]))
marg_crea <- estimate_parametric_marginal(data_clean$CREA, "normMix", param = param_crea)
marg_age <- estimate_parametric_marginal(data_clean$age, "nbinom")
marg_BW <- estimate_parametric_marginal(data_clean$BW, "lnorm")
# Draw sample
set.seed(seed_nr)
sim_marginal <- data.frame(CREA = marg_crea$rdist(n_sim*m),
                           age = marg_age$rdist(n_sim*m),
                           BW = marg_BW$rdist(n_sim*m),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim))
# Get result statistics
marginal_I_results <- get_statistics_multiple_sims(sim_marginal, m, n_statistics, type = "marginal_I")

#### - Marginals II: splines - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_spline_marginal(data_clean$age)
marg_BW <- estimate_spline_marginal(data_clean$BW)
# Draw sample
set.seed(seed_nr)
sim_marginal <- data.frame(CREA = marg_crea$rdist(n_sim*m),
                           age = marg_age$rdist(n_sim*m),
                           BW = marg_BW$rdist(n_sim*m),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim))
# Get result statistics
marginal_II_results <- get_statistics_multiple_sims(sim_marginal, m, n_statistics, type = "marginal_II")

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
# Draw sample
set.seed(seed_nr)
sim_cd <- simCovMICE(m = m, data_clean, catCovs = NULL)
names(sim_cd)[names(sim_cd) == "NSIM"] <- "simulation_nr"
# Get result statistics
cd_results <- get_statistics_multiple_sims(sim_cd, m, n_statistics, type = "conditional distribution")

##### - Plot Results - #####
all_statistics <- copula_I_results %>% 
  bind_rows(copula_II_results) %>% 
  bind_rows(marginal_I_results) %>% 
  bind_rows(marginal_II_results) %>% 
  bind_rows(mvnorm_results) %>% 
  bind_rows(cd_results)

write.csv(all_statistics, file = paste0("results/comparison/results_", simulation_type, ".csv"), row.names = F)

results_plot <-  all_statistics %>% 
  filter(statistic %in% c("mean", "median", "sd") | covariate == "all") %>% 
  ggplot(aes(y = value, x = type, color = covariate)) +
  geom_boxplot() +
  geom_hline(data = obs_results %>% mutate(type = "observed") %>% 
               filter(statistic %in% c("mean", "median", "sd") | covariate == "all") , 
             aes(yintercept = value, linetype = type)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(statistic ~ covariate , scales = "free_y", nrow = 3, dir = "v") +
  scale_linetype_manual(values = 2, name = NULL) +
  theme_bw()

pdf("results/figures/simulation_statistics", simulation_type,".pdf", height = 7, width = 10)
print(results_plot)
dev.off()

bias_variance <- all_statistics %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(rdiff = (value - observed)/observed) %>% 
  group_by(statistic, covariate, type) %>% 
  summarize(rBias = abs(mean(rdiff)), rRMSE = sqrt(mean(rdiff^2)), .groups = "keep") %>% 
  ungroup()

colors_sim <- c("#D41212", "#2E9BFB", "#1C69AD", "#E2943C", "#B16919", "#73BB2D")

bias_variance_plot <- bias_variance %>% 
  filter(statistic %in% c("mean", "median", "sd") | covariate == "all") %>% 
  ggplot(aes(y = rBias, x = rRMSE, color = type, shape = type)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
  geom_point() +
  scale_color_manual(values = colors_sim) +
  facet_wrap(statistic ~ covariate, nrow = 3, dir = "v") +
  theme_bw()

pdf(paste0("results/figures/simulation_statistics_bias_var", simulation_type,".pdf"), height = 5, width = 8)
print(bias_variance_plot)
dev.off()