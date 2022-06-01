#Comparison of different simulation methods
##### - Packages - #####
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(EnvStats)
library(mvtnorm)
library(kde1d)
select <- dplyr::select

##### - Data - #####
# Load data
pediatric_data_raw <- read.csv("data/GTV_final_roos_final2.csv", na.strings = ".", colClasses = "numeric")

# Clean data
data_t0 <- pediatric_data_raw %>% 
  filter(TIME == 0) %>% 
  dplyr::select(-c(TIME, DV, MDV)) %>% 
  mutate(AGED = round(AGEW*7)) %>% 
  mutate(age = ifelse(AGED == 0, AGED, AGED - 1))

age_no_crea <- data_t0 %>% group_by(age) %>% 
  summarize(crea_imputed = sum(CREA == 72)/n() == 1) %>% 
  ungroup() %>% 
  filter(crea_imputed) %>% unlist() %>% max()

data_full <- data_t0 %>% 
  filter(!ID %in% c(2004, 2175, 2141)) %>% 
  filter(GRP %in% c(1, 2)) %>%
  filter(age > age_no_crea)

data.table::fwrite(data_full, file = "data/clean/pediatric_data.csv")

data_total <- data_full %>%
  dplyr::select(age, BW, CREA)


##### - Functions - #####
source("scripts/functions/functions.R")
source("scripts/functions/Smania_Jonsson_MICE_simulation.R")
source("scripts/functions/plot_distributions.R")
color_palette <- create_colors(c("Observed", "Copula", "Marginal", "CD", "MVN"), 
                               selected = c("grey", "turquoise", "dark yellow", "pink", "dark green"))

##### - Simulations - #####
# Simulation settings
#set type of simulation study: "full_data" or "random_split" 
simulation_type <- "full_data"


if (simulation_type == "full_data") {
  data_clean <- data_total
  data_test <- data_total
} else if (simulation_type == "random_split") {
  set.seed(7910462)
  ind_train <- sample(1:nrow(data_total), round(nrow(data_total)/2), replace = F)
  ind_test <- -ind_train
  data_clean <- data_total[ind_train, ]
  data_test <- data_total[ind_test, ]
}

#other options
n_sim <- nrow(data_test)
m <- 100
n_statistics <- ncol(data_clean)*5  + choose(ncol(data_clean), 2)
seed_nr <- 794594



# Statistics for observed distribution (truth)
obs_results <- get_statistics(data_test, columns = c("age", "BW", "CREA"))

write.csv(obs_results, file = paste0("results/comparison/obs_results_", simulation_type, ".csv"), row.names = F)

#### - Copulas I: parametric - ####
# Fit distribution
param_crea_raw <- mixtools::normalmixEM(data_clean$CREA, mu = c(20, 60))
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
sim_copula_I <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                           age = marg_age$pdf(sim_unif$age),
                           BW = marg_BW$pdf(sim_unif$BW),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim))

# Get result statistics
copula_I_results <- get_statistics_multiple_sims(sim_copula_I, m, n_statistics, type = "copula_I")

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
sim_copula_II <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                           age = marg_age$pdf(sim_unif$age),
                           BW = marg_BW$pdf(sim_unif$BW),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim))
# Get result statistics
copula_II_results <- get_statistics_multiple_sims(sim_copula_II, m, n_statistics, type = "copula_II")

#### - Copulas III: mix spline marginal + parametric copula - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_spline_marginal(data_clean$age)
marg_BW <- estimate_spline_marginal(data_clean$BW)

data_unif <- data.frame(CREA = tranform_to_uniform(data_clean$CREA, marg_crea),
                        age = tranform_to_uniform(data_clean$age, marg_age),
                        BW = tranform_to_uniform(data_clean$BW, marg_BW))

vine_fit <- vinecop(data_unif, family_set = "parametric", var_types = c("c", "c", "c"))
vine_distribution <- vinecop_dist(vine_fit$pair_copulas, vine_fit$structure, var_types =  c("c", "c", "c"))
# Draw sample
set.seed(seed_nr)
sim_unif <- as.data.frame(rvinecop(n_sim*m, vine_distribution))
names(sim_unif) <- names(data_unif)[1:3]
sim_copula_III <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                             age = marg_age$pdf(sim_unif$age),
                             BW = marg_BW$pdf(sim_unif$BW),
                             type = "simulated", 
                             simulation_nr = rep(1:m, each = n_sim))

copula_III_results <- get_statistics_multiple_sims(sim_copula_III, m, n_statistics, type = "copula_III")


if (simulation_type == "full_data") {
  example_data <- sim_copula_III %>% 
    filter(simulation_nr == 1) %>% 
    dplyr::select(CREA, age, BW) %>% 
    rownames_to_column("ID")
  write.csv(example_data, file = "example_data/pediatrics_example_data.csv", row.names = FALSE, quote = FALSE)
}


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
  dplyr::select(log_age, log_BW, CREA)
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
all_statistics <- copula_III_results %>% 
  #bind_rows(copula_I_results) %>% 
  #bind_rows(copula_II_results) %>% 
  #bind_rows(marginal_I_results) %>% 
  bind_rows(marginal_II_results) %>% 
  bind_rows(mvnorm_results) %>% 
  bind_rows(cd_results)

write.csv(all_statistics, file = paste0("results/comparison/results_", simulation_type, ".csv"), row.names = F)

all_sim <- sim_copula_III %>% 
  mutate(simulation = "copula") %>% 
  bind_rows(sim_marginal %>% 
              mutate(simulation = "marginal") ) %>% 
  bind_rows(sim_mvnorm %>% 
              mutate(simulation = "mvnorm") ) %>% 
  bind_rows(sim_cd %>% 
              mutate(simulation = "cd")) %>% 
  bind_rows(data_clean %>% mutate(simulation = "observed"))


write.csv(all_sim, file = paste0("results/comparison/simulated_values_", simulation_type, ".csv"), row.names = F)

results_plot <-  all_statistics %>% 
  filter(statistic %in% c("mean", "median", "sd", "covariance")) %>% 
  ggplot(aes(y = value, x = statistic, fill = type)) +
  geom_boxplot() +
  geom_hline(data = obs_results %>% mutate(type = "observed") %>% 
               filter(statistic %in% c("mean", "median", "sd", "covariance")) , 
             aes(yintercept = value, linetype = type)) +

  scale_fill_manual(values = create_colors(c("copula_III", "marginal_II", "multivariate normal", "conditional distribution"), 
                                           c("turquoise", "dark yellow", "dark green", "pink"))) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap( ~ covariate , scales = "free", nrow = 2, dir = "v") +
  scale_linetype_manual(values = 2, name = NULL) +
  theme_bw()


results_plot_relative <- all_statistics %>% 
  filter(statistic %in% c("mean", "median", "sd", "covariance")) %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  ggplot(aes(y = rel_value, x = covariate, fill = type)) +
  geom_boxplot() +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_fill_manual(values = create_colors(c("copula_III", "marginal_II", "multivariate normal", "conditional distribution"), 
                                           c("turquoise", "dark yellow", "dark green", "pink"))) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(statistic ~ . , scales = "free", nrow = 2, dir = "v") +
  theme_bw()

pdf(paste0("results/figures/simulation_statistics_", simulation_type,".pdf"), height = 7, width = 10)
print(results_plot)
print(results_plot_relative)
plot_comparison_distribution_sim_obs_generic(sim_cd, data_test, sim_nr = 2, variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Conditional Distributions", pick_color = color_palette[c("CD", "Observed")])
plot_comparison_distribution_sim_obs_generic(sim_copula_III, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Copula", pick_color = color_palette[c("Copula", "Observed")])
plot_comparison_distribution_sim_obs_generic(sim_mvnorm, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Multivariate normal", pick_color = color_palette[c("MVN", "Observed")])
plot_comparison_distribution_sim_obs_generic(sim_marginal, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Marginal", pick_color = color_palette[c("Marginal", "Observed")])
dev.off()

bias_variance <- all_statistics %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(rdiff = (value - observed)/observed) %>% 
  group_by(statistic, covariate, type) %>% 
  summarize(rBias = abs(mean(rdiff)), rRMSE = sqrt(mean(rdiff^2)), 
            .groups = "keep") %>% 
  ungroup()

colors_sim <- create_colors(c("copula_III", "marginal_II", "multivariate normal", "conditional distribution"), 
                            c("turquoise", "dark yellow", "dark green", "pink"))

bias_variance_plot <- bias_variance %>% 
  filter(statistic %in% c("mean", "median", "sd", "covariance")) %>% 
  filter(type != "marginal_II") %>% 
  ggplot(aes(y = rBias, x = rRMSE, color = type, shape = type)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
  geom_point() +
  scale_color_manual(values = colors_sim[order(names(colors_sim))]) +
  facet_wrap(statistic ~ covariate, nrow = 3, dir = "v") +
  theme_bw()

pdf(paste0("results/figures/simulation_statistics_bias_var_", simulation_type,".pdf"), height = 5, width = 8)
print(bias_variance_plot)
dev.off()


#More viz
all_statistics <- read.csv("results/comparison/results_full_data.csv")
all_statistics_rs <- read.csv("results/comparison/results_random_split.csv")

obs_results <- read.csv("results/comparison/obs_results_full_data.csv")
obs_results_rs <- read.csv("results/comparison/obs_results_random_split.csv")


all_stats <- rbind.data.frame(all_statistics %>% mutate(data = "Full data"),
                              all_statistics_rs %>% mutate(data = "Random split"))

obs <- rbind.data.frame(obs_results %>% mutate(data = "Full data"),
                        obs_results_rs %>% mutate(data = "Random split"))

simulation_performance <- all_stats %>% 
  filter(statistic %in% c("mean", "covariance", "sd")) %>% 
  mutate(key = paste(statistic, covariate, sep = "_")) %>% 
  left_join(obs %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(diff_ratio = value/observed) %>% 
  
  filter(data == "Full data") %>% 
  
  ggplot(aes(y = diff_ratio, x = key, fill = type)) +
  geom_boxplot() +
  geom_hline(yintercept = c(1 - 0.2, 1 + 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 1, linetype = 1, color = "black") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_y_continuous(limits = c(NA, 3)) +
 # facet_wrap( ~ data, scales = "free_x", nrow = 2) +
  labs(x = NULL, y = "Relative error (simulated/observed)") +
  scale_fill_manual(values = create_colors(c("copula_III", "marginal_II", "multivariate normal", "conditional distribution"), 
                                           c("turquoise", "dark yellow", "dark green", "pink"))) +
  theme_bw()

pdf("results/figures/simulation_performance_3vars.pdf", height = 5, width = 6)
print(simulation_performance)
dev.off()

results_plot_relative <- all_statistics %>% 
  filter(statistic %in% c("mean", "median", "sd", "covariance")) %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  ggplot(aes(y = rel_value, x = covariate, fill = type)) +
  geom_boxplot() +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap(statistic ~ . , scales = "free", nrow = 2, dir = "v") +
  theme_bw()

all_stats %>% 
  filter(statistic %in% c("mean", "covariance", "sd")) %>% 
  mutate(key = paste(statistic, covariate, sep = "_")) %>% 
  left_join(obs %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(diff_ratio = value/observed) %>% 
  group_by(key, type, data) %>% 
  summarize(rel_value  = mean(rel_value), diff_ratio = mean(diff_ratio)) %>%  ungroup() %>% 
  ggplot(aes(y = diff_ratio, x = key, color = type)) +
  geom_point(aes(shape = type)) +
  geom_hline(yintercept = c(1 - 0.2, 1 + 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 1, linetype = 1, color = "black") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  facet_wrap( ~ data, scales = "free_x", nrow = 2) +
  theme_bw()
