#Comparison of different simulation methods in 12 dimensions
##### - Libraries - #####
library(tidyverse)
library(mvtnorm)
library(patchwork)
library(kde1d)
library(rvinecopulib)
mvnorm.mle <- Rfast::mvnorm.mle

##### - Functions - #####
source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/functions.R")
source("scripts/functions/Smania_Jonsson_MICE_simulation.R")

##### - Data - #####
# Load data from simulation_comparison.R
data_total <- read.csv("data/clean/pediatric_data.csv") %>%
  select(age, BW, CREA, LENG, CREF, FRCR, FCRE, CG, SCHW, MDRD, CGSW, MDSW)


##### - Simulations - #####
# Simulation settings
n_sim <- nrow(data_total)
m <- 100
seed_nr <- 89126460

#### - Copulas: spline marginal + parametric copula - ####
set.seed(seed_nr)
large_cop <- estimate_vinecopula_from_data(data_total, family_set = "parametric")
sim_cop <- simulate(large_cop, n = n_sim*m) %>%
  mutate(simulation_nr = rep(1:m, each = n_sim))

#### - Marginals: splines - ####
set.seed(seed_nr)
large_marg <- estimate_vinecopula_from_data(data_total, family_set = "indep")
sim_marg <- simulate(large_marg, n = n_sim*m) %>% 
  mutate(simulation_nr = rep(1:m, each = n_sim))

#### - Multivariate Normal Distribution - ####
set.seed(seed_nr)
mvnorm_param <- mvnorm.mle(as.matrix(data_total))
mvnorm_fit <- list(pdf = function(u) qmvnorm(u, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   pit = function(x) pmvnorm(x, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   rdist = function(n) rmvnorm(n, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma))
sim_mvnorm <- as.data.frame(mvnorm_fit$rdist(m*n_sim)) %>% 
  mutate(simulation_nr = rep(1:m, each = n_sim))
names(sim_mvnorm)[1:ncol(data_total)] <- names(data_total)

#### - Conditional Distribution - ####
set.seed(seed_nr)
sim_cd <- simCovMICE(m = m, data_total, catCovs = NULL) %>% select(-NSIM) %>% 
  mutate(simulation_nr = rep(1:m, each = n_sim))

#### - Bootstrap - ####
set.seed(seed_nr)
sim_bootstrap <- data_total[sample(1:nrow(data_total), size = n_sim*m, replace = TRUE), ] %>% 
  mutate(simulation_nr = rep(1:m, each = n_sim))

#save all simulations
save(data_total, sim_cop, sim_mvnorm, sim_cd, sim_marg, sim_bootstrap, file = "results/comparison/simulated_values_12.Rdata")

### - Visualization - ###
load("results/comparison/simulated_values_12.Rdata")
plot_dat <- get_statistics_multiple_sims(sim_cop, m = m) %>% mutate(simulation = "copula") %>% 
  bind_rows(get_statistics_multiple_sims(sim_cd, m = m) %>% mutate(simulation = "conditional distribution")) %>% 
  bind_rows(get_statistics_multiple_sims(sim_marg, m = m) %>% mutate(simulation = "marginal distribution")) %>% 
  bind_rows(get_statistics_multiple_sims(sim_mvnorm, m = m) %>% mutate(simulation = "mvnormal distribution")) %>% 
  bind_rows(get_statistics_multiple_sims(sim_bootstrap, m = m) %>% mutate(simulation = "bootstrap")) %>% 
  mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
  left_join(get_statistics(data_total) %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate)) %>% 
  group_by(statistic, covariate, simulation) %>% 
  mutate(RMSE = sqrt(mean((value - observed)^2))) %>% ungroup()
plot_dat[, c("covA", "covB")] <- t(apply(as.data.frame(plot_dat[, c("cov1", "cov2")]), 1, sort))

plot_differences <- plot_dat %>% 
  filter(statistic == "correlation") %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE),
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  group_by(statistic, covariate, simulation) %>% 
  summarize(med_rel_error = median(rel_value)) %>% 
  ungroup() %>% 
  filter(med_rel_error > -3) %>% 
  ggplot(aes(y = med_rel_error, x = covariate , color = simulation, shape = simulation)) +
  geom_hline(yintercept = 0, color = "grey35") +  
  geom_hline(yintercept = c(-0.2, 0.2), color = "grey35", linetype = 2) +  
  geom_point() +
  scale_color_manual(values =  create_colors(c("observed", "copula", "marginal distribution", "conditional distribution", "mvnormal distribution", "bootstrap"), 
                                             selected = c("grey", "turquoise", "dark yellow", "pink", "dark green", "midlight blue")), limits = force) +
  labs(x = "Covariate combinations", y = "Median relative error", color = "Method", shape = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_shape_manual(values = c(3, 15:18)) +
  facet_wrap(~ statistic) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6),  strip.background = element_rect(fill = "white"),
        strip.text = element_text( margin = margin( b = 1, t = 1), size = 12))

load("results/figures/manuscript/F2A_performance_3_dimensions.Rdata")

pdf(file = "results/figures/manuscript/F2_performance.pdf", width = 8, height = 7)
(results_plot_relative + labs(x = NULL, tag = "(a)") )/ (plot_differences+ labs(x = "Covariates", tag = "(b)") )
dev.off()


plot_dat %>% filter(statistic == "correlation") %>% 
  group_by(simulation) %>% 
  summarize(med_rel_error = median(rel_value), rmse = mean(rel_value^2)) %>% 
  ungroup()

#### Save copula object ####
save(large_cop, file = "copulas/pediatric_copula.Rdata")

#### Create Suppl table 4 ####
plot_dat %>% 
  filter(statistic %in% c("correlation", "mean")) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE),
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  group_by(statistic, covariate, simulation, RMSE) %>% 
  summarize(med_rel_error = median(rel_value)) %>% 
  ungroup() %>% 
  pivot_longer(c(med_rel_error, RMSE), names_to = "error metric", values_to = "error") %>% 
  pivot_wider(names_from = simulation, values_from = error) %>% 
  arrange(`error metric`, covariate) %>%write_csv(file = "results/12d_performance.csv")
