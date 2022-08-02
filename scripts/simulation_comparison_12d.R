#Comparison of different simulation methods in 12 dimensions
##### - Libraries - #####
library(tidyverse)
library(mvtnorm)

##### - Functions - #####
source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/plot_distributions.R")
source("scripts/functions/functions.R")
source("scripts/functions/Smania_Jonsson_MICE_simulation.R")

##### - Data - #####
# Load data from simulation_comparison.R
data_total <- read.csv("data/clean/pediatric_data.csv") %>%
  select(age, BW, CREA, LENG, CREF, FRCR, FCRE, CG, SCHW, MDRD, CGSW, MDSW)


##### - Simulations - #####
# Simulation settings
n_sim <- 5000
seed_nr <- 89126460

#### - Copulas: spline marginal + parametric copula - ####
set.seed(seed_nr)
large_cop <- estimate_vinecopula_from_data(data_total, family_set = "parametric")
sim_cop <- simulate(large_cop, n = n_sim)

#### - Marginals: splines - ####
set.seed(seed_nr)
large_marg <- estimate_vinecopula_from_data(data_total, family_set = "indep")
sim_marg <- simulate(large_marg, n = n_sim)

#### - Multivariate Normal Distribution - ####
set.seed(seed_nr)
mvnorm_param <- Rfast::mvnorm.mle(as.matrix(data_total))
mvnorm_fit <- list(pdf = function(u) qmvnorm(u, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   pit = function(x) pmvnorm(x, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   rdist = function(n) rmvnorm(n, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma))
sim_mvnorm <- as.data.frame(mvnorm_fit$rdist(n_sim))
names(sim_mvnorm) <- names(data_total)

#### - Conditional Distribution - ####
set.seed(seed_nr)
sim_cd <- simCovMICE(m = 7, data_total, catCovs = NULL) %>% select(-NSIM)

#### - Bootstrap - ####
set.seed(seed_nr)
sim_bootstrap <- data_total[sample(1:nrow(data_total), size = n_sim, replace = TRUE), ]

#save all simulations
#save(data_total, sim_cop, sim_mvnorm, sim_cd, sim_marg, sim_bootstrap, file = "results/comparison/simulated_values_12cov.Rdata")

### - Visualization - ###
#load("results/comparison/simulated_values_12cov.Rdata")

plot_dat <- get_statistics(sim_cop) %>% mutate(simulation = "copula") %>% 
  bind_rows(get_statistics(sim_cd) %>% mutate(simulation = "conditional distribution")) %>% 
  bind_rows(get_statistics(sim_marg) %>% mutate(simulation = "marginal distribution")) %>% 
  bind_rows(get_statistics(sim_mvnorm) %>% mutate(simulation = "mvnormal distribution")) %>% 
  bind_rows(get_statistics(sim_bootstrap) %>% mutate(simulation = "bootstrap")) %>% 
  mutate(value = ifelse(abs(value) == Inf, NA, value)) %>% 
  left_join(get_statistics(data_total) %>% rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(rel_error = abs(value - observed)/abs(observed),
         ratio_diff = value/observed,
         abs_diff = value - observed,
         RMSE = abs(rel_value),
         bias = rel_value) %>% 
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))
plot_dat[, c("covA", "covB")] <- t(apply(as.data.frame(plot_dat[, c("cov1", "cov2")]), 1, sort))

plot_differences <- plot_dat %>% 
  filter(ratio_diff > -3) %>% 
  filter(statistic == "covariance") %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE),
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  ggplot(aes(y = rel_value, x = covariate , color = simulation, shape = simulation)) +
  geom_hline(yintercept = 0, color = "grey35") +  
  geom_hline(yintercept = c(-0.2, 0.2), color = "grey35", linetype = 2) +  
  geom_point() +
  scale_color_manual(values =  create_colors(c("observed", "copula", "marginal distribution", "conditional distribution", "mvnormal distribution", "bootstrap"), 
                                             selected = c("grey", "turquoise", "dark yellow", "pink", "dark green", "midlight blue")), limits = force) +
  labs(x = "Covariate combinations", y = "Relative error", color = "Method", shape = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_shape_manual(values = c(3, 15:18)) +
  facet_wrap(~ statistic) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6),  strip.background = element_rect(fill = "white"),
        strip.text = element_text( margin = margin( b = 1, t = 1), size = 12))


pdf(file = "results/figures/manuscript/F1_performance.pdf", width = 8, height = 3.5)
print(plot_differences)
dev.off()


#results quantified
plot_dat %>% filter(statistic == "covariance") %>% 
  group_by(simulation) %>% 
  summarize(median(rel_value))
