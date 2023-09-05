#Comparison of different simulation methods in 3 dimensions
##### - Libraries - #####
library(tidyverse)
library(rvinecopulib)
library(mvtnorm)
library(kde1d)
library(patchwork)
select <- dplyr::select
mvnorm.mle <- Rfast::mvnorm.mle

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
color_palette <- create_colors(c("observed", "copula", "marginal distribution", "conditional distribution", "mvnormal distribution", "bootstrap"), 
                               selected = c("grey", "turquoise", "dark yellow", "pink", "dark green", "midlight blue"))

##### - Simulations - #####
# Simulation settings
n_sim <- nrow(data_total)
m <- 100
seed_nr <- 794594

# Statistics for observed distribution (truth)
obs_results <- get_statistics(data_total, columns = c("age", "BW", "CREA"))

write.csv(obs_results, file = paste0("results/comparison/obs_results.csv"), row.names = F)

#### - Copulas: spline marginal + parametric copula - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_total$CREA)
marg_age <- estimate_spline_marginal(data_total$age)
marg_BW <- estimate_spline_marginal(data_total$BW)

data_unif <- data.frame(CREA = tranform_to_uniform(data_total$CREA, marg_crea),
                        age = tranform_to_uniform(data_total$age, marg_age),
                        BW = tranform_to_uniform(data_total$BW, marg_BW))

vine_fit <- vinecop(data_unif, family_set = "parametric", var_types = c("c", "c", "c"))
vine_distribution <- vinecop_dist(vine_fit$pair_copulas, vine_fit$structure, var_types =  c("c", "c", "c"))
# Draw sample
set.seed(seed_nr)
sim_unif <- as.data.frame(rvinecop(n_sim*m, vine_distribution))
names(sim_unif) <- names(data_unif)[1:3]
sim_copula <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                         age = marg_age$pdf(sim_unif$age),
                         BW = marg_BW$pdf(sim_unif$BW),
                         simulation_nr = rep(1:m, each = n_sim),
                         simulation = "copula")

copula_results <- get_statistics_multiple_sims(sim_copula, m, type = "copula", columns = c("age", "BW", "CREA"))


#### - Marginals: splines - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_total$CREA)
marg_age <- estimate_spline_marginal(data_total$age)
marg_BW <- estimate_spline_marginal(data_total$BW)

# Draw sample
set.seed(seed_nr)
sim_marginal <- data.frame(CREA = marg_crea$rdist(n_sim*m),
                           age = marg_age$rdist(n_sim*m),
                           BW = marg_BW$rdist(n_sim*m),
                           simulation_nr = rep(1:m, each = n_sim),
                           simulation = "marginal distribution")
# Get result statistics
marginal_results <- get_statistics_multiple_sims(sim_marginal, m, type = "marginal distribution", columns = c("age", "BW", "CREA"))

#### - Multivariate Normal Distribution - ####
# Fit distribution
mvnorm_param <- mvnorm.mle(as.matrix(data_total))
mvnorm_fit <- list(pdf = function(u) qmvnorm(u, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   pit = function(x) pmvnorm(x, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   rdist = function(n) rmvnorm(n, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma))
# Draw sample
set.seed(seed_nr)
sim_mvnorm_raw <- as.data.frame(mvnorm_fit$rdist(n_sim*m))
names(sim_mvnorm_raw) <- names(data_total)
sim_mvnorm <- data.frame(sim_mvnorm_raw,
                         simulation_nr = rep(1:m, each = n_sim),
                         simulation = "mvnormal distribution")

# Get result statistics
mvnorm_results <- get_statistics_multiple_sims(sim_mvnorm, m, type = "mvnormal distribution", columns = c("age", "BW", "CREA"))

#### - Conditional Distribution - ####
# Draw sample
set.seed(seed_nr)
sim_cd <- simCovMICE(m = m, data_total, catCovs = NULL)
names(sim_cd)[names(sim_cd) == "NSIM"] <- "simulation_nr"
sim_cd$simulation <- "conditional distribution"
# Get result statistics
cd_results <- get_statistics_multiple_sims(sim_cd, m, type = "conditional distribution", columns = c("age", "BW", "CREA"))

#### - Bootstrap - ####
# Draw sample
set.seed(seed_nr)
sim_bootstrap_raw <- data_total[sample(1:n_sim, size = n_sim*m, replace = TRUE), ]
sim_bootstrap <- data.frame(sim_bootstrap_raw,
                            simulation_nr = rep(1:m, each = n_sim),
                            simulation = "bootstrap")

# Get result statistics
bootstrap_results <- get_statistics_multiple_sims(sim_bootstrap, m, type = "bootstrap", columns = c("age", "BW", "CREA"))

##### - Plot Results - #####
all_statistics <- copula_results %>% 
  bind_rows(marginal_results) %>% 
  bind_rows(mvnorm_results) %>% 
  bind_rows(cd_results) %>% 
  bind_rows(bootstrap_results)

write.csv(all_statistics, file = paste0("results/comparison/results.csv"), row.names = F)

all_sim <- sim_copula %>% 
  bind_rows(sim_marginal) %>% 
  bind_rows(sim_mvnorm) %>% 
  bind_rows(sim_cd) %>% 
  bind_rows(sim_bootstrap) %>% 
  bind_rows(data_total %>% mutate(simulation = "observed")) 

write.csv(all_sim, file = paste0("results/comparison/simulated_values.csv"), row.names = F)


### - Visualization - ###
all_statistics <- read.csv(file = paste0("results/comparison/results.csv"))
obs_results <- read.csv(file = "results/comparison/obs_results.csv")

# simulation performance Figure 1a
plot_dat_3d <- all_statistics %>% 
  filter(statistic %in% c("mean", "sd", "correlation")) %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  mutate(rel_value = (value - observed)/observed,
         covariate = gsub("_", " - ", covariate, fixed = TRUE))

results_plot_relative <- plot_dat_3d %>% 
  ggplot(aes(y = rel_value, x = covariate, color = type)) +
  geom_vline(xintercept = seq(0.5, 8.5, by = 1), color = "grey95") +
  geom_boxplot(fill = "white") +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_color_manual(values = color_palette, limits = force) +
  labs(x = "Covariates", y = "Relative error", color = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.25, 0.25))) +
  facet_grid(~ statistic , scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

save(results_plot_relative, file = "results/figures/manuscript/F2A_performance_3_dimensions.Rdata")

#values of simulations
plot_dat_3d %>% filter(statistic == "correlation") %>% 
  group_by(type, covariate) %>% 
  summarize(med_rel_error = median(rel_value), rmse = mean(rel_value^2)) 

#Densities figure S2
all_sim <- read.csv("results/comparison/simulated_values.csv")
simulated_data <- all_sim %>% 
  filter(simulation != "observed") %>% 
  mutate(`Body weight` = BW/1000) %>% 
  rename(SCr = CREA, Age = age) %>% 
  select(Age, BW, SCr, `Body weight`, simulation, simulation_nr)

observed_data <- all_sim %>% 
  filter(simulation == "observed") %>% 
  mutate(`Body weight` = BW/1000) %>% 
  rename(SCr = CREA, Age = age) %>% 
  select(Age, BW, SCr, `Body weight`)

plot_data_3d <- simulated_data %>% 
  filter(simulation_nr %in% 1:4) %>% 
  filter(!is.na(simulation))

plot_chars <- list(geom_density2d(aes(color = simulation), bins = 20),
                   geom_density2d(data = observed_data,  alpha = 1, linetype = 2, color = "grey35", bins = 20),
                   scale_y_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, NA)),
                   scale_x_continuous(expand = expansion(mult = c(0, 0.02)), limits = c(0, NA)),
                   facet_grid(~ simulation),
                   scale_color_manual(values = color_palette, limits = force),
                   theme_bw(),
                   theme(legend.position = "none", strip.text = element_text(margin = margin(b = 2, t = 2), size = 14),
                         plot.margin = unit(c(0, 0, 0.1, 0), units = "cm"), strip.background = element_rect(fill = "white")))
age_bw <- plot_data_3d %>% 
  ggplot(mapping = aes(y = `Body weight`, x = Age)) +
  plot_chars

age_scr <-plot_data_3d %>% 
  ggplot(mapping = aes(y = SCr, x = Age)) +
  plot_chars +
  theme(strip.text.x = element_blank())

bw_scr <- plot_data_3d %>% 
  ggplot(mapping = aes(y = SCr, x = `Body weight`)) +
  plot_chars +
  labs(caption = paste0("\nFigure S2: Densities of the three covariate simulations. ",
                        "Grey dashed lines show the observed joint density for each ",
                        "pair of covariates. The solid lines \nrepresent the joint ",
                        "density of a simulated population for each of the five ", 
                        "simulation methods: bootstrap (blue), conditional distributions ", 
                        "(pink), \ncopula (turquoise), marginal distribution (yellow) ", 
                        "and multivariate normal distribution (green).")) +
  theme(strip.text.x = element_blank(), plot.caption = element_text(size = 12, hjust = 0))


pdf("results/figures/manuscript/FS2_comparison_3_dimensions_densities.pdf", width = 12, height = 8.2)
age_bw/age_scr/bw_scr
dev.off()
