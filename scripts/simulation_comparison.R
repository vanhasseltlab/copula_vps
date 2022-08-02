#Comparison of different simulation methods
##### - Packages - #####
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(EnvStats)
library(mvtnorm)
library(kde1d)
library(patchwork)
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
color_palette <- create_colors(c("observed", "copula", "marginal distribution", "conditional distribution", "mvnormal distribution", "bootstrap"), 
                               selected = c("grey", "turquoise", "dark yellow", "pink", "dark green", "midlight blue"))

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

#### - Copulas: spline marginal + parametric copula - ####
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
                             simulation_nr = rep(1:m, each = n_sim),
                             simulation = "copula")

copula_results <- get_statistics_multiple_sims(sim_copula_III, m, n_statistics, type = "copula")


if (simulation_type == "full_data") {
  example_data <- sim_copula_III %>% 
    filter(simulation_nr == 1) %>% 
    dplyr::select(CREA, age, BW) %>% 
    rownames_to_column("ID")
  write.csv(example_data, file = "example_data/pediatrics_example_data.csv", row.names = FALSE, quote = FALSE)
}


#### - Marginals: splines - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_spline_marginal(data_clean$age)
marg_BW <- estimate_spline_marginal(data_clean$BW)

# Draw sample
set.seed(seed_nr)
sim_marginal <- data.frame(CREA = marg_crea$rdist(n_sim*m),
                           age = marg_age$rdist(n_sim*m),
                           BW = marg_BW$rdist(n_sim*m),
                           simulation_nr = rep(1:m, each = n_sim),
                           simulation = "marginal distribution")
# Get result statistics
marginal_results <- get_statistics_multiple_sims(sim_marginal, m, n_statistics, type = "marginal distribution")

#### - Multivariate Normal Distribution - ####
# Fit distribution
mvnorm_param <- Rfast::mvnorm.mle(as.matrix(data_clean))
mvnorm_fit <- list(pdf = function(u) qmvnorm(u, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   pit = function(x) pmvnorm(x, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma),
                   rdist = function(n) rmvnorm(n, mean = mvnorm_param$mu, sigma = mvnorm_param$sigma))
# Draw sample
set.seed(seed_nr)
sim_mvnorm_raw <- as.data.frame(mvnorm_fit$rdist(n_sim*m))
names(sim_mvnorm_raw) <- names(data_clean)
sim_mvnorm <- data.frame(sim_mvnorm_raw,
                         simulation_nr = rep(1:m, each = n_sim),
                         simulation = "mvnormal distribution")

# Get result statistics
mvnorm_results <- get_statistics_multiple_sims(sim_mvnorm, m, n_statistics, type = "mvnormal distribution")



#### - Conditional Distribution - ####
# Draw sample
set.seed(seed_nr)
sim_cd <- simCovMICE(m = m, data_clean, catCovs = NULL)
names(sim_cd)[names(sim_cd) == "NSIM"] <- "simulation_nr"
sim_cd$simulation <- "conditional distribution"
# Get result statistics
cd_results <- get_statistics_multiple_sims(sim_cd, m, n_statistics, type = "conditional distribution")

#### - Bootstrap - ####
# Draw sample
set.seed(seed_nr)
sim_bootstrap_raw <- data_clean[sample(1:n_sim, size = n_sim*m, replace = TRUE), ]
sim_bootstrap <- data.frame(sim_bootstrap_raw,
                            simulation_nr = rep(1:m, each = n_sim),
                            simulation = "bootstrap")

# Get result statistics
bootstrap_results <- get_statistics_multiple_sims(sim_bootstrap, m, n_statistics, type = "bootstrap")

##### - Plot Results - #####
all_statistics <- copula_results %>% 
  bind_rows(marginal_results) %>% 
  bind_rows(mvnorm_results) %>% 
  bind_rows(cd_results) %>% 
  bind_rows(bootstrap_results)

write.csv(all_statistics, file = paste0("results/comparison/results_", simulation_type, ".csv"), row.names = F)

all_sim <- sim_copula_III %>% 
  bind_rows(sim_marginal) %>% 
  bind_rows(sim_mvnorm) %>% 
  bind_rows(sim_cd) %>% 
  bind_rows(sim_bootstrap) %>% 
  bind_rows(data_clean %>% mutate(simulation = "observed")) 


write.csv(all_sim, file = paste0("results/comparison/simulated_values_", simulation_type, ".csv"), row.names = F)


##plots
all_statistics <- read.csv(file = paste0("results/comparison/results_", simulation_type, ".csv"))
all_sim <- read.csv(file = paste0("results/comparison/simulated_values_", simulation_type, ".csv"))


results_plot_relative <- all_statistics %>% 
  filter(statistic %in% c("mean", "sd", "covariance")) %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  mutate(rel_value = (value - observed)/observed,
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  ggplot(aes(y = rel_value, x = covariate, fill = type)) +
  geom_vline(xintercept = seq(0.5, 8.5, by = 1), color = "grey95") +
  geom_boxplot() +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_fill_manual(values = color_palette, limits = force) +
  labs(x = "Covariates", y = "Relative error", fill = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90), expand = expansion(mult = c(0.25, 0.25))) +
  facet_grid(~ statistic , scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank())

pdf("results/figures/simulation3d_performance.pdf", height = 4, width =7)
print(results_plot_relative)
dev.off()

pdf(paste0("results/figures/simulation_statistics_", simulation_type,".pdf"), height = 7, width = 10)
print(results_plot_relative)
plot_comparison_distribution_sim_obs_generic(sim_cd, data_test, sim_nr = 2, variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Conditional Distributions", pick_color = color_palette[c("conditional distribution", "observed")])
plot_comparison_distribution_sim_obs_generic(sim_copula_III, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Copula", pick_color = color_palette[c("copula", "observed")])
plot_comparison_distribution_sim_obs_generic(sim_mvnorm, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Multivariate normal", pick_color = color_palette[c("mvnormal distribution", "observed")])
plot_comparison_distribution_sim_obs_generic(sim_marginal, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                     plot_type = "both", title = "Marginal", pick_color = color_palette[c("marginal distribution", "observed")])
plot_comparison_distribution_sim_obs_generic(sim_bootstrap, data_test, sim_nr = 2,  variables = c("age", "BW", "CREA"), 
                                             plot_type = "both", title = "Bootstrap", pick_color = color_palette[c("bootstrap", "observed")])
dev.off()

#alternative
results_plot_relative <- all_statistics %>% 
  filter(statistic %in% c("mean", "sd", "covariance")) %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  mutate(rel_value = (value - observed)/observed,
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
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

#table with results
summary_stats <- all_statistics %>% 
  left_join(obs_results %>% rename(observed = value)) %>%
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(diff_ratio = value/observed) %>% 
  group_by(statistic, covariate, type) %>% 
  summarize(median_value = median(value), rel_error = mean(rel_value), 
            diff_ratio = median(diff_ratio), RMSE = sqrt(mean(rel_value^2)),
            bias = mean(rel_value)) %>% 
  ungroup() %>% 
  mutate(error_prop = 1 - diff_ratio) %>% 
  arrange(type, covariate) %>% 
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub("\\w+.\\_", "", covariate))


#results similar to Smania and Jonsson
plot_covariance_SJ <- summary_stats %>% 
  pivot_longer(cols = c(bias, RMSE), names_to = "error_type", values_to = "error") %>% 
  filter(statistic == "covariance") %>% 
  ggplot(aes(y = cov1, x = cov2)) +
  geom_tile(aes(fill = error)) +
  geom_text(aes(label = round(error, 2))) +
  scale_fill_gradientn(colours = c("#001158", "white", "#001158"), limits = c(-0.3, 0.3)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  facet_grid(error_type ~ type) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid = element_blank())

plot_marginal_SJ <- summary_stats %>% 
  pivot_longer(cols = c(bias, RMSE), names_to = "error_type", values_to = "error") %>% 
  filter(statistic != "covariance") %>% 
  ggplot(aes(y = statistic, x = covariate)) +
  geom_tile(aes(fill = error)) +
  geom_text(aes(label = round(error, 2))) +
  scale_fill_gradientn(colours = c("#001158", "white", "#001158"), limits = c(-1, 1)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  facet_grid(error_type ~ type) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid = element_blank())

pdf("results/figures/SJ_plots.pdf", width = 10, height = 6)
print(plot_marginal_SJ)
print(plot_covariance_SJ)
dev.off()


summary_stats %>% 
  filter(statistic != "covariance") %>% 
  ggplot(aes(y = statistic, x = covariate)) +
  geom_tile(aes(fill = bias)) +
  geom_text(aes(label = round(bias, 2))) +
  scale_fill_gradientn(colours = c("#f46e32", "white", "#f46e32"), limits = c(-1.1, 1.1)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  facet_grid( ~ type) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid = element_blank())

summary_stats %>% 
  filter(statistic == "covariance") %>% 
  ggplot(aes(y = cov1, x = cov2)) +
  geom_tile(aes(fill = RMSE)) +
  geom_text(aes(label = round(RMSE, 2))) +
  scale_fill_gradientn(colours = c("#f46e32", "white", "#f46e32"), limits = c(-1.1, 1.1)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  facet_grid(~ type) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid = element_blank())


#################


bias_variance <- all_statistics %>% 
  left_join(obs_results %>% rename(observed = value)) %>% 
  mutate(rdiff = (value - observed)/observed) %>% 
  group_by(statistic, covariate, type) %>% 
  summarize(rBias = abs(mean(rdiff)), rRMSE = sqrt(mean(rdiff^2)), 
            .groups = "keep") %>% 
  ungroup()

bias_variance_plot <- bias_variance %>% 
  filter(statistic %in% c("mean", "sd", "covariance")) %>% 
  ggplot(aes(y = rBias, x = rRMSE, color = type, shape = type)) +
  geom_vline(xintercept = 0, linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 2, color = "grey65") +
  geom_point() +
  scale_color_manual(values = color_palette, limits = force) +
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
  #scale_y_continuous(limits = c(NA, 3)) +
 # facet_wrap( ~ data, scales = "free_x", nrow = 2) +
  labs(x = NULL, y = "Relative error (simulated/observed)") +
  scale_fill_manual(values = color_palette, limits = force) +
  theme_bw()

pdf("results/figures/simulation_performance_3vars.pdf", height = 5, width = 6)
print(simulation_performance)
dev.off()


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
  scale_color_manual(values = color_palette, limits = force) +
  facet_wrap( ~ data, scales = "free_x", nrow = 2) +
  theme_bw()




#alternative
simulations_3d <- read.csv("results/comparison/simulated_values_full_data.csv")
simulated_data <- simulations_3d %>% 
  filter(simulation != "observed") %>% 
  mutate(`Body weight` = BW/1000) %>% 
  rename(SCr = CREA, Age = age) %>% 
  select(Age, BW, SCr, `Body weight`, simulation, simulation_nr)

observed_data <- simulations_3d %>% 
  filter(simulation == "observed") %>% 
  mutate(`Body weight` = BW/1000) %>% 
  rename(SCr = CREA, Age = age) %>% 
  select(Age, BW, SCr, `Body weight`)

plot_data_3d <- simulated_data %>% 
  filter(simulation_nr %in% 1:4) %>% 
  filter(!is.na(simulation))

plot_chars <- list(geom_density2d(aes(color = simulation), bins = 20),
                   #geom_point(aes(color = simulation), alpha = 0.5, stroke = 0, shape = 16),
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
  theme(strip.text.x = element_blank())


pdf("results/figures/comparison_density_methods.pdf", width = 12, height = 8)
age_bw/age_scr/bw_scr
dev.off()



abw <- observed_data %>% 
  ggplot(aes(x = Age, y = `Body weight`)) +
  geom_density2d(color = "black", bins = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(-1, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(-50, NA)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())
abw

asc <- observed_data %>% 
  ggplot(aes(x = Age, y = `SCr`)) +
  geom_density2d(color = "black", bins = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 1.2))) +
  scale_x_continuous(expand = expansion(mult = c(0, 0))) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())

scbw <- observed_data %>% 
  ggplot(aes(x = `Body weight`, y = `SCr`)) +
  geom_density2d(color = "black", bins = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0)), limits = c(-10, 100)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), limits = c(-1, 9)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())

abw + asc + scbw


rs <- data_full %>% 
  ggplot(aes(x = MDRD, y = `CG`)) +
  geom_density2d(color = "black", bins = 10) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(-20, NA)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.01)), limits = c(0, NA)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  theme(panel.grid = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.border = element_blank())


pdf("presentation/figures/densities_conceptual.pdf", width = 4.5, height = 1.6)
abw + scbw + rs


abw +  geom_point(color = "blue") + scbw + rs


dev.off()


