#Simulation comparison ranges

# Simulation settings
simulation_type <- "range_split"
range_test <- c(min(data_total$age), median(data_total$age) + 1)
range_train_direction <- ">"

# range_test <- c(28, max(data_total$age))
# range_train_direction <- "<"
range_variable <- "age"

if (range_train_direction == ">") {
  ind_train <- which(data_total[, range_variable] > range_test[2])
} else if (range_train_direction == "<") {
  ind_train <- which(data_total[, range_variable] < range_test[1])
}
ind_test <- which(data_total[, range_variable] >= range_test[1] & data_total[, range_variable] <= range_test[2])
simulation_type <- paste(simulation_type, range_variable,  range_test[1],  range_test[2], sep = "_")

data_clean <- data_total[ind_train, ]
data_test <- data_total[ind_test, ]


obs_results <- get_statistics(data_test, columns = c("age", "BW", "CREA"))
#three methods:
# Fit distribution
param_crea_raw <- mixtools::normalmixEM(data_clean$CREA)
param_crea <- list(estimate = c(mean1 = param_crea_raw$mu[1], sd1 = param_crea_raw$sigma[1], 
                                mean2 = param_crea_raw$mu[2], sd2 = param_crea_raw$sigma[2], 
                                p.mix = param_crea_raw$lambda[2]))
marg_crea <- estimate_parametric_marginal(data_clean$CREA, "normMix", param = param_crea)
marg_age <- estimate_parametric_marginal(c(data_clean$age, data_test$age), "nbinom")
marg_BW <- estimate_parametric_marginal(data_clean$BW, "lnorm")

data_unif <- data.frame(CREA = tranform_to_uniform(data_clean$CREA, marg_crea),
                        age = tranform_to_uniform(data_clean$age, marg_age),
                        BW = tranform_to_uniform(data_clean$BW, marg_BW),
                        age_d = tranform_to_uniform(data_clean$age - 1, marg_age))

vine_fit <- vinecop(data_unif, family_set = "parametric", var_types = c("c", "d", "c"))
vine_distribution <- vinecop_dist(vine_fit$pair_copulas, vine_fit$structure, var_types =  c("c", "d", "c"))

unif_range <- marg_age$pit(range_test)
# Draw sample
set.seed(seed_nr)
sim_unif_raw <- as.data.frame(rvinecop(n_sim*m/(unif_range[2] - unif_range[1]), vine_distribution))
names(sim_unif_raw) <- names(data_unif)[1:3]
#filter range
sim_unif <- sim_unif_raw %>% 
  filter(age > unif_range[1] & age < unif_range[2]) %>% 
  slice(1:(n_sim*m))

sim_copula_I <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                           age = marg_age$pdf(sim_unif$age),
                           BW = marg_BW$pdf(sim_unif$BW),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim)[1:nrow(sim_unif)])

# Get result statistics
copula_I_results <- get_statistics_multiple_sims(sim_copula_I, m, n_statistics, type = "copula_I")

#### - Copulas II: non-parametric - ####
# Fit distribution
marg_crea <- estimate_spline_marginal(data_clean$CREA)
marg_age <- estimate_spline_marginal(c(data_clean$age, data_test$age))
marg_BW <- estimate_spline_marginal(data_clean$BW)

data_unif <- data.frame(CREA = tranform_to_uniform(data_clean$CREA, marg_crea),
                        age = tranform_to_uniform(data_clean$age, marg_age),
                        BW = tranform_to_uniform(data_clean$BW, marg_BW))

vine_fit <- vinecop(data_unif, family_set = "all", var_types = c("c", "c", "c"))
vine_distribution <- vinecop_dist(vine_fit$pair_copulas, vine_fit$structure, var_types =  c("c", "c", "c"))
unif_range <- marg_age$pit(range_test)
# Draw sample
set.seed(seed_nr)
sim_unif_raw <- as.data.frame(rvinecop(n_sim*m/(unif_range[2] - unif_range[1]), vine_distribution))
names(sim_unif_raw) <- names(data_unif)[1:3]
#filter range
sim_unif <- sim_unif_raw %>% 
  filter(age > unif_range[1] & age < unif_range[2]) %>% 
  slice(1:(n_sim*m))

sim_copula_II <- data.frame(CREA = marg_crea$pdf(sim_unif$CREA),
                           age = marg_age$pdf(sim_unif$age),
                           BW = marg_BW$pdf(sim_unif$BW),
                           type = "simulated", 
                           simulation_nr = rep(1:m, each = n_sim)[1:nrow(sim_unif)])
# Get result statistics
copula_II_results <- get_statistics_multiple_sims(sim_copula_II, m, n_statistics, type = "copula_II")

#### - Conditional Distribution - ####
# Draw sample

data_test_var <- data_test
data_test_var[, names(data_test_var) != range_variable] <- NA
data_train_with_seed <- rbind.data.frame(data_test_var, data_clean)
set.seed(seed_nr)
sim_cd_raw <- simCovMICE(m = m, orgCovs = data_train_with_seed, catCovs = NULL, seedCovs = "age", 
                     targetRangeSeedCovs = range_test)
names(sim_cd_raw)[names(sim_cd_raw) == "NSIM"] <- "simulation_nr"

sim_cd <- sim_cd_raw %>% 
  filter(age > range_test[1] & age < range_test[2]) %>% group_by(simulation_nr) %>% 
  slice(1:n_sim) %>% ungroup()

# Get result statistics
cd_results <- get_statistics_multiple_sims(sim_cd, m, n_statistics, type = "conditional distribution")


##### - Plot Results - #####
all_statistics <- copula_I_results %>% 
  bind_rows(copula_II_results) %>%
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

pdf(paste0("results/figures/simulation_statistics_", simulation_type,".pdf"), height = 7, width = 10)
print(results_plot)
plot_comparison_distribution_sim_obs(sim_cd, data_test, sim_nr = 1, 
                                     plot_type = "both", title = "Conditional Distributions")

plot_comparison_distribution_sim_obs(sim_copula_I, data_test, sim_nr = 1, 
                                     plot_type = "both", title = "Copula I")

plot_comparison_distribution_sim_obs(sim_copula_II, data_test, sim_nr = 1, 
                                     plot_type = "both", title = "Copula II")
dev.off()