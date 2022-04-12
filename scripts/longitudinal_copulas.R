#Longitudinal analysis
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(EnvStats)
library(mvtnorm)
library(kde1d)
library(lme4)
select <- dplyr::select

source("scripts/functions/functions.R")
source("scripts/functions/plot_distributions.R")
source("scripts/functions/estimate_vinecopula_from_data.R")
#read in pregnancy data
data_pregnancy_raw <- read.csv("data/Jig - pregnancy file August 2012 minus coag info.csv", 
                               row.names = NULL, na.strings = c("")) %>% 
  mutate(Neutrophils = as.numeric(str_replace(Neutrophils, ":", ".")))


#Data preparation: Remove (most) imputed values
data_reduced <- data_pregnancy_raw
remove_rows <- numeric()
index_gest <- c(grep("gest", colnames(data_reduced)), grep("dat1", colnames(data_reduced)))
previous_meas <- as.character(data_reduced[1, -index_gest])
previous_meas[is.na(previous_meas)] <- "-999"
for (i in 2:nrow(data_reduced)) {
  current_meas <- as.character(data_reduced[i, -index_gest])
  current_meas[is.na(current_meas)] <- "-999"
  if (all(current_meas == previous_meas)) {
    remove_rows <- c(remove_rows, i)
  }
  previous_meas <- current_meas
}
data_reduced <- data_reduced[-remove_rows, ]

write.csv(data_reduced, file = "data/clean/pregnancy_reduced.csv", row.names = FALSE,
          quote = FALSE)

###
#continue with data_reduced
data_reduced <- read.csv("data/clean/pregnancy_reduced.csv", row.names = NULL)

### Using polynomials
#remove all observations with baby (only looking at pregnancy)
data_clean <- data_reduced %>% 
  filter(BABY == 0)


#Polynomials for multiple variables
variables_of_interest <- c("Platelets", "SCr", "Albumin", "Bilirubin", "Lymphocytes", "Neutrophils")
copula_long <- estimate_vinecopula_from_data(data_clean, polynomial = TRUE, ID_name = "ID", time_name = "gest",
                                             variables_of_interest = variables_of_interest, family_set = "parametric", keep_data = TRUE)
set.seed(68924730)
df_sim <- simulate(copula_long, n = length(unique(data_clean$ID)), value_only = FALSE)

plot(copula_long, var_names = "use", tree = 1:2)
pdf("results/figures/longitudinal_parameter_copula.pdf", width = 20, height = 20)
contour(copula_long)
pairs_copula_data(copula_long$uniform_data)
plot_comparison_distribution_sim_obs_generic(sim_data = df_sim$parameters, obs_data = copula_long$original_data, plot_type = "both")
dev.off()


#transform estimated parameters of observed population to time curves-> df_poly
gest_times <- seq(min(data_clean$gest), max(data_clean$gest), by = 1/7)
df_poly <- copula_long$original_data %>% 
  rownames_to_column("ID") %>%
  right_join(as.data.frame(expand_grid(gest = gest_times, ID = rownames(copula_long$original_data))))
for (variable in variables_of_interest) {
  col_ind <- grep(paste0("_", variable), names(df_poly))
  df_poly[, variable] <- df_poly[, col_ind[1]] + df_poly[, col_ind[2]]*df_poly$gest + df_poly[, col_ind[3]]*df_poly$gest^2
}

#Plotting results
#observed data
plot_explore <- data_clean %>% 
  select(all_of(c("ID", "gest", variables_of_interest))) %>% 
  pivot_longer(-c(ID, gest), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw()

#side-by-side comparison of the estimated and simulated curves
df_comparing_sim_org <- df_sim$values %>% 
  mutate(type = "Simulated") %>% 
  bind_rows(df_poly %>% mutate(type = "Observed")) %>% 
  dplyr::select(all_of(c("ID", "gest", variables_of_interest, "type"))) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey15", alpha = 0.5) +
  geom_smooth(color = "red", se = F) +
  labs(y = NULL, x = "Gestational age (weeks)") +
  facet_grid(variable ~ type, scales = "free_y", switch = "y") + 
  theme_bw() +
  theme(strip.placement = "outside")

pdf("results/figures/longitudinal_polynomials.pdf", width = 6, height = 9)
print(df_comparing_sim_org)
dev.off()

###
#performance on parameter scale
obs_statistics <- get_statistics(copula_long$original_data)
large_sim <-  simulate(copula_long, n = 5000, value_only = FALSE)$parameters


large_sim_statistics <- get_statistics(large_sim)

sim_error <- large_sim_statistics %>%   
  left_join(obs_statistics %>% rename(observed = value)) %>% 
  mutate(rel_error = (value - observed)/observed,
         ratio_diff = value/observed,
         abs_diff = value - observed) %>% 
  mutate(stat = ifelse(covariate == "all", "cov", statistic))


sim_error %>% 
  ggplot(aes(x = stat, y = ratio_diff)) +
  geom_boxplot() +
  geom_hline(yintercept = 1, color = "grey35") +  
  geom_hline(yintercept = c(1-0.2, 1+0.2), color = "grey35", linetype = 2) +  
  theme_bw()


sim_error %>% group_by(stat) %>% 
  summarize(median = median(ratio_diff), mean = mean(ratio_diff),
            median_error = 1 - median(ratio_diff), mean_error = 1 - mean(ratio_diff))



#check weird marginal profiles
marginal_long <- estimate_vinecopula_from_data(data_clean, polynomial = TRUE, ID_name = "ID", time_name = "gest",
                                             variables_of_interest = variables_of_interest, family_set = "indep")
set.seed(68924730)
df_sim_marg <- simulate(marginal_long, n = length(unique(data_clean$ID)), value_only = TRUE)

df_comparing_sim_org_marg <- df_sim_marg %>% 
  mutate(type = "Marginal") %>% 
  bind_rows(df_poly %>% mutate(type = "Observed")) %>% 
  bind_rows(df_sim$values %>% mutate(type = "Copula")) %>% 
  dplyr::select(all_of(c("ID", "gest", variables_of_interest, "type"))) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey15", alpha = 0.5) +
  geom_smooth(color = "red", se = F) +
  labs(y = NULL, x = "Gestational age (weeks)") +
  facet_grid(variable ~ type, scales = "free_y", switch = "y") + 
  theme_bw() +
  theme(strip.placement = "outside")



df_long_summary <- df_sim_marg %>% 
  mutate(type = "Marginal") %>% 
  bind_rows(df_poly %>% select(which(!grepl("b\\d_", names(df_poly)))) %>% mutate(type = "Observed")) %>% 
  bind_rows(df_sim$values %>% mutate(type = "Copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc") %>% 
  group_by(gest, type, biomarker) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = max(conc), min = min(conc)) %>% ungroup()

plot_df_long_summary <- df_long_summary %>% 
  ggplot(aes(x = gest, fill = type)) +
  geom_line(aes(y = median)) +
  geom_line(aes(y = p_25), linetype = 3) +
  geom_line(aes(y = p_75), linetype = 3) +
  geom_line(aes(y = p_high), linetype = 2) +
  geom_line(aes(y = p_low), linetype = 2) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  facet_grid(biomarker ~ type, scales = "free_y") +
  theme_classic()

plot_combination_long <-  df_sim_marg %>% 
  mutate(type = "Marginal") %>% 
  bind_rows(df_sim$values %>% mutate(type = "Copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc") %>% 
  ggplot(aes(x = gest, y = conc)) +
  geom_line(aes(group = ID), alpha = 0.3, color = "red") +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)") +
  
  geom_line(data = df_long_summary %>% filter(type == "Observed") %>% select(-type), aes(y = median)) +
  geom_line(data = df_long_summary %>% filter(type == "Observed") %>% select(-type), aes(y = p_25), linetype = 3) +
  geom_line(data = df_long_summary %>% filter(type == "Observed") %>% select(-type), aes(y = p_75), linetype = 3) +
  geom_line(data = df_long_summary %>% filter(type == "Observed") %>% select(-type), aes(y = p_high), linetype = 2) +
  geom_line(data = df_long_summary %>% filter(type == "Observed") %>% select(-type), aes(y = p_low), linetype = 2) +
  facet_grid(biomarker ~ type, scales = "free_y") +
  theme_classic()



AUC_df_long <- df_sim_marg %>% 
  mutate(type = "Marginal") %>% 
  bind_rows(df_poly %>% select(which(!grepl("b\\d_", names(df_poly)))) %>% mutate(type = "Observed")) %>% 
  bind_rows(df_sim$values %>% mutate(type = "Copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc") %>% 
  group_by(ID, type, biomarker) %>% 
  summarise(AUC = DescTools::AUC(x = gest, y = conc, 
                          method = "trapezoid", na.rm = FALSE, absolutearea = TRUE),
            AUC_late = DescTools::AUC(x = gest, y = conc, from = 32,
                                      method = "trapezoid", na.rm = FALSE, absolutearea = TRUE)) %>% 
  ungroup()


AUC_df_long %>% filter(type != "Observed") %>% 
  ggplot(aes(y = AUC, fill = type, x = biomarker)) +
  geom_boxplot(data = AUC_df_long %>% filter(type == "Observed"), alpha = 1) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("#FC6257", "#00AF2A", "grey65")) +
  facet_wrap(~ biomarker, scales = "free", nrow = 1) +
  theme_bw()


AUC_summary <- AUC_df_long %>% group_by(type, biomarker) %>% 
  summarize(sum_AUC = sum(AUC), sd_AUC = sd(AUC)) %>% ungroup()


AUC_summary %>% filter(type != "Observed") %>% 
  ggplot(aes(y = sum_AUC, color = type, x = biomarker)) +
  geom_boxplot(data = AUC_summary %>% filter(type == "Observed"), alpha = 1) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("#FC6257", "#00AF2A", "grey65")) +
  facet_wrap(~ biomarker, scales = "free", nrow = 1) +
  theme_bw()


AUC_summary %>% filter(type != "Observed") %>% 
  left_join(AUC_summary %>% filter(type == "Observed") %>% 
              rename(obs_sum = sum_AUC, obs_sd = sd_AUC) %>% select(-type)) %>% 
  mutate(ratio_diff = (sum_AUC - obs_sum)/obs_sum, ratio_diff_sd = sd_AUC/obs_sd)


AUC_df_long %>% group_by(type, biomarker) %>% 
  summarize(median_AUC = median(AUC), sd_AUC = sd(AUC)/median(AUC)) %>% ungroup() %>% filter(type != "Observed") %>% 
  ggplot(aes(y = sd_AUC, color = type, x = biomarker)) +
  geom_boxplot(data = AUC_df_long %>% group_by(type, biomarker) %>% 
                 summarize(median_AUC = median(AUC), sd_AUC = sd(AUC)/median(AUC)) %>% ungroup() %>% filter(type == "Observed"), alpha = 1) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("#FC6257", "#00AF2A", "grey65")) +
  facet_wrap(~ biomarker, scales = "free_x", nrow = 1) +
  theme_bw()



#larger simulations

df_sim_marg %>% 
  mutate(type = "Marginal") %>% 
  bind_rows(df_poly %>% select(which(!grepl("b\\d_", names(df_poly)))) %>% mutate(type = "Observed")) %>% 
  bind_rows(df_sim$values %>% mutate(type = "Copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc") %>% 
  group_by(gest, type, biomarker) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = max(conc), min = min(conc)) %>% ungroup()

plot_df_long_summary <- df_long_summary %>% 
  ggplot(aes(x = gest, fill = type)) +
  geom_line(aes(y = median)) +
  geom_line(aes(y = p_25), linetype = 3) +
  geom_line(aes(y = p_75), linetype = 3) +
  geom_line(aes(y = p_high), linetype = 2) +
  geom_line(aes(y = p_low), linetype = 2) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  facet_grid(biomarker ~ type, scales = "free_y") +
  theme_classic()




####### Summary metrics AUC bigger simulation #########
set.seed(68924730)
large_sim_cop <- simulate(copula_long, n = 5000, value_only = FALSE)
large_sim_marg <- simulate(marginal_long, n = 5000, value_only = FALSE)

#AUCs
AUC_df_long <- large_sim_marg$values %>% 
  mutate(type = "Marginal") %>% 
  bind_rows(df_poly %>% select(which(!grepl("b\\d_", names(df_poly)))) %>% mutate(type = "Observed")) %>% 
  bind_rows(large_sim_cop$values %>% mutate(type = "Copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc") %>% 
  group_by(ID, type, biomarker) %>% 
  summarise(AUC = DescTools::AUC(x = gest, y = conc, 
                                 method = "trapezoid", na.rm = FALSE, absolutearea = TRUE)) %>% 
  ungroup()


AUC_summary <- AUC_df_long %>% group_by(type, biomarker) %>% 
  #summarize(p_25 = quantile(AUC, 0.25), p_75 = quantile(AUC, 0.75)) %>% 
  #mutate(IQR = p_75 - p_25) %>% ungroup()
  summarize(IQR = sd(AUC)) %>% ungroup()

AUC_IQR <- AUC_summary %>% filter(type != "Observed") %>% 
  left_join(AUC_summary %>% filter(type == "Observed") %>% 
              rename(obs = IQR) %>% select(c(biomarker, obs))) %>% 
  mutate(ratio_diff = abs(IQR - obs)/IQR)

AUC_IQR %>% ggplot(aes(x = type, y = ratio_diff, color = biomarker, group = biomarker)) +
  geom_point() +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  labs(y = "Interquartile range(AUC)") +
  theme_bw()

AUC_IQR %>% group_by(type) %>% 
  summarize(mean_spread = mean(ratio_diff))

