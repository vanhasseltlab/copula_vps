#compare PK profiles Grimsley's model
library(DescTools)
library(tidyverse)
library(RxODE)
source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/functions.R")
source("scripts/functions/run_Grimsley.R")
source("scripts/functions/plot_distributions.R")
source("scripts/functions/Smania_Jonsson_MICE_simulation.R")
color_palette <- create_colors(c("Observed", "Copula", "Marginal", "CD", "MVN"), 
                               selected = c("grey", "turquoise", "dark yellow", "pink", "dark green"))

#cleaned data from simulation_comparison.R
data_grimsley <- read.csv("data/clean/pediatric_data.csv") %>%
  select(age, BW, CREA) %>% 
  mutate(BW = BW/1000) #convert to kilograms

set.seed(23856103+1)
copula_grimsley <- estimate_vinecopula_from_data(data_grimsley, family_set = "parametric")
sim_cop_grimsley <- simulate(copula_grimsley, n = nrow(data_grimsley))

#marginal
marg_grimsley <- estimate_vinecopula_from_data(data_grimsley, family_set ="indep")
sim_marg_grimsley <- simulate(marg_grimsley, n =  nrow(data_grimsley))

#cd
sim_cd_grimsley <- simCovMICE(m = 1, data_grimsley, catCovs = NULL) %>% select(-NSIM)

# run ODE model (solve for simulated patients)
df_true <- run_grimsley(n = nrow(data_grimsley), wgt = data_grimsley$BW, 
                        scr = data_grimsley$CREA, other_covariates = data_grimsley[, "age", drop = FALSE])
df_copula <- run_grimsley(n = nrow(sim_cop_grimsley), wgt = sim_cop_grimsley$BW, 
                          scr = sim_cop_grimsley$CREA, other_covariates = sim_cop_grimsley[, "age", drop = FALSE])
df_marg <- run_grimsley(n = nrow(sim_marg_grimsley), wgt = sim_marg_grimsley$BW,
                        scr = sim_marg_grimsley$CREA, other_covariates = sim_marg_grimsley[, "age", drop = FALSE])
df_cd <- run_grimsley(n = nrow(sim_cd_grimsley), wgt = sim_cd_grimsley$BW,
                        scr = sim_cd_grimsley$CREA, other_covariates = sim_cd_grimsley[, "age", drop = FALSE])
df_pk <- df_copula %>% mutate(type = "copula") %>% 
  bind_rows(df_true %>% mutate(type = "observed")) %>% 
  bind_rows(df_cd %>% mutate(type = "conditional")) %>% 
  bind_rows(df_marg %>% mutate(type = "marginal"))

# plot
plot_lines_weight <- df_pk %>% 
  filter(id %in% 1:nrow(data_grimsley)) %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = conc, group = id, color = wgt), alpha = 0.3) +
  scale_color_viridis_c(trans = "log10") +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  labs(color = "Weight (kg)") +
  facet_wrap(~ type) +
  theme_classic()


df_pk_summary <- df_pk %>% 
  group_by(time, type) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = quantile(conc, 0.99), min = quantile(conc, 0.01)) %>% ungroup()

#95% intervals
plot_pk_summary <- df_pk_summary %>% 
  ggplot(aes(x = time/60, fill = type)) +
  geom_line(aes(y = median)) +
  geom_line(aes(y = p_25), linetype = 3) +
  geom_line(aes(y = p_75), linetype = 3) +
  geom_line(aes(y = p_high), linetype = 2) +
  geom_line(aes(y = p_low), linetype = 2) +
  geom_ribbon(aes(ymin = min, ymax = max), alpha = 0.3) +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  facet_grid(~ type) +
  theme_classic()

pdf("results/figures/PK_model_assessment.pdf", width = 8, height = 5)

print(plot_pk_summary)
print(plot_lines_weight)

print(plot_lines_weight + geom_line(data = df_pk_summary, aes(y = median)) +
        geom_line(data = df_pk_summary, aes(y = p_25), linetype = 3) +
        geom_line(data = df_pk_summary, aes(y = p_75), linetype = 3) +
        geom_line(data = df_pk_summary, aes(y = p_high), linetype = 2) +
        geom_line(data = df_pk_summary, aes(y = p_low), linetype = 2))

dev.off()


df_pk_summary <- df_pk %>% 
  group_by(time, type) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = quantile(conc, 0.99), min = quantile(conc, 0.01)) %>% ungroup() %>% 
  mutate(Type = factor(str_to_title(type), levels = c("Observed", "Marginal", "Copula", "Conditional")))

plot_lines_weight_full <- df_pk %>% 
  mutate(Type = factor(str_to_title(type), levels = c("Observed", "Marginal", "Copula"))) %>% 
  filter(id %in% 1:nrow(data_grimsley)) %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = conc, group = id, color = wgt), alpha = 0.3) +
  #scale_color_viridis_c(trans = "log10", option = "plasma",direction = -1) +
  #scale_color_gradient(low = "#E4AB01", high = "#3ABAC1", trans = "log10") +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  labs(color = "Weight (kg)") +
  facet_wrap(~ Type) +
  theme_bw() + 
  geom_line(data = df_pk_summary, aes(y = median, linetype = "Median")) +
  geom_line(data = df_pk_summary, aes(y = p_25, linetype = "Quartiles")) +
  geom_line(data = df_pk_summary, aes(y = p_75, linetype = "Quartiles")) +
  geom_line(data = df_pk_summary, aes(y = p_high, linetype = "95% Quantiles")) +
  geom_line(data = df_pk_summary, aes(y = p_low, linetype = "95% Quantiles")) +
  scale_linetype_manual(values = c(1, 3, 2), breaks = c("Median", "Quartiles", "95% Quantiles"), name = NULL) + 
  theme(strip.background = element_rect(fill = "white"))


pdf("results/figures/manuscript/PK_model_assessment.pdf", width = 7, height = 4)
#print(plot_lines_weight_full + scale_color_viridis_c(trans = "log10", option = "plasma",direction = -1))
#print(plot_lines_weight_full + scale_color_gradient2(high = "#001158", mid ="#f46e32", low = "yellow", trans = "log10", midpoint = log10(0.3)))
print(plot_lines_weight_full + scale_color_gradient2(high = color_palette["MVN"], mid = color_palette["Copula"], low = color_palette["Marginal"], trans = "log10", midpoint = log10(3)))
print(plot_lines_weight_full + scale_color_gradient2(high = color_palette["Copula"], mid = color_palette["CD"], low = color_palette["Marginal"], trans = "log10", midpoint = log10(2)))

print(plot_lines_weight_full + scale_color_gradient2(high = color_palette["Copula"], mid = color_palette["CD"], low = color_palette["Marginal"], trans = "log10", midpoint = log10(1)))
print(plot_lines_weight_full + scale_color_gradient2(high = color_palette["MVN"], mid = color_palette["Marginal"], low = color_palette["CD"], trans = "log10", midpoint = log10(3)))

dev.off()


pdf("results/figures/PK_simulation_covariates.pdf", width = 10, height = 10)
plot_comparison_distribution_sim_obs_generic(sim_data = sim_cop_grimsley[1:nrow(data_grimsley), ], obs_data = data_grimsley, variables = names(sim_cop_grimsley),
                                             plot_type = "both", title = "Copula")
plot_comparison_distribution_sim_obs_generic(sim_data = sim_marg_grimsley[1:nrow(data_grimsley), ], obs_data = data_grimsley, variables = names(sim_marg_grimsley),
                                             plot_type = "both", title = "Marginal")


dev.off()



#Plot presentations
plot_data <- df_pk %>% 
  
  filter(time/60 <= 24) %>% 
  
  filter(id %in% 1:nrow(data_grimsley)) %>% 
  mutate(Type = factor(str_to_title(type), levels = c("Observed", "Marginal", "Copula", "Conditional"))) %>% 
  mutate(Type = recode(type, marginal = "marginal\ndistribution", conditional = "conditional\ndistribution"))

df_pk_summary_24 <-  df_pk_summary %>% 
  filter(time/60 <= 24) %>% 
  mutate(Type = recode(type, marginal = "marginal\ndistribution", conditional = "conditional\ndistribution"))

plot_lines_weight_pres <- plot_data %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = conc, group = id, color = wgt), alpha = 0.3)  +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  labs(color = "Weight (kg)") +
  facet_grid(~ Type) +
  theme_bw() + 
  geom_line(data = df_pk_summary_24, aes(y = median, linetype = "Median")) +
  geom_line(data = df_pk_summary_24, aes(y = p_25, linetype = "Quartiles")) +
  geom_line(data = df_pk_summary_24, aes(y = p_75, linetype = "Quartiles")) +
  geom_line(data = df_pk_summary_24, aes(y = p_high, linetype = "95% Quantiles")) +
  geom_line(data = df_pk_summary_24, aes(y = p_low, linetype = "95% Quantiles")) +
  scale_linetype_manual(values = c(1, 3, 2), breaks = c("Median", "Quartiles", "95% Quantiles"), name = NULL) + 
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 10))


pdf("presentation/figures/lines_PK_model_assessment.pdf", width = 7, height = 3.6)
print(plot_lines_weight_pres +
        scale_color_gradient2(low = "#f46e32", high = "#001158", trans = "log10", mid = "#243676", midpoint = log10(7)))
print(plot_lines_weight_pres +
        scale_color_gradient2(low = "#F54C00", mid = "#3652B7", trans = "log10", high = "#001158", midpoint = log10(5)))
dev.off()

#AUC comparison

AUC_df <- df_pk %>% 
  group_by(id, type, wgt, scr) %>% 
  summarise(AUC0_24 = AUC(x = time, y = conc, from = 0, to = 24*60, 
                          method = "trapezoid", na.rm = FALSE),
            AUC24_48 = AUC(x = time, y = conc, from = 24*60, to = 48*60, 
                           method = "trapezoid", na.rm = FALSE),
            AUC48_72 = AUC(x = time, y = conc, from = 48*60, to = 72*60, 
                           method = "trapezoid", na.rm = FALSE)) %>% 
  ungroup() %>% 
  pivot_longer(c(AUC0_24, AUC24_48, AUC48_72), names_to = "dosing_time", values_to = "AUC", names_prefix ="AUC")
AUC_df %>% 
  ggplot(aes(x = AUC, fill = type)) +
  geom_density(alpha = 0.3) +
  facet_grid(~ dosing_time) +
  theme_bw()

AUC_df %>% filter(dosing_time == "0_24") %>% group_by(type) %>% 
  summarize(cor_wgt_AUC = round(cor(wgt, AUC), 3),
            cor_scr_AUC = round(cor(scr, AUC), 3)) %>% as.data.frame

AUC_df %>% filter(dosing_time == "0_24") %>% 
  ggplot(aes(x = wgt, y = AUC)) +
  geom_point() +
  geom_smooth(se = F, color = "red") +
  facet_grid(~ type) +
  theme_bw()

AUC_df %>% filter(dosing_time == "0_24") %>% 
  ggplot(aes(x = scr, y = AUC)) +
  geom_point() +
  geom_smooth(se = F, color = "red") +
  facet_grid(~ type) +
  theme_bw()


AUC_df %>% filter(type != "observed") %>% 
  ggplot(aes(y = AUC, fill = type, x = dosing_time)) +
  geom_boxplot(data = AUC_df %>% filter(type == "observed"), alpha = 1) +
  geom_boxplot(alpha = 0.5) +
  scale_fill_manual(values = c("#FC6257", "#00AF2A", "grey65")) +
  
  theme_bw()

#repeat sampling for correlation AUC - BW and AUC - SCR
df_true <- run_grimsley(n = nrow(data_grimsley), wgt = data_grimsley$BW, 
                        scr = data_grimsley$CREA, other_covariates = data_grimsley[, "age", drop = FALSE])
AUC_observed <- df_true %>% 
  group_by(id, wgt, scr) %>% 
  summarise(AUC0_24 = AUC(x = time, y = conc, from = 0, to = 24*60, 
                          method = "trapezoid", na.rm = FALSE))

true_corr <- cor(AUC_observed[, -1])[3, 1:2]


copula_grimsley <- estimate_vinecopula_from_data(data_grimsley, family_set = "parametric")
marg_grimsley <- estimate_vinecopula_from_data(data_grimsley, family_set ="indep")
m <- 100
df_corr <- NULL

for (b in 1:m) {
  set.seed(9462056 + b)
  sim_cop_grimsley <- simulate(copula_grimsley, n = nrow(data_grimsley))
  #marginal
  sim_marg_grimsley <- simulate(marg_grimsley, n =  nrow(data_grimsley))
  
  # run ODE model (solve for simulated patients)
  df_copula <- run_grimsley(n = nrow(sim_cop_grimsley), wgt = sim_cop_grimsley$BW, 
                            scr = sim_cop_grimsley$CREA, other_covariates = sim_cop_grimsley[, "age", drop = FALSE])
  df_marg <- run_grimsley(n = nrow(sim_marg_grimsley), wgt = sim_marg_grimsley$BW,
                          scr = sim_marg_grimsley$CREA, other_covariates = sim_marg_grimsley[, "age", drop = FALSE])
  df_pk <- df_copula %>% mutate(type = "copula") %>% 
    bind_rows(df_marg %>% mutate(type = "marginal"))
  
  
  AUC_df <- df_pk %>% 
    group_by(id, type, wgt, scr) %>% 
    summarise(AUC0_24 = AUC(x = time, y = conc, from = 0, to = 24*60, 
                            method = "trapezoid", na.rm = FALSE)) %>% 
    ungroup %>% group_by(type) %>% 
    summarize(corr_weight = cor(wgt, AUC0_24),
              corr_scr = cor(scr, AUC0_24)) %>% ungroup() %>% as.data.frame()
  
  correlations <- data.frame(AUC_df)
  
  df_corr <- rbind.data.frame(df_corr, AUC_df)
}



#write.csv(df_corr, file = "results/PK_AUC_correlations.csv", row.names = FALSE)
df_corr <- read.csv("results/PK_AUC_correlations.csv")
df_cor_true <- data.frame(correlation = true_corr, covariate = c("weight", "scr"))

df_corr %>% 
  pivot_longer(-type, names_prefix = "corr_", names_to = "covariate", values_to = "correlation") %>% 
  ggplot(aes(y = correlation, fill = type, x = type)) +
  geom_hline(yintercept = 0, color = "grey65") +
  geom_boxplot() +
  geom_hline(data = df_cor_true, aes(yintercept = correlation), linetype = 2) +
  facet_wrap(~ covariate, scale = "free_x") +
  labs(x = NULL) +
  theme_bw()

df_corr %>% group_by(type) %>% 
  summarize(corr_weight = mean(corr_weight), corr_scr = mean(corr_scr))
  
