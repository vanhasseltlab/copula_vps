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

df_pk_summary <- df_pk %>% 
  group_by(time, type) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = quantile(conc, 0.99), 
            min = quantile(conc, 0.01), mean = mean(conc)) %>% ungroup() %>% 
  mutate(Type = factor(str_to_title(type), levels = c("Observed", "Marginal", "Copula", "Conditional")))

# plot
plot_lines_weight_full <- df_pk %>% 
  mutate(Type = factor(str_to_title(type), levels = c("Observed", "Marginal", "Copula", "Conditional"))) %>% 
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
  mutate(Type = recode(type, marginal = "marginal\ndistribution", conditional = "conditional\ndistribution")) %>% 
  mutate(Type = factor(Type, levels = c("observed", "copula", "marginal\ndistribution", "conditional\ndistribution"))) %>% 
  filter(Type != "conditional\ndistribution")

df_pk_summary_24 <-  df_pk_summary %>% 
  filter(time/60 <= 24) %>% 
  mutate(Type = recode(type, marginal = "marginal\ndistribution", conditional = "conditional\ndistribution")) %>% 
  mutate(Type = factor(Type, levels = c("observed", "copula", "marginal\ndistribution", "conditional\ndistribution"))) %>% 
  filter(Type != "conditional\ndistribution")

plot_lines_weight_pres <- plot_data %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = conc, group = id, color = wgt), alpha = 0.3)  +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  labs(color = "Weight (kg)") +
  facet_grid(~ Type) +
  theme_bw() + 
  geom_line(data = df_pk_summary_24, aes(y = median, linetype = "Median"), show.legend = F) +
  geom_line(data = df_pk_summary_24, aes(y = p_25, linetype = "Quartiles"), show.legend = F) +
  geom_line(data = df_pk_summary_24, aes(y = p_75, linetype = "Quartiles"), show.legend = F) +
  geom_line(data = df_pk_summary_24, aes(y = p_high, linetype = "95% Quantiles"), show.legend = F) +
  geom_line(data = df_pk_summary_24, aes(y = p_low, linetype = "95% Quantiles"), show.legend = F) +
  scale_linetype_manual(values = c(1, 3, 2), breaks = c("Median", "Quartiles", "95% Quantiles"), name = NULL) + 
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 10))


pdf("presentation/figures/lines_PK_model_assessment.pdf", width = 7, height = 3.6)
print(plot_lines_weight_pres +
        scale_color_gradientn(colours = c("#f46e32", "#f46e32", "#8592BC", "#001158"), trans = "log10"))
dev.off()

pdf("presentation/figures/lines_PK_model_assessment.pdf", width = 6, height = 3.6)
print(plot_lines_weight_pres +
        scale_color_gradientn(colours = c("#f46e32", "#f46e32", "#8592BC", "#001158"), trans = "log10"))
dev.off()


#Typical PK
plot_typical_PK <- df_pk_summary_24 %>% 
  filter(type == "observed") %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = median), color = "#8592BC", size = 1.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time (hour)", y = "Vancomycin concentration (mg/L)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(file = "presentation/figures/line_typical_pk_Grimsey.pdf", width = 2.7, height = 2.7)
print(plot_typical_PK)
dev.off()

plot_one_PK <- df_pk_summary_24 %>% 
  filter(type == "observed") %>% 
  ggplot(aes(x = time/60)) +
  geom_ribbon(aes(ymin = p_25,ymax = p_75), fill = "grey80") +
  geom_line(aes(y = median), color = "grey36", size = 1.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24), limits = c(0, 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time", y = "Concentration") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

plot_data_PK <- plot_data %>% 
  filter(time %in% c(60, 120, 600, 1400)) %>% 
  filter(id %in% 1:10) %>% 
  filter(Type == "copula") %>% 
  ggplot(aes(y = conc, x = time/60)) +
  stat_summary(fun = mean, geom = "line", color = "grey49", size = 1.5, linetype = 1) +
  geom_point()  +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24), limits = c(0, 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time", y = "Concentration") +
  labs(color = "Weight (kg)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


col_pats <- read.csv("example/simulated_data_conceptual_figure.csv") %>% 
  arrange(size) %>% 
  mutate(id = 1:20)

plot_model_PK_variation <- plot_data %>% 
  filter(id %in% 1:20) %>% 
  left_join(col_pats) %>% 
  filter(Type == "copula") %>% 
  ggplot(aes(y = conc, x = time/60, color = color_hexa2)) +
  geom_ribbon(data = df_pk_summary_24 %>% filter(type == "copula"), 
              aes(ymin = p_low, ymax = p_high, x = time/60), fill = "grey80", inherit.aes = F) +
  geom_line(aes(group = id),size = 1.5)  +
  scale_color_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24), limits = c(0, 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time", y = "Concentration") +
  labs(color = "Weight (kg)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

pdf(file = "presentation/figures/line_pks_conceptual.pdf", width = 3.57, height = 2)
print(plot_one_PK + labs(x = NULL, y = NULL))
print(plot_data_PK + labs(x = NULL, y = NULL))
print(plot_model_PK_variation + labs(x = NULL, y = NULL))

print(plot_one_PK)
print(plot_data_PK)
print(plot_model_PK_variation)
dev.off()




#######
#AUC vs weight pointplots
AUC_df <- df_pk %>% 
  group_by(id, type, wgt, scr) %>% 
  summarise(AUC0_24 = AUC(x = time, y = conc, from = 0, to = 24*60, 
                          method = "trapezoid", na.rm = FALSE),
            AUC24_48 = AUC(x = time, y = conc, from = 24*60, to = 48*60, 
                           method = "trapezoid", na.rm = FALSE),
            AUC48_72 = AUC(x = time, y = conc, from = 48*60, to = 72*60, 
                           method = "trapezoid", na.rm = FALSE)) %>% 
  ungroup() %>% 
  pivot_longer(c(AUC0_24, AUC24_48, AUC48_72), names_to = "dosing_time", values_to = "AUC", names_prefix ="AUC") %>% 
  mutate(Type = recode(type, marginal = "marginal\ndistribution", conditional = "conditional\ndistribution")) %>% 
  mutate(Type = factor(Type, levels = c("observed", "copula", "marginal\ndistribution", "conditional\ndistribution")))


plot_auc_weight <- AUC_df %>% 
  filter(dosing_time == "0_24") %>% 
  filter(type != "conditional") %>% 
  ggplot(aes(x = wgt, y = AUC, color = Type)) +
  geom_point(alpha = 0.7, shape = 16) +
  labs(x = "Body weight (kg)") +
  scale_color_manual(values = create_colors(c("copula", "marginal\ndistribution", "observed", "conditional\ndistribution"),
                                            c("turquoise", "dark yellow", "grey", "pink")), limits = force) +
  facet_grid(~ Type) +
  theme_bw()

pdf(file = "results/figures/PK_AUC_weight_scatter.pdf", width = 7, height = 3)
print(plot_auc_weight)
dev.off()
###
# geom_ribbon(data=subset(x, 2 <= x & x <= 3), 
#             aes(ymin=twox,ymax=x2), fill="blue", alpha=0.5) 
try_new_plot <- df_pk_summary_24 %>% 
  filter(Type != "observed") %>% 
  ggplot(aes(x = time/60, color = Type)) +
  geom_line(aes(y = median, linetype = "Median"), show.legend = F) +
 
  geom_ribbon(data = df_pk_summary_24 %>% 
                filter(Type == "observed"), aes(ymin = p_low, ymax = p_high), fill = "grey95", color = "white") +
  geom_ribbon(data = df_pk_summary_24 %>% 
                filter(Type == "observed"), aes(ymin = p_25, ymax = p_75), fill = "grey85", color = "white") +
  geom_line(data = df_pk_summary_24 %>% 
              filter(Type == "observed"), aes(y = median), color = "black") +
  geom_line(aes(y = p_25, linetype = "Quartiles"), show.legend = F) +
  geom_line(aes(y = p_75, linetype = "Quartiles"), show.legend = F) +
  geom_line(aes(y = p_high, linetype = "95% Quantiles"), show.legend = F) +
  geom_line(aes(y = p_low, linetype = "95% Quantiles"), show.legend = F) +
  geom_line(aes(y = median, linetype = "Median"), show.legend = F) +
  scale_linetype_manual(values = c(1, 2, 3), breaks = c("Median", "Quartiles", "95% Quantiles"), name = NULL) + 
  scale_color_manual(values = create_colors(c("copula", "marginal\ndistribution", "observed"),
                                            c("turquoise", "dark yellow", "grey"))) +
  theme_bw() 

pdf("results/figures/PK_shaded_area.pdf", width = 4, height = 4)
print(try_new_plot)
dev.off()

#AUC comparison


AUC_df %>% 
  ggplot(aes(x = AUC, fill = type)) +
  geom_density(alpha = 0.3) +
  facet_grid(~ dosing_time) +
  theme_bw()

AUC_df %>% filter(dosing_time == "0_24") %>% group_by(type) %>% 
  summarize(cor_wgt_AUC = round(cor(wgt, AUC), 3),
            cor_scr_AUC = round(cor(scr, AUC), 3)) %>% as.data.frame

AUC_df %>% filter(dosing_time == "0_24") %>% 
  ggplot(aes(x = wgt, y = AUC, color = type)) +
  geom_point() +
  scale_color_manual(values = create_colors(c("copula", "marginal", "observed", "conditional"),
                                            c("turquoise", "dark yellow", "grey", "pink"))) +
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
  scale_fill_manual(values = c("#FC6257", "#00AF2A",  "blue", "grey65")) +
  
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
  summarize(corr_weight = mean(corr_weight), corr_scr = mean(corr_scr)) %>% View
  
