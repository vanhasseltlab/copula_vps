# PK profiles from simulated covariates (Grimsley)
##### - Libraries - #####
library(DescTools)
library(tidyverse)
library(RxODE)
library(patchwork)
library(kde1d)
library(rvinecopulib)
library(Rfast)

##### - Functions - #####
source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/functions.R")
source("scripts/functions/run_Grimsley.R")
source("scripts/functions/Smania_Jonsson_MICE_simulation.R")

##### - Data - #####
# Load data from simulation_comparison.R
data_grimsley <- read.csv("data/clean/pediatric_data.csv") %>%
  select(age, BW, CREA) %>% 
  mutate(BW = BW/1000) #convert to kilograms

##### - Simulations - #####
# Simulation settings
n_sim <- nrow(data_grimsley)
seed_nr <- 23856104

#### - Copulas: spline marginal + parametric copula - ####
set.seed(seed_nr)
copula_grimsley <- estimate_vinecopula_from_data(data_grimsley, family_set = "parametric")
sim_cop_grimsley <- simulate(copula_grimsley, n = n_sim)

#### - Marginals: splines - ####
set.seed(seed_nr)
marg_grimsley <- estimate_vinecopula_from_data(data_grimsley, family_set ="indep")
sim_marg_grimsley <- simulate(marg_grimsley, n =  n_sim)

#### - Multivariate (log-)Normal Distribution - ####
set.seed(seed_nr)
data_grimsley_log <- data_grimsley
for (j in 1:ncol(data_grimsley_log)) {
  data_grimsley_log[, j] <- log(data_grimsley_log[, j])
}
mvnorm_param <- mvnorm.mle(as.matrix(data_grimsley_log))
sim_mvnorm_grimsley <- as.data.frame(rmvnorm(n_sim, mu = mvnorm_param$mu, sigma = mvnorm_param$sigma, seed = seed_nr))
for (j in 1:ncol(sim_mvnorm_grimsley)) {
  sim_mvnorm_grimsley[, j] <- exp(sim_mvnorm_grimsley[, j])
}


##### - ODE model - #####
df_true <- run_grimsley(n = n_sim, wgt = data_grimsley$BW, 
                        scr = data_grimsley$CREA, other_covariates = data_grimsley[, "age", drop = FALSE])
df_copula <- run_grimsley(n = n_sim, wgt = sim_cop_grimsley$BW, 
                          scr = sim_cop_grimsley$CREA, other_covariates = sim_cop_grimsley[, "age", drop = FALSE])
df_marg <- run_grimsley(n = n_sim, wgt = sim_marg_grimsley$BW,
                        scr = sim_marg_grimsley$CREA, other_covariates = sim_marg_grimsley[, "age", drop = FALSE])
df_mvnorm <- run_grimsley(n = n_sim, wgt = sim_mvnorm_grimsley$BW,
                          scr = sim_mvnorm_grimsley$CREA, other_covariates = sim_mvnorm_grimsley[, "age", drop = FALSE])
##### - Visualize - #####
# Data prepration
df_pk <- df_copula %>% mutate(type = "copula") %>% 
  bind_rows(df_true %>% mutate(type = "observed")) %>% 
  bind_rows(df_marg %>% mutate(type = "marginal\ndistribution")) %>% 
  bind_rows(df_mvnorm %>% mutate(type = "mvnormal\ndistribution") %>% 
              mutate(conc = ifelse(conc < 0, 0, conc))) %>% 
  mutate(Type = factor(type, levels = c("observed", "copula", "mvnormal\ndistribution", "marginal\ndistribution")))


df_pk_summary <- df_pk %>% 
  group_by(time, Type) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = quantile(conc, 0.99), 
            min = quantile(conc, 0.01), mean = mean(conc)) %>% ungroup()

plot_data <- df_pk %>% 
  filter(time/60 <= 24) %>% 
  filter(id %in% 1:nrow(data_grimsley))

df_pk_summary_24 <- df_pk_summary %>% 
  filter(time/60 <= 24)

#save(df_pk_summary_24, plot_data, file = "results/PK_24h.Rdata")

# Plot PK curves
plot_lines_weight_pres <- plot_data %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = conc, group = id, color = wgt), alpha = 0.3)  +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  labs(color = "Weight (kg)", tag = "(a)") +
  facet_grid(~ Type) +
  theme_bw() + 
  geom_line(data = df_pk_summary_24, aes(y = median, linetype = "Median"), show.legend = T) +
  geom_line(data = df_pk_summary_24, aes(y = p_25, linetype = "Quartiles"), show.legend = T) +
  geom_line(data = df_pk_summary_24, aes(y = p_75, linetype = "Quartiles"), show.legend = T) +
  geom_line(data = df_pk_summary_24, aes(y = p_high, linetype = "95% quantiles"), show.legend = T) +
  geom_line(data = df_pk_summary_24, aes(y = p_low, linetype = "95% quantiles"), show.legend = T) +
  scale_linetype_manual(values = c(1, 3, 2), breaks = c("Median", "Quartiles", "95% quantiles"), name = NULL, labels = c("Median", "Quartiles", "95% quantiles")) + 
  scale_color_gradientn(colours = c("#f46e32", "#f46e32",  "#8592BC","#001158"), trans = "log10") +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 10))


# Plot AUC vs weight
AUC_df <- df_pk %>% 
  group_by(id, Type, wgt, scr) %>% 
  summarise(AUC0_24 = AUC(x = time, y = conc, from = 0, to = 24*60, 
                          method = "trapezoid", na.rm = FALSE),
            AUC24_48 = AUC(x = time, y = conc, from = 24*60, to = 48*60, 
                           method = "trapezoid", na.rm = FALSE),
            AUC48_72 = AUC(x = time, y = conc, from = 48*60, to = 72*60, 
                           method = "trapezoid", na.rm = FALSE)) %>% 
  ungroup() %>% 
  pivot_longer(c(AUC0_24, AUC24_48, AUC48_72), names_to = "dosing_time", values_to = "AUC", names_prefix ="AUC")

plot_auc_weight <- AUC_df %>% 
  filter(dosing_time == "0_24") %>% 
  ggplot(aes(x = wgt, y = AUC, color = Type)) +
  geom_point(alpha = 0.7, shape = 16, show.legend = F) +
  labs(x = "Body weight (kg)", tag = "(b)") +
  scale_color_manual(values = create_colors(c("copula", "marginal\ndistribution", "observed", "mvnormal\ndistribution"),
                                            c("turquoise", "dark yellow", "grey", "dark green")), limits = force) +
  facet_grid(~ Type) +
  scale_x_log10(breaks =  c(0.3, 1.0, 3.0, 10, 30), labels = c(0.3, 1.0, 3.0, 10, 30)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), strip.text = element_text(size = 10))


pdf(file = "results/figures/manuscript/F3_PK_and_weight_mvn.pdf", width = 8, height = 6.6)
plot_lines_weight_pres / plot_auc_weight + 
  plot_layout(heights = unit(c(2, 1), c('null', 'null')))
dev.off()

AUC_df %>% filter(dosing_time == "0_24") %>% 
  group_by(Type) %>% 
  summarize(correlation = cor(wgt, AUC)) %>% 
  ungroup()
