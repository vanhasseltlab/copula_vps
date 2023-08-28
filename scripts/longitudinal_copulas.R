# Longitudinal covariate simulation
##### - Libraries - #####
library(tidyverse)
library(rvinecopulib)
library(kde1d)
library(lme4)
select <- dplyr::select

##### - Functions - #####
source("scripts/functions/functions.R")
source("scripts/functions/estimate_vinecopula_from_data.R")
color_palette <- create_colors(c("observed", "copula", "marginal\ndistribution", "conditional\ndistribution"), 
                               selected = c("grey", "turquoise", "dark yellow", "pink"))

#### Data ####
#data_reduced from data_preparation_pregnancy.R
data_reduced <- read.csv("data/clean/pregnancy_reduced.csv", row.names = NULL)
#remove all observations with baby (only looking at pregnancy)
data_clean <- data_reduced %>% 
  filter(BABY == 0)

##### - Simulations - #####
# Simulation settings
seed_nr <- 68924730


##### - Estimate polynomials and fit copula - ####
variables_of_interest <- c("Platelets", "SCr", "Albumin", "Bilirubin", "Lymphocytes", "Neutrophils")
copula_long <- estimate_vinecopula_from_data(data_clean, polynomial = TRUE, ID_name = "ID", time_name = "gest",
                                             variables_of_interest = variables_of_interest, family_set = "parametric", keep_data = TRUE)
set.seed(seed_nr)
df_sim <- simulate(copula_long, n = length(unique(data_clean$ID)), value_only = FALSE)


#### Marginal distributions ####
#check marginal profiles
marginal_long <- estimate_vinecopula_from_data(data_clean, polynomial = TRUE, ID_name = "ID", time_name = "gest",
                                               variables_of_interest = variables_of_interest, family_set = "indep")
set.seed(seed_nr)
df_sim_marg <- simulate(marginal_long, n = length(unique(data_clean$ID)), value_only = FALSE)

#### - Create observed curves from polynomial estimation - ####
# transform estimated parameters of observed population to time curves
gest_times <- seq(min(data_clean$gest), max(data_clean$gest), by = 1/7)
df_poly <- copula_long$original_data %>% 
  rownames_to_column("ID") %>%
  right_join(as.data.frame(expand_grid(gest = gest_times, ID = rownames(copula_long$original_data))))
for (variable in variables_of_interest) {
  col_ind <- grep(paste0("_", variable), names(df_poly))
  df_poly[, variable] <- df_poly[, col_ind[1]] + df_poly[, col_ind[2]]*df_poly$gest + df_poly[, col_ind[3]]*df_poly$gest^2
}

#### Combine results copula and marginal ####
df_long <- df_sim_marg$values %>% 
  mutate(type = "marginal\ndistribution") %>% 
  bind_rows(df_poly %>% select(which(!grepl("b\\d_", names(df_poly)))) %>% mutate(type = "observed")) %>% 
  bind_rows(df_sim$values %>% mutate(type = "copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc")

df_long_summary <- df_long %>% 
  group_by(gest, type, biomarker) %>% 
  summarize(median = median(conc), p_high = quantile(conc, 0.975), 
            p_low = quantile(conc, 0.025), p_25 = quantile(conc, 0.25), 
            p_75 = quantile(conc, 0.75), max = max(conc), min = min(conc)) %>% ungroup()


####Visualization####
# Plot data
long_all <- df_sim_marg$values %>% 
  mutate(type = "marginal\ndistribution") %>% 
  bind_rows(df_poly %>% select(which(!grepl("b\\d_", names(df_poly)))) %>% mutate(type = "observed")) %>% 
  bind_rows(df_sim$values %>% mutate(type = "copula")) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "biomarker", values_to = "conc")

names(color_palette)[1:3] <- c("observed", "copula", "marginal\ndistribution")

# Plot biomarkers over time
plot_df_long_summary_lines <- df_long_summary %>% 
  mutate(type = factor(type, levels = c("observed", "copula", "marginal\ndistribution"))) %>% 
  ggplot(aes(x = gest)) +
  geom_line(data = long_all %>% 
              mutate(type = factor(type, levels = c("observed", "copula", "marginal\ndistribution"))), 
            aes(y = conc, group = ID, color = type), alpha = 0.7, show.legend = F) +
  geom_line(aes(y = median, linetype = "Median")) +
  geom_line(aes(y = p_25, linetype = "Quartiles")) +
  geom_line(aes(y = p_75, linetype = "Quartiles")) +
  geom_line(aes(y = p_high, linetype = "95% quantiles")) +
  geom_line(aes(y = p_low, linetype = "95% quantiles")) +
  scale_linetype_manual(values = c(1, 3, 2), breaks = c("Median", "Quartiles", "95% quantiles"), name = NULL) + 
  scale_color_manual(values = color_palette, limits = force) +
  scale_x_continuous(name = "Time (weeks)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)", expand = expansion(mult = c(0.01, 0.05))) +
  facet_grid(biomarker ~ type, scales = "free_y")+
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))



pdf("results/figures/manuscript/F4_longitudinal_covariates.pdf", height = 8, width = 7)
print(plot_df_long_summary_lines)
dev.off()

#### Performance on parameter scale ####
perf_dat <- get_statistics(df_sim$parameters) %>% mutate(simulation = "copula") %>% 
  bind_rows(get_statistics(df_sim_marg$parameters) %>% mutate(simulation = "marginal distribution")) %>% 
  left_join(get_statistics(copula_long$original_data) %>% rename(observed = value)) %>% 
  mutate(rel_error = (value - observed)/observed,
         ratio_diff = value/observed,
         abs_diff = value - observed)

perf_dat %>%
  group_by(simulation, statistic) %>% 
  summarize(median = median(rel_error), mean = mean(rel_error))

#### Save copula object ####
copula_long$original_data <- NULL
copula_long$uniform_data <- NULL
copula_long$time_range <- c(0.3, 41)
save(copula_long, file = "copulas/longitudinal_copula.Rdata")
