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
source("scripts/functions/estimate_vinecopula_from_data.R")

variables_of_interest <- c("Platelets", "SCr", "Albumin", "Bilirubin", "Lymphocytes", "Neutrophils")
copula_long <- estimate_vinecopula_from_data(data_clean, polynomial = TRUE, ID_name = "ID", time_name = "gest",
                                             variables_of_interest = variables_of_interest, family_set = "parametric")
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


# create noisy subjects
## create copula for gest and residual noise
variables_poly <- paste0("b", 0:2, "_", rep(variables_of_interest, each = 3))
df_clean_pred <- data_clean %>% 
  left_join(df_poly %>% select(all_of(c("ID", variables_poly))) %>% 
              mutate(ID = as.integer(ID)) %>% distinct())

for (variable in variables_of_interest) {
  col_ind <- grep(paste0("_", variable), names(df_clean_pred))
  df_clean_pred[, paste0(variable, "_pred")] <- df_clean_pred[, col_ind[1]] + df_clean_pred[, col_ind[2]]*df_clean_pred$gest + df_clean_pred[, col_ind[3]]*df_clean_pred$gest^2
}

