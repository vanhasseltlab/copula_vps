#mimic data
library(tidyverse)
load("data/MIMIC/data_wide_127vars.Rdata")
names(reduced_rows4)
str(reduced_rows4)
item_selection <- read.csv("data/MIMIC/labs cvh.csv")

item_selection[, "select"] <- ifelse(is.na(item_selection[, "select"] ), 0, item_selection[, "select"])
chosen_items <- setdiff(item_selection$label[item_selection$select == 1], c("Absolute Neutrophil Count", "Absolute Count - Retic"))
clean_MIMIC <- as.data.frame(apply(reduced_rows4[, c(3, 7:134)], 2, function(x) {
  x[x == "999999"] <- NA
  x <- as.numeric(x)
  return(x)})) %>%
  select(all_of(c(chosen_items, "weight", "anchor_age"))) 


cop_data_MIMIC <- clean_MIMIC %>% 
  filter(rowSums(!is.na(clean_MIMIC)) > 1)

rm(clean_MIMIC)
rm(reduced_rows4)

source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/functions.R")
source("scripts/functions/plot_distributions.R")
####!RUN ON SERVER!###
start_time <- Sys.time()
cop_mimic <- estimate_vinecopula_from_data(dat = cop_data_MIMIC, 
                                           variables_of_interest = c(chosen_items, "weight", "anchor_age"), 
                                           keep_data = FALSE, family_set = "parametric")

#save(cop_mimic, file = "copula_mimic_12var_5kpatients.Rdata")
end_time <- Sys.time()
print(end_time - start_time)


###!END RUN ON SERVER!###

load("results/MIMIC/copula_mimic_allvar_53kpatients_filtered.Rdata")
cop_mimic$variables_of_interest

set.seed(123)
mimic_simulated <- simulate(cop_mimic, 1000)

pdf("results/MIMIC/plot_contour_32vars.pdf", width = 50, height = 50)
#contour(cop_mimic)
#plot(cop_mimic, tree = 1:2, var_names = "use")
plot_comparison_distribution_sim_obs_generic(sim_data = mimic_simulated, obs_data = cop_data_MIMIC, plot_type = "density", 
                                             variables = cop_mimic$variables_of_interest[1:10])
dev.off()

cors <- cor(cop_data_MIMIC, use = "pairwise.complete.obs")
cors_df <- data.frame(row=rownames(cors)[row(cors)[upper.tri(cors)]], 
           col=colnames(cors)[col(cors)[upper.tri(cors)]], 
           corr=cors[upper.tri(cors)])
vars_corr <- c("Albumin", "C Reactive Protein (CRP)", "Absolute Neutrophil Count", "HDL", "INR", "Prothrombin time", "Hemoglobin", "Hematocrit (serum)", "ALT", "AST")


pdf("results/MIMIC/plot_contour_32vars.pdf", width = 50, height = 50)
#contour(cop_mimic)
#plot(cop_mimic, tree = 1:2, var_names = "use")
plot_comparison_distribution_sim_obs_generic(sim_data = mimic_simulated, obs_data = cop_data_MIMIC, plot_type = "density", 
                                             variables = vars_corr)
dev.off()

nr_pairwise_obs <- psych::pairwiseCount(cop_data_MIMIC)



set.seed(123)
mimic_simulated_large <- simulate(cop_mimic, 10000)
stats_copula <- get_statistics(mimic_simulated_large)
stats_obs <- get_statistics(cop_data_MIMIC)


stats_both <- stats_copula %>% left_join(stats_obs %>% rename(observed = value)) %>% 
  mutate(rel_error = (value - observed)/observed,
         ratio_diff = value/observed,
         abs_diff = value - observed)


stats_both %>% filter(statistic == "covariance") %>% 
  summarize(1 - median(ratio_diff, na.rm =  T), median(rel_error, na.rm =  T))

stats_both %>% filter(statistic == "covariance") %>% 
  ggplot(aes(y = rel_error)) +
  #geom_point() +
  geom_boxplot(outlier.shape = NA)  +
  scale_y_continuous(limits = quantile(stats_both$rel_error, c(0.1, 0.9), na.rm = TRUE)) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  theme_bw()

stats_both %>% filter(statistic == "covariance") %>% 
  ggplot(aes(y = ratio_diff)) +
  geom_boxplot(outlier.shape = NA)  +
  scale_y_continuous(limits = quantile(stats_both$ratio_diff, c(0.1, 0.9), na.rm = TRUE)) +
  geom_hline(yintercept = c(1-0.2, 1+0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 1, linetype = 2, color = "black") +
  theme_bw()


stats_both %>% filter(statistic == "covariance") %>% 
  ggplot(aes(y = rel_error, x = covariate)) +
  geom_point() +
  scale_y_continuous(limits = quantile(stats_both$rel_error, c(0.05, 0.95), na.rm = TRUE)) +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 2, color = "black") +
  theme_bw()


dens_plot <- cop_data_MIMIC %>% 
  filter(INR < 5) %>% 
  ggplot(mapping = aes(y = `Prothrombin time`, x = INR)) +
  geom_density2d(alpha = 1, linetype = 2) +
  geom_point(color = "red", alpha = 0.1) +
  theme_bw() +
  theme(legend.position = "none")


str(cop_data_MIMIC)

#Copula on discrete data

check_integer <- function(x) {
  all(x%%1 == 0, na.rm = TRUE)
}

discrete_variables <- sapply(cop_data_MIMIC, check_integer)

#load results from server
load("results/MIMIC/copula_mimic_allvar_53kpatients_filtered_disc.Rdata")

set.seed(234973580)
mimic_simulated_large <- simulate(cop_mimic_disc, 10000)

stats_copula <- get_statistics(mimic_simulated_large)
stats_obs <- get_statistics(cop_data_MIMIC)


stats_both <- stats_copula %>% left_join(stats_obs %>% rename(observed = value)) %>% 
  mutate(rel_error = (value - observed)/observed,
         ratio_diff = value/observed,
         abs_diff = value - observed)
