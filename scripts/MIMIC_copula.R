#mimic data
library(tidyverse)
library(ggside)
library(patchwork)


source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/functions.R")
source("scripts/functions/plot_distributions.R")


#Load and clean data
load("data/MIMIC/data_wide_127vars.Rdata")
item_selection <- read.csv("data/MIMIC/labs cvh.csv")
item_selection[, "select"] <- ifelse(is.na(item_selection[, "select"] ), 0, item_selection[, "select"])

#Remove unused items
chosen_items <- setdiff(item_selection$label[item_selection$select == 1], c("Absolute Neutrophil Count", "Absolute Count - Retic"))
clean_MIMIC <- as.data.frame(apply(reduced_rows4[, c(3, 7:134)], 2, function(x) {
  x[x == "999999"] <- NA
  x <- as.numeric(x)
  return(x)})) %>%
  select(all_of(c(chosen_items, "weight", "anchor_age"))) 

#filter and transform variables
cop_data_MIMIC <- clean_MIMIC %>% 
  filter(rowSums(!is.na(clean_MIMIC)) > 1) %>% 
  mutate(logAST = ifelse(AST == 0, NA, log(AST))) %>% 
  mutate(logBNP = ifelse(`Brain Natiuretic Peptide (BNP)` == 0, NA, log(`Brain Natiuretic Peptide (BNP)`))) %>% 
  select(-c(AST, `Brain Natiuretic Peptide (BNP)`)) %>% 
  mutate(weight = ifelse(weight > 500, weight/10, weight))

#remove large objects
rm(clean_MIMIC)
rm(reduced_rows4)

#save cleaned data
save(cop_data_MIMIC, file = "data/clean/data_copula_mimic_clean.Rdata")


#Fit copula to MIMIC data
####!RUN ON SERVER!###
start_time <- Sys.time()
cop_mimic <- estimate_vinecopula_from_data(dat = cop_data_MIMIC, 
                                           variables_of_interest = names(cop_data_MIMIC), 
                                           keep_data = FALSE, family_set = "parametric", cores = 40)
save(cop_mimic, file = "results/MIMIC/copula_mimic_allvar_53kpatients_filtered_logtransformed.Rdata")
end_time <- Sys.time()
print(end_time - start_time)

###!END RUN ON SERVER!###


#explore copula
load("results/MIMIC/copula_mimic_allvar_53kpatients_filtered_logtransformed.Rdata") #cop_mimic
load("data/clean/data_copula_mimic_clean.Rdata") #cop_data_MIMIC

cat(paste(cop_mimic$vine_copula$names, collapse = ","))

set.seed(83520583)
large_sim <- simulate(cop_mimic, 5000)

plot(cop_mimic$marginals$Albumin$dist$grid_points, cop_mimic$marginals$Albumin$dist$values)

cop_mimic$marginals$Albumin$dist


#test without the data
try_cov <- cop_data_MIMIC$Albumin[!is.na(cop_data_MIMIC$Albumin)]
try_kernel <- kde1d(try_cov)

try_kernel$x <- NULL

rkde1d(10, try_kernel)

try_cop <- large_cop$vine_copula$pair_copulas[[1]][[1]]

large_cop

append(list(1, 5, 3), list(4), after = 2)

pdf("results/figures/all_densities_clean_log.pdf", width = 60, height = 60)
plot_comparison_distribution_sim_obs_generic(sim_data = large_sim, obs_data = cop_data_MIMIC, pick_color = c("#3ABAC1", "#969696"),
                                             plot_type = "density", variables = names(large_sim))
dev.off()

##
#visualize 6 examples
all_MIMIC <- large_sim %>% mutate(type = "copula") %>% 
  bind_rows(cop_data_MIMIC %>% mutate(type = "observed")) %>% 
  rename(Age = anchor_age, Weight = weight)

sets_of_interest <- matrix(c(c("BUN", "`C Reactive Protein (CRP)`", "BUN", "C Reactive Protein (CRP)"),
                             c("HDL", "Triglyceride", "HDL", "Triglyceride"),
                             c("`C Reactive Protein (CRP)`", "`Glucose (serum)`", "C Reactive Protein (CRP)", "Glucose (serum)"),
                             c("Albumin", "`C Reactive Protein (CRP)`", "Albumin", "C Reactive Protein (CRP)"),
                             c("Weight", "Age", "Weight", "Age"),
                             c("`Hematocrit (serum)`", "Hemoglobin", "Hematocrit (serum)", "Hemoglobin")), 
                           ncol = 4, byrow = T)

MIMIC_plot_list <- list()
for (i in 1:nrow(sets_of_interest)) {
  
  plot_dat_ind <- !is.na(all_MIMIC[, sets_of_interest[i, 3]]) & !is.na(all_MIMIC[, sets_of_interest[i, 4]])
  MIMIC_plot_list[[i]] <- all_MIMIC[plot_dat_ind, ] %>% 
    ggplot(aes_string(x = sets_of_interest[i, 1], y = sets_of_interest[i, 2], color = "type", linetype = "type")) +
    geom_density_2d(bins = 10, show.legend = F) +
    geom_point(shape = NA, show.legend = F) +
    scale_linetype_manual(values = c(1, 5)) +
    scale_color_manual(values =  create_colors(c("observed", "copula", "marginal distribution", "conditional distribution", "mvnormal distribution", "bootstrap"),
                                               selected = c("grey", "turquoise", "dark yellow", "pink", "dark green", "midlight blue")), limits = force) +
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       limits = quantile(all_MIMIC[, sets_of_interest[i, 3]], probs = c(0.01, 0.95), na.rm = TRUE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0)),
                       limits = quantile(all_MIMIC[, sets_of_interest[i, 4]], probs = c(0.01, 0.95), na.rm = TRUE)) +
    geom_xsidedensity(show.legend = F) +
    geom_ysidedensity(show.legend = F) +
    scale_ysidex_continuous(minor_breaks = NULL, limits = c(0,NA), expand = expansion(mult = c(0, 0.05))) +
    scale_xsidey_continuous(minor_breaks = NULL, limits = c(0,NA), expand = expansion(mult = c(0, 0.05))) +
    theme_bw() +
    theme(aspect.ratio = 1, ggside.panel.grid = element_blank(), ggside.axis.line = element_line(color = "white"),
          ggside.axis.text = element_blank(), ggside.axis.ticks = element_blank(), 
          ggside.panel.border = element_rect(colour = "white"), ggside.panel.scale = .2)
}

pdf(file = "results/MIMIC/selected_densities_patch.pdf", width = 0.85*9, height = 0.85*6)
(MIMIC_plot_list[[1]] + MIMIC_plot_list[[2]] + MIMIC_plot_list[[3]]) / (MIMIC_plot_list[[4]] + MIMIC_plot_list[[5]] + MIMIC_plot_list[[6]])
dev.off()




###################
cors <- cor(cop_data_MIMIC, use = "pairwise.complete.obs")
cors_df <- data.frame(row=rownames(cors)[row(cors)[upper.tri(cors)]], 
           col=colnames(cors)[col(cors)[upper.tri(cors)]], 
           corr=cors[upper.tri(cors)])
vars_corr <- c("Albumin", "C Reactive Protein (CRP)", "HDL", "INR", "Prothrombin time", "Hemoglobin", "Hematocrit (serum)", "ALT", "AST")


pdf("results/MIMIC/plot_contour_32vars_v2.pdf", width = 20, height = 20)
#contour(cop_mimic)
#plot(cop_mimic, tree = 1:2, var_names = "use")
plot_comparison_distribution_sim_obs_generic(sim_data = mimic_simulated, obs_data = cop_data_MIMIC, plot_type = "density", 
                                             variables = vars_corr)
dev.off()

nr_pairwise_obs <- psych::pairwiseCount(cop_data_MIMIC)



stats_copula <- get_statistics(large_sim)
stats_obs <- get_statistics(cop_data_MIMIC)


stats_both <- stats_copula %>% left_join(stats_obs %>% rename(observed = value)) %>% 
  mutate(covariate = gsub("anchor_age", "age", covariate, fixed = TRUE)) %>% 
  mutate(rel_error = (value - observed)/observed,
         ratio_diff = value/observed,
         abs_diff = value - observed) %>% 
  mutate(cov1 = gsub("\\_.*", "", covariate),
         cov2 = gsub(".*\\_", "", covariate))
stats_both[, c("covA", "covB")] <- t(apply(as.data.frame(stats_both[, c("cov1", "cov2")]), 1, sort))

vars_off <- c("Brain Natiuretic Peptide (BNP)", "Total Protein", "Absolute Count - Basos", "HDL", "C Reactive Protein (CRP)", "Potassium (serum)", "Hematocrit (serum)")


pdf("results/MIMIC/plot_contour_32vars_off.pdf", width = 15, height = 15)
#contour(cop_mimic)
#plot(cop_mimic, tree = 1:2, var_names = "use")
plot_comparison_distribution_sim_obs_generic(sim_data = mimic_simulated, obs_data = cop_data_MIMIC, plot_type = "density", 
                                             variables = vars_off)
dev.off()


stats_both %>% filter(statistic == "covariance") %>% 
  summarize(1 - median(ratio_diff, na.rm =  T), median(rel_error, na.rm =  T), 
            RMSE = sqrt(mean(rel_error^2)))

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



plot_SJ_30d <- stats_both %>% 
  filter(statistic == "covariance") %>% 
  ggplot(aes(y = covA, x = covB)) +
  
  geom_tile(aes(fill = abs(rel_error))) +
  geom_hline(yintercept = seq(0.5, 11.5, by = 1), color = "white") +
  geom_vline(xintercept = seq(0.5, 11.5, by = 1), color = "white") +
  geom_text(aes(label = round(rel_error, 2))) +
  scale_fill_gradientn(colours = c("white", "#f46e32"), limits = c(0, 1)) +
  scale_x_discrete(expand = expansion(mult = c(0, 0)), guide = guide_axis(angle = 90)) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"), panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(), panel.grid = element_blank())


stats_both %>% 
  filter(abs(rel_error) < 13) %>% 
  filter(statistic == "covariance") %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE),
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  ggplot(aes(y = rel_error, x = covariate)) +
  geom_hline(yintercept = 0, color = "grey35") +  
  geom_hline(yintercept = c(-0.2, 0.2), color = "grey35", linetype = 2) +  
  geom_point(color = "#3ABAC1") +
  labs(x = "Covariate combinations", y = "Relative error", color = "Method", shape = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 6),  strip.background = element_rect(fill = "white"),
        strip.text = element_text( margin = margin( b = 1, t = 1), size = 12))


stats_both %>% 
  filter(abs(rel_error) < 13) %>% 
  filter(statistic == "covariance") %>% 
  mutate(statistic = gsub("sd", "standard deviation", statistic, fixed = TRUE),
         covariate = gsub("_", " - ", covariate, fixed = TRUE)) %>% 
  mutate(covariate = gsub("CREA", "SCr", covariate, fixed = TRUE)) %>% 
  ggplot(aes(y = rel_error, x = observed)) +
  geom_point(color = "#3ABAC1") +
  #labs(x = "Covariate combinations", y = "Relative error", color = "Method", shape = "Method") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text( margin = margin( b = 1, t = 1), size = 12))
