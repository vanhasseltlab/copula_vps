#MIMIC data analysis -> run on machine with large memory!
cat("start script \n")
##### - Libraries - #####
library(tidyverse)
library(ggside)
library(patchwork)
library(kde1d)
library(rvinecopulib)

##### - Functions - #####
source("scripts/functions/estimate_vinecopula_from_data.R")
source("scripts/functions/functions.R")
source("scripts/functions/plot_distributions.R")


##### - Data - #####
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
save(cop_mimic, file = "copulas/mimic_copula.Rdata")
end_time <- Sys.time()
cat("time of copula fitting: ", end_time - start_time, "\n")
###!END RUN ON SERVER!###


#explore copula
load("copulas/mimic_copula.Rdata") #cop_mimic
load("data/clean/data_copula_mimic_clean.Rdata") #cop_data_MIMIC

#### - Simulation - ####
set.seed(83520583)
large_sim <- simulate(cop_mimic, 5000)

cat("Visualize all densities \n")
plot_distributions_mimic <- plot_comparison_distribution_sim_obs_generic(sim_data = large_sim, obs_data = cop_data_MIMIC, pick_color = c("#3ABAC1", "#969696"),
                                             plot_type = "density", variables = names(large_sim), grob = TRUE, 
                                             caption = "\n  Figure S2: densities of all covariates in MIMIC which were modeled using a copula. The density of the simulations from the copula (blue solid line) and the \n  observed density (grey dashed line) show overlap for many covariate combinations.")

ggsave(file = "results/figures/manuscript/FS2_full_copula_density_MIMIC.pdf", plot_distributions_mimic, width = 12, height = 12, units = "in", scale = 70/12, limitsize = FALSE)


##### - Visualization - #####
# visualize 6 examples
cat("Visualize set of 6 densities\n")
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

pdf(file = "results/figures/manuscript/F5_six_densities_MIMIC.pdf", width = 0.85*9, height = 0.85*6)
(MIMIC_plot_list[[1]] + MIMIC_plot_list[[2]] + MIMIC_plot_list[[3]]) / (MIMIC_plot_list[[4]] + MIMIC_plot_list[[5]] + MIMIC_plot_list[[6]])
dev.off()


# calculate performance
cop_sim <- get_statistics(large_sim)
obs <- get_statistics(cop_data_MIMIC)

performance_data <- cop_sim %>%
  filter(statistic == "correlation" ) %>% 
  mutate(key = paste(statistic, covariate, sep = "_")) %>% 
  left_join(obs %>% dplyr::rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed)

median(performance_data$rel_value)
