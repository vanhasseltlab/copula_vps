library(tidyverse)
library(data.table)

#unzip MIMIC files (need to do only once)
mimic_path <- "data/MIMIC/mimic-iv-1.0/"
for (filename in list.files(mimic_path, recursive = TRUE, pattern = "*\\.gz$")) {
  R.utils::gunzip(paste0(mimic_path, filename), remove = TRUE)
}


icu.d_items <- fread("data/MIMIC/mimic-iv-1.0/icu/d_items.csv")
id_serum_creatinine <- hops.d_labitems$itemid[(grepl("Creatinine", hops.d_labitems$label) & grepl("Blood", hops.d_labitems$fluid))]
weight_ids <- icu.d_items$itemid[grep("Weight", icu.d_items$label)]
icu.d_items[grep("Weight", icu.d_items$label), ]

core.patient <- fread("data/MIMIC/mimic-iv-1.0/core/patients.csv")

icu.chartevents_weight <- fread("data/MIMIC/reduced/icu_chartevents_weight.csv")

icu.chartevents_weight_first <- icu.chartevents_weight %>% group_by(subject_id) %>% 
  filter(charttime == min(charttime)) %>% ungroup() %>% 
  rename(weight = value) %>% 
  select(subject_id, weight)

rm(icu.chartevents_weight)

icu.chartevents_scr <- fread("data/MIMIC/reduced/icu_chartevents_scr.csv")

icu.chartevents_first <- icu.chartevents_scr %>% group_by(subject_id) %>% 
  filter(charttime == min(charttime)) %>% ungroup() %>% 
  rename(scr = value) %>% 
  select(subject_id, scr)

rm(icu.chartevents_scr)


total_data <- core.patient %>% left_join(icu.chartevents_weight_first) %>% 
  left_join(icu.chartevents_first) %>% as.data.frame()


w_age_data <- total_data %>% filter(!is.na(weight)) %>% 
  filter(!is.na(scr))# %>% slice_head(n = 10000)
######################
# try copula fit

library(rvinecopulib)
source("/home/laura/copulas/copulas_pharmacology/scripts/functions/estimate_vinecopula_from_data.R")
source("/home/laura/copulas/copulas_pharmacology/scripts/functions/functions.R")
library(kde1d)

mimic_copula <- estimate_vinecopula_from_data(w_age_data, variables_of_interest = c("anchor_age", "weight", "scr"), keep_data = FALSE, family_set = "parametric")

contour(mimic_copula)

new_patients <- simulate(mimic_copula, n = 1000)


#try to read the full data at once
icu.chartevents <- fread("data/MIMIC/mimic-iv-1.0/icu/chartevents.csv")

reduced_chartevents <- icu.chartevents[, c(1, 4, 6, 7)]
rm(icu.chartevents)
head(reduced_chartevents)


########
core.patient <- fread("data/MIMIC/mimic-iv-1.0/core/patients.csv")
icu.chartevents_weight <- fread("data/MIMIC/reduced/icu_chartevents_weight.csv")
icu.d_items <- fread("data/MIMIC/mimic-iv-1.0/icu/d_items.csv")

icu.chartevents_weight_first <- icu.chartevents_weight %>% group_by(subject_id) %>% 
  filter(charttime == min(charttime)) %>% ungroup() %>% 
  rename(weight = value) %>% 
  select(subject_id, weight)
rm(icu.chartevents_weight)

total_data <- core.patient %>% left_join(icu.chartevents_weight_first)

lab_ids <- icu.d_items %>% filter(category == "Labs" & !grepl("^Z", label)) %>% select(itemid) %>% unlist
names(lab_ids) <- NULL

for (item_id in lab_ids[1:10]) {
  one_event_only <- reduced_chartevents %>% filter(itemid == item_id) %>% 
    group_by(subject_id) %>% 
    filter(charttime == min(charttime)) %>% ungroup() %>%  
    select(!!quo_name(icu.d_items$label[icu.d_items$itemid == item_id]) := value, subject_id)%>% 
    distinct(subject_id, .keep_all = TRUE)
  
  total_data <- total_data %>% left_join(one_event_only, all.x = TRUE, all.y = FALSE)
}
number_na <- rowSums(is.na(total_data[, 7:17]))


reduced_rows <- total_data %>% filter(number_na < 11) %>% 
  distinct(subject_id, .keep_all = T)
save(reduced_rows, file = "data/MIMIC/data_wide_11vars.Rdata")
load("data/MIMIC/data_wide_11vars.Rdata")

for (item_id in lab_ids[11:30]) {
  one_event_only <- reduced_chartevents %>% filter(itemid == item_id) %>% 
    group_by(subject_id) %>% 
    filter(charttime == min(charttime)) %>% ungroup() %>%  
    select(!!quo_name(icu.d_items$label[icu.d_items$itemid == item_id]) := value, subject_id) %>% 
    distinct(subject_id, .keep_all = TRUE)
  reduced_rows <- reduced_rows %>% left_join(one_event_only, )
  print(dim(reduced_rows))
}

reduced_rows2 <- reduced_rows %>% 
  distinct(subject_id, .keep_all = T)

save(reduced_rows2, file = "data/MIMIC/data_wide_31vars.Rdata")

for (item_id in lab_ids[31:100]) {
  one_event_only <- reduced_chartevents %>% filter(itemid == item_id) %>% 
    group_by(subject_id) %>% 
    filter(charttime == min(charttime)) %>% ungroup() %>%  
    select(!!quo_name(icu.d_items$label[icu.d_items$itemid == item_id]) := value, subject_id) %>% 
    distinct(subject_id, .keep_all = TRUE)
  reduced_rows2 <- reduced_rows2 %>% left_join(one_event_only, )
  print(dim(reduced_rows2))
}

reduced_rows3 <- reduced_rows2 %>% 
  distinct(subject_id, .keep_all = T)

save(reduced_rows3, file = "data/MIMIC/data_wide_100vars.Rdata")

for (item_id in lab_ids[101:127]) {
  one_event_only <- reduced_chartevents %>% filter(itemid == item_id) %>% 
    group_by(subject_id) %>% 
    filter(charttime == min(charttime)) %>% ungroup() %>%  
    select(!!quo_name(icu.d_items$label[icu.d_items$itemid == item_id]) := value, subject_id) %>% 
    distinct(subject_id, .keep_all = TRUE)
  reduced_rows3 <- reduced_rows3 %>% left_join(one_event_only, )
  print(dim(reduced_rows3))
}

reduced_rows4 <- reduced_rows3 %>% 
  distinct(subject_id, .keep_all = T)

save(reduced_rows4, file = "data/MIMIC/data_wide_127vars.Rdata")

###
#copula estimation
#mimic data
library(tidyverse)
load("/home/laura/copulas/copulas_pharmacology/data/MIMIC/data_wide_127vars.Rdata")
names(reduced_rows4)
str(reduced_rows4)
item_selection <- read.csv("/home/laura/copulas/copulas_pharmacology/data/MIMIC/labs cvh.csv")

item_selection[, "select"] <- ifelse(is.na(item_selection[, "select"] ), 0, item_selection[, "select"])
chosen_items <- setdiff(item_selection$label[item_selection$select == 1], c("Absolute Neutrophil Count", "Absolute Count - Retic"))
clean_MIMIC <- as.data.frame(apply(reduced_rows4[, c(3, 7:134)], 2, function(x) {
  x[x == "999999"] <- NA
  x <- as.numeric(x)
  return(x)})) %>%
  select(all_of(c(chosen_items, "weight", "anchor_age"))) 
  

cop_data_MIMIC <- clean_MIMIC %>% 
  filter(rowSums(!is.na(clean_MIMIC)) > 1)
  


source("/home/laura/copulas/copulas_pharmacology/scripts/functions/estimate_vinecopula_from_data.R")
source("/home/laura/copulas/copulas_pharmacology/scripts/functions/functions.R")
source("/home/laura/copulas/copulas_pharmacology/scripts/functions/plot_distributions.R")

start_time <- Sys.time()
cop_mimic <- estimate_vinecopula_from_data(dat = cop_data_MIMIC, 
                                           variables_of_interest = c(chosen_items, "weight", "anchor_age"), 
                                           keep_data = FALSE, family_set = "parametric", cores = 40)

save(cop_mimic, file = "copula_mimic_allvar_53kpatients_filtered.Rdata")
end_time <- Sys.time()
print(end_time - start_time)

load("copula_mimic_allvar_53kpatients_filtered.Rdata")
cop_mimic$variables_of_interest


#Integer variables
check_integer <- function(x) {
  all(x%%1 == 0, na.rm = TRUE)
}

discrete_variables <- sapply(cop_data_MIMIC, check_integer)

var_types_mimic <- ifelse(discrete_variables, "d", "c")

start_time <- Sys.time()
cop_mimic_disc <- estimate_vinecopula_from_data(dat = cop_data_MIMIC,
                                                var_types = var_types_mimic,
                                                keep_data = FALSE, family_set = "parametric", cores = 40)

save(cop_mimic_disc, file = "copula_mimic_allvar_53kpatients_filtered_disc.Rdata")
end_time <- Sys.time()
print(end_time - start_time)

test_sim <- simulate(cop_mimic_disc, 100)
set.seed(123)
large_sim <- simulate(cop_mimic_disc, 1000)
cop_sim <- get_statistics(large_sim)
obs <- get_statistics(cop_data_MIMIC)


plot_data <- cop_sim %>%
  filter(covariate == "all" ) %>% 
  mutate(key = paste(statistic, covariate, sep = "_")) %>% 
  left_join(obs %>% dplyr::rename(observed = value)) %>% 
  mutate(rel_value = (value - observed)/observed) %>% 
  mutate(diff_ratio = value/observed)

simulation_performance <- plot_data %>% 
  ggplot(aes(y = rel_value, x = key)) +
  geom_point() +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = NULL, y = "Relative error (simulated/observed)") +
  theme_bw()

simulation_performance2 <- plot_data %>% 
  ggplot(aes(y = rel_value)) +
  geom_boxplot() +
  geom_hline(yintercept = c(-0.2, 0.2), linetype = 2, color = "grey65") +
  geom_hline(yintercept = 0, linetype = 1, color = "black") +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = NULL, y = "Relative error (simulated/observed)") +
  scale_y_continuous(limits = c(-3, 3)) +
  theme_bw()



vars_of_interest <- c("`Hematocrit (serum)`", "`Potassium (serum)`")

dens_plot <- cop_data_MIMIC %>% 
  ggplot(mapping = aes_string(y = vars_of_interest[1], x = vars_of_interest[2])) +
  geom_density2d(alpha = 1, linetype = 2) +
  geom_density2d(data = large_sim, alpha = 1, linetype = 2, color = "red") +
  #geom_point(color = "red", alpha = 0.1) +
  theme_bw() +
  theme(legend.position = "none")

vars_of_interest <- c("`Brain Natiuretic Peptide (BNP)`", "`Total Protein`")

dens_plot2 <- cop_data_MIMIC %>% 
  ggplot(mapping = aes_string(y = vars_of_interest[1], x = vars_of_interest[2])) +
  geom_density2d(alpha = 1, linetype = 2) +
  geom_density2d(data = large_sim, alpha = 1, linetype = 2, color = "red") +
  #geom_point(color = "red", alpha = 0.1) +
  theme_bw() +
  theme(legend.position = "none")

pdf("/home/laura/copulas/copulas_pharmacology/results/figures/density_interesting.pdf")
print(dens_plot)
print(dens_plot2)
dev.off()

pdf("/home/laura/copulas/copulas_pharmacology/results/figures/all_densities.pdf", width = 50, height = 50)
plot_comparison_distribution_sim_obs_generic(sim_data = large_sim, obs_data = cop_data_MIMIC, 
                                             plot_type = "density", variables = paste0("`", names(large_sim), "`"))
dev.off()

##############################

icu.d_items %>% filter(category == "Labs") %>% arrange(label) %>% 
  write.csv(file = "data/d_labitems.csv", row.names = FALSE)


first_measurement <- icu.charterevents %>% 
  filter(V6 %in% weight_ids) %>% 
  mutate(starttime = as.POSIXct(starttime)) %>% 
  group_by(subject_id) %>%
  filter(starttime == min(starttime)) %>% ungroup()

rm(icu.inputevents)

information <- list()
for (filename in list.files(mimic_path, recursive = TRUE, pattern = "*\\.csv$")) {
  information[[filename]][["names"]] <- names(read.csv(paste0(mimic_path, filename), nrows = 1))
  lines_raw <- system(paste0("wc -l ", mimic_path, filename), intern = TRUE)
  #number of lines (removing header line)
  information[[filename]][["nr_lines"]] <- as.numeric(gsub(" .*", "", lines_raw)) - 1
}

save(information, file = "data/MIMIC/information.Rdata")
load("data/MIMIC/information.Rdata")
unlist(lapply(information, function(x) grep("weight", x[["names"]], value = TRUE)))


all_names <- lapply(information, function(x) x[["names"]])




icu.procedureevents <- read.csv("data/MIMIC/mimic-iv-1.0/icu/procedureevents.csv")

length(unique(icu.procedureevents$subject_id))

#files of interest


#load in dynamically
core.patient <- read.csv("data/MIMIC/mimic-iv-1.0/core/patients.csv")

hosp.labevents <- read.csv("data/MIMIC/mimic-iv-1.0/hosp/labevents.csv", nrows = 10, skip = 10, header = FALSE)

hops.d_labitems <- read.csv("data/MIMIC/mimic-iv-1.0/hosp/d_labitems.csv")
id_serum_creatinine <-  hops.d_labitems$itemid[(grepl("Creatinine", hops.d_labitems$label) & grepl("Blood", hops.d_labitems$fluid))]
id_creatinine_clearance <-  hops.d_labitems$itemid[grepl("Creatinine Clearance", hops.d_labitems$label)]



#add to core.patient
core.patient
hosp.labevents_names <- names(read.csv("data/MIMIC/mimic-iv-1.0/hosp/labevents.csv", nrows = 1))
lines_raw <- system(paste0("wc -l ", "data/MIMIC/mimic-iv-1.0/hosp/labevents.csv"), intern = TRUE)

#number of lines: 122103668 

all_serum <- NULL
n_file <- 100000
library(data.table)
#for (j in seq(0, 122103668, by = n_file)) {
for (j in seq(0, 1000000, by = n_file)) {
  only_serum <- NULL
  
  #for (i in seq(1, 122103668 - n_file, by = n_file)) {
  for (i in seq(1, 1000000 - n_file, by = n_file)) {
    hosp.labevents <- fread("data/MIMIC/mimic-iv-1.0/hosp/labevents.csv", nrows = n_file, skip = i + j, 
                               header = FALSE, stringsAsFactors = FALSE)
    creatinine_events <- hosp.labevents[hosp.labevents[, 5] %in% c(id_serum_creatinine, id_creatinine_clearance), ]
    only_serum <- rbind(only_serum, creatinine_events)
  }
  names(only_serum) <- hosp.labevents_names
  first_measurements <- only_serum %>% group_by(subject_id) %>%
    mutate(charttime = as.POSIXct(charttime)) %>% 
    filter(charttime == min(charttime))
  all_serum <- rbind.data.frame(all_serum, first_measurements)
}

all_serum <- all_serum %>% left_join(core.patient)
#itemid -> 5
#subject.id -> 2

hosp.labevents[, ]

creatinine_events <- hosp.labevents[hosp.labevents[, 5] %in% c(id_serum_creatinine, id_creatinine_clearance), ]
creatinine_events[, 2]


core.patient[unique(as.character(creatinine_events[, 2])), ] <- creatinine_events

rownames(core.patient) <- core.patient$subject_id

length(unique(core.patient$subject_id)) == length((core.patient$subject_id)) 
