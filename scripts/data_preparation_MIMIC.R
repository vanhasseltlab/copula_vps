# preparation MIMIC data
##### - Libraries - #####
library(tidyverse)
library(data.table)

##### - Data - #####
# unzip MIMIC files (runs only if files still unzipped (.gz))
mimic_path <- "data/MIMIC/mimic-iv-1.0/"
for (filename in list.files(mimic_path, recursive = TRUE, pattern = "*\\.gz$")) {
  R.utils::gunzip(paste0(mimic_path, filename), remove = TRUE)
}

# Read data (large data)
cat("start reading data \n")
icu.chartevents <- fread("data/MIMIC/mimic-iv-1.0/icu/chartevents.csv")
reduced_chartevents <- icu.chartevents[, c(1, 4, 6, 7)]
rm(icu.chartevents)

core.patient <- fread("data/MIMIC/mimic-iv-1.0/core/patients.csv")
icu.chartevents_weight <- fread("data/MIMIC/reduced/icu_chartevents_weight.csv") #from data_reduction.py
icu.d_items <- fread("data/MIMIC/mimic-iv-1.0/icu/d_items.csv")

##### - Data preparation - #####
icu.chartevents_weight_first <- icu.chartevents_weight %>% group_by(subject_id) %>% 
  filter(charttime == min(charttime)) %>% ungroup() %>% 
  rename(weight = value) %>% 
  select(subject_id, weight)
rm(icu.chartevents_weight)

total_data <- core.patient %>% left_join(icu.chartevents_weight_first)

lab_ids <- icu.d_items %>% filter(category == "Labs" & !grepl("^Z", label)) %>% select(itemid) %>% unlist
names(lab_ids) <- NULL

#loop through all lab_ids to add covariates in columns
cat("number of lab ids = ", length(lab_ids), "\n")
for (item_id in lab_ids) {
  one_event_only <- reduced_chartevents[reduced_chartevents$itemid == item_id, ] %>% 
    group_by(subject_id) %>% 
    filter(charttime == min(charttime)) %>% ungroup() %>%  
    select(!!quo_name(icu.d_items$label[icu.d_items$itemid == item_id]) := value, subject_id) %>% 
    distinct(subject_id, .keep_all = TRUE)
  
  total_data <- total_data %>% left_join(one_event_only, all.x = TRUE, all.y = FALSE, by = "subject_id")
  
  #Remove patients with lacking most data after adding 10 items
  columns_added <- ncol(total_data)
  if (columns_added == 17) {
    number_na <- rowSums(is.na(total_data[, 7:17]))
    total_data <- total_data %>% filter(number_na < 11) %>% 
      distinct(subject_id, .keep_all = T)
  }
  cat("nr of columns = ", columns_added, "\r")
}
cat("\n")

reduced_rows4 <- total_data %>% 
  distinct(subject_id, .keep_all = T)

cat("saving data_wide_127vars.Rdata")
save(reduced_rows4, file = "data/MIMIC/data_wide_127vars.Rdata")


# Save information about each file in MIMIC database
information <- list()
for (filename in list.files(mimic_path, recursive = TRUE, pattern = "*\\.csv$")) {
  information[[filename]][["names"]] <- names(read.csv(paste0(mimic_path, filename), nrows = 1))
  lines_raw <- system(paste0("wc -l ", mimic_path, filename), intern = TRUE)
  #number of lines (removing header line)
  information[[filename]][["nr_lines"]] <- as.numeric(gsub(" .*", "", lines_raw)) - 1
}
cat("saving information.Rdata")
save(information, file = "data/MIMIC/information.Rdata")
