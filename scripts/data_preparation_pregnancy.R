## Data preparation ##
# Pregnancy data
# Patel, J. P. et al. Circulation 128, 1462–1469 (2013).

##### - Libraries - #####
library(tidyverse)

#### Data ####
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
