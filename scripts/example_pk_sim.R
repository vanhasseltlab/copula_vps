library(tidyverse)
library(RxODE)
library(DescTools)

# load function
source("./scripts/functions/run_Grimsley.R", local = T)

# import covariate data
dat_covar <- read_csv("./example_data/pediatrics_example_data.csv")

# run simulation
df_test <- run_grimsley(
  n = length(dat_covar$ID), # number of subjects
  wgt = dat_covar$BW/1000, # body weight input (kg)
  scr = dat_covar$CREA # serum creatinine (umol/l)
)

# plot
df_test %>% 
  ggplot(aes(x = time/60, y = conc, group = id)) +
  geom_line(color = "darkblue", alpha = 0.3) +
  scale_x_continuous(name = "Time (hours)", breaks = c(0:9)*8) +
  scale_y_continuous(name = "Concentration (mg/L)") +
  theme_classic()

# AUC 0-24 hours 
df_test %>% 
  group_by(id) %>% 
  summarise(
    AUC24 = AUC(
      x = time,
      y = conc,
      from = 0,
      to = 24*60,
      method = "trapezoid",
      na.rm = FALSE
    )
  ) %>% 
  ungroup() %>% 
  view()
