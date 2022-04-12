#Longitudinal analysis plotting
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(EnvStats)
library(mvtnorm)
library(kde1d)
library(lme4)

source("scripts/functions/functions.R")
#read in pregnancy data
data_pregnancy_raw <- read.csv("data/Jig - pregnancy file August 2012 minus coag info.csv", 
                               row.names = NULL, na.strings = c("")) %>% 
  mutate(Neutrophils = as.numeric(str_replace(Neutrophils, ":", ".")))
data_reduced <- read.csv("data/clean/pregnancy_reduced.csv", row.names = NULL)

data_clean <- data_reduced %>% 
  filter(BABY == 0)

#explore some patients
data_head <- data_pregnancy_raw %>% filter(ID %in% 1:5)
data_head %>% 
  ggplot(aes(x = Platelets, y = SCr)) +
  geom_point(aes(color = as.factor(ID))) +
  theme_bw()

variables_of_interest <- c("Platelets", "SCr", "Albumin", "Bilirubin", "Lymphocytes", "Neutrophils")
plot_explore <- data_head %>% 
  select(all_of(c("ID", "gest", variables_of_interest))) %>% 
  pivot_longer(-c(ID, gest), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value, group = as.factor(ID), color = as.factor(ID))) +
  geom_line() +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw()

#plot reduced data - patient profiles
variables_of_interest <- c("Platelets", "SCr", "Albumin", "Bilirubin", "Lymphocytes", "Neutrophils")
plot_explore <- data_reduced %>% 
  select(all_of(c("ID", "gest", variables_of_interest))) %>% 
  pivot_longer(-c(ID, gest), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw()

#select one variable for over time copula: Albumin
#plot Albumin
plot_Albumin_orig <- data_pregnancy_raw %>% 
  ggplot(aes(x = gest, y = Albumin)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  theme_bw()

plot_Albumin <- data_reduced %>% 
  ggplot(aes(x = gest, y = Albumin)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  theme_bw()

pdf("results/figures/Albumin_over_time.pdf", width = 5, height = 4)
plot_Albumin_orig
plot_Albumin
dev.off()


#visualize polynomials for subset
data_clean %>% 
  filter(ID %in% c(1:3, 5:13)) %>%
  ggplot(aes(x = gest, y = Albumin)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
  facet_wrap(~ ID) +
  theme_bw()

fixed_lm_function <- function(x, lm_obj) {
  coefs <- coef(lm_obj)
  y <- coefs[1] + coefs[2]*x + coefs[3]*x^2
  return(y)
}

with(data_clean, plot(Albumin ~ gest))
lines(fixed_lm_function(sort(unique(data_clean$gest)), fixed_lm_poly) ~ sort(unique(data_clean$gest)), col = "red")


#individual polynomial predictions:
#mixed effect model
re_lm <- lmer(data = data_clean,
              Albumin ~ poly(gest, 2, raw = TRUE) + (1 + poly(gest,2, raw = TRUE)|ID), REML = TRUE)
sum_model <- summary(re_lm)
individual_coefs <- coef(re_lm)$ID
names(individual_coefs) <- paste0("b", 0:2)

#plot results
poly_df <- individual_coefs %>% 
  rownames_to_column("ID") %>%
  mutate(ID = as.integer(ID)) %>%
  right_join(data_clean %>% select(ID, Albumin, gest)) %>% 
  mutate(pred =  b0 + b1*gest + b2*gest^2) %>% 
  mutate(res = Albumin - pred)


poly_df %>% 
  filter(ID %in% 1:12) %>% 
  ggplot(aes(x = gest)) +
  geom_line(aes(y = pred), color = "red") +
  geom_point(aes(y = Albumin)) +
  facet_wrap(~ ID) +
  theme_bw()

pdf("results/figures/longitudinal_alb_polynomials.pdf", width = 20, height = 20)
print(poly_df %>% 
        ggplot(aes(x = gest)) +
        geom_line(aes(y = pred), color = "red") +
        geom_point(aes(y = Albumin)) +
        facet_wrap(~ ID) +
        theme_bw())
dev.off()

something = "b"
test_exp <- expression(paste(hat("y"), "=", something))


test_exp <- bquote(italic(hat("y"))~"="~.(something)[0]~"+"~.(round(something)[1]*"\u00B7"*age~"+"~.(something)[2]*"\u00B7"*age^2))


plot(1:10, 1:10, main = test_exp)

eq <- ddply(df,.(group),lm_eqn)


poly_plot_df <- poly_df %>% 
  filter(ID %in% c(23, 45, 93, 15, 2, 9)) %>% 
  mutate(ID = paste0("ID", as.numeric(as.factor(ID)))) %>% 
  group_by(ID) %>% 
  mutate(form = as.character(as.expression(bquote(italic(hat("y"))~"="~.(round(b0[1], 2))~"+"~.(round(b1[1], 2))*"\u00B7"*age~"+"~.(formatC(b2[1], 1, format = "e"))*"\u00B7"*age^2)))) %>% 
  ungroup()


poly_plot_df_dist <- poly_plot_df %>% distinct(ID, form)

pdf("results/figures/longitudinal_alb_polynomials_example.pdf", width = 7, height = 4)
print(poly_plot_df %>% 
        ggplot(aes(x = gest)) +
        geom_text(aes(label = form), parse = TRUE, x = 22, y = 49, size = 2.3) +
        geom_line(aes(y = pred), color = "red") +
        geom_point(aes(y = Albumin)) +
        labs(y = "Albumin concentration", x = "Gestational age") +
        facet_wrap(~ ID) +
        theme_bw())
dev.off()

citation("rvinecopulib")
#copula

hist(individual_coefs$b0)
hist(individual_coefs$b1)
hist(individual_coefs$b2)