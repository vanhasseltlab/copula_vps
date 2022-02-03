#Longitudinal analysis
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

#plot data
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

### Using polynomials
#remove all observations with baby (only looking at pregnancy)
data_clean <- data_reduced %>% 
  filter(BABY == 0)


#visualize polynomials for subset
data_clean %>% 
  filter(ID %in% c(1:3, 5:13)) %>%
  ggplot(aes(x = gest, y = Albumin)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
  facet_wrap(~ ID) +
  theme_bw()

fixed_lm_poly <- lm(data = data_clean, Albumin ~ poly(gest, 2, raw = TRUE))

fixed_lm_function <- function(x, lm_obj) {
  coefs <- coef(lm_obj)
  y <- coefs[1] + coefs[2]*x + coefs[3]*x^2
  return(y)
}

with(data_clean, plot(Albumin ~ gest))
lines(fixed_lm_function(sort(unique(data_clean$gest)), fixed_lm_poly) ~ sort(unique(data_clean$gest)), col = "red")


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

#fit copula curves

hist(individual_coefs$b0)
hist(individual_coefs$b1)
hist(individual_coefs$b2)
marg_beta0 <- estimate_parametric_marginal(individual_coefs$b0, "norm")
marg_beta1 <- estimate_parametric_marginal(individual_coefs$b1, "norm")
marg_beta2 <- estimate_parametric_marginal(individual_coefs$b2, "norm")

uniform_coefs <- data.frame(b0 = marg_beta0$pit(individual_coefs$b0),
                            b1 = marg_beta1$pit(individual_coefs$b1),
                            b2 = marg_beta2$pit(individual_coefs$b2))
vine_coefs <- vinecop(uniform_coefs)
contour(vine_coefs)
pairs_copula_data(uniform_coefs)
cor(uniform_coefs, method = "kendall")

pdf("results/figures/longitudinal_parameter_copula.pdf")
contour(vine_coefs)
pairs_copula_data(uniform_coefs)
dev.off()


#Simulate polynomials using copula
set.seed(123)
m <- 123
sim_poly_unif <- as.data.frame(rvinecop(n = m, vine_coefs))

sim_poly <- data.frame(b0 = marg_beta0$pdf(sim_poly_unif$b0),
                       b1 = marg_beta1$pdf(sim_poly_unif$b1),
                       b2 = marg_beta2$pdf(sim_poly_unif$b2))

times <- seq(min(data_clean$gest), max(data_clean$gest), by = 1/7)
poly_df_sim <- sim_poly %>% 
  rownames_to_column("ID") %>%
  right_join(as.data.frame(expand_grid(gest = times, ID = rownames(sim_poly)))) %>% 
  mutate(pred =  b0 + b1*gest + b2*gest^2)

poly_df_sim %>%
  filter(ID %in% 1:12) %>% 
  ggplot(aes(x = gest)) +
  geom_line(aes(y = pred), color = "red") +
  facet_wrap(~ ID) +
  theme_bw()




## create copula for gest and residual noise
poly_df %>% 
  ggplot(aes(x = gest, y = res)) +
  geom_point() +
  theme_bw()
hist(poly_df$gest)
marg_gest <- estimate_spline_marginal(poly_df$gest)
check_fit_plot(poly_df$gest, marg_gest$density)

marg_res <- estimate_spline_marginal(poly_df$res[!is.na(poly_df$res)])
check_fit_plot(poly_df$res, marg_res$density)


uniform_data_res <- data.frame(res = marg_res$pit(poly_df$res[!is.na(poly_df$res)]),
                               gest = marg_gest$pit(poly_df$gest[!is.na(poly_df$res)]))
cop_res <- bicop(uniform_data_res)
contour(cop_res)

n <- 2000
m <- 123
set.seed(123)
test_sim <- rbicop(n = n, cop_res)
sim_org <- data.frame(gest = marg_gest$pdf(test_sim[, "gest"]),
                           res = marg_res$pdf(test_sim[, "res"]),
                           ID = rep(1:m, length.out = n)) %>% 
  arrange(ID, gest)

#combine residuals and polynomials

poly_df_sim <- sim_poly %>% 
  rownames_to_column("ID") %>%
  mutate(ID = as.integer(ID)) %>%
  right_join(sim_org) %>% 
  mutate(pred =  b0 + b1*gest + b2*gest^2) %>% 
  mutate(pred_res = pred + res)


poly_df_sim  %>% 
  filter(ID %in% 1:12) %>% 
  ggplot(aes(x = gest)) +
  geom_line(aes(y = pred), color = "red") +
  geom_point(aes(y = pred_res)) +
  facet_wrap(~ ID) +
  theme_bw()

pdf("results/figures/longitudinal_alb_polynomials_simulated.pdf", width = 20, height = 20)
print(poly_df_sim %>% 
        ggplot(aes(x = gest)) +
        geom_line(aes(y = Albumin), color = "red") +
        geom_point(aes(y = pred_res)) +
        facet_wrap(~ ID) +
        theme_bw())
dev.off()

poly_df %>% 
  ggplot(aes(x = gest, y = pred, color = as.factor(ID), group = ID)) +
  geom_line(show.legend = F) +
  theme_bw()

poly_df_sim %>% 
  ggplot(aes(x = gest, y = pred, color = as.factor(ID), group = ID)) +
  geom_line(show.legend = F) +
  theme_bw()


plot_simulated_curves <- poly_df %>% 
  select(ID, gest, pred) %>% 
  mutate(type = "observed") %>% 
  bind_rows(poly_df_sim %>% 
              select(ID, gest, pred) %>% 
              mutate(ID = as.numeric(ID) + 200,
                     type = "simulated")) %>% 
  ggplot(aes(x = gest, y = pred, color = as.factor(ID), group = ID)) +
  geom_line(show.legend = F, alpha = 0.5) +
  facet_grid(~ type) +
  theme_bw()

pdf("results/figures/longitudinal_alb_polynomials_both.pdf", width = 6, height = 4)
print(plot_simulated_curves)
dev.off()