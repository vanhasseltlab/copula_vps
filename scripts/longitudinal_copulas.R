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
source("scripts/functions/plot_distributions.R")
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

#mixed effect model
re_lm <- lmer(data = data_clean,
              Albumin ~ poly(gest, 2, raw = TRUE) + (1 + poly(gest,2, raw = TRUE)|ID), REML = TRUE)
sum_model <- summary(re_lm)
individual_coefs <- coef(re_lm)$ID
names(individual_coefs) <- paste0("b", 0:2)

#fit copula curves
marg_beta0 <- estimate_parametric_marginal(individual_coefs$b0, "norm")
marg_beta1 <- estimate_parametric_marginal(individual_coefs$b1, "norm")
marg_beta2 <- estimate_parametric_marginal(individual_coefs$b2, "norm")

uniform_coefs <- data.frame(b0 = marg_beta0$pit(individual_coefs$b0),
                            b1 = marg_beta1$pit(individual_coefs$b1),
                            b2 = marg_beta2$pit(individual_coefs$b2))
vine_coefs_alb <- vinecop(uniform_coefs)

pdf("results/figures/longitudinal_parameter_copula_albumin.pdf", height = 4, width = 4)
contour(vine_coefs_alb)
pairs_copula_data(uniform_coefs)
dev.off()

#fit copula splinecurves
marg_beta0_spl <- estimate_spline_marginal(individual_coefs$b0)
marg_beta1_spl <- estimate_spline_marginal(individual_coefs$b1)
marg_beta2_spl <- estimate_spline_marginal(individual_coefs$b2)

uniform_coefs_spl <- data.frame(b0 = marg_beta0_spl$pit(individual_coefs$b0),
                            b1 = marg_beta1_spl$pit(individual_coefs$b1),
                            b2 = marg_beta2_spl$pit(individual_coefs$b2))
vine_coefs_alb_spl <- vinecop(uniform_coefs_spl)

pdf("results/figures/longitudinal_parameter_copula_albumin_spl.pdf", height = 4, width = 4)
contour(vine_coefs_alb_spl)
pairs_copula_data(uniform_coefs_spl)
dev.off()



#Simulate polynomials using copula
set.seed(123)
m <- 1000
sim_poly_unif <- as.data.frame(rvinecop(n = m, vine_coefs_alb))

sim_poly <- data.frame(b0 = marg_beta0$pdf(sim_poly_unif$b0),
                       b1 = marg_beta1$pdf(sim_poly_unif$b1),
                       b2 = marg_beta2$pdf(sim_poly_unif$b2))

times <- seq(min(data_clean$gest), max(data_clean$gest), by = 1/7)
poly_df_sim <- sim_poly %>% 
  rownames_to_column("ID") %>%
  right_join(as.data.frame(expand_grid(gest = times, ID = rownames(sim_poly)))) %>% 
  mutate(pred =  b0 + b1*gest + b2*gest^2)

plot_simulated_curves <- poly_df %>% 
  select(ID, gest, pred) %>% 
  mutate(type = "observed") %>% 
  bind_rows(poly_df_sim %>% 
              select(ID, gest, pred) %>% 
              mutate(ID = as.numeric(ID) + 200,
                     type = "simulated")) %>% 
  ggplot(aes(x = gest, y = pred, group = ID)) +
  geom_line(show.legend = F, alpha = 0.5, color = "grey15") +
  facet_grid(~ type) +
  labs(x = "Gestational age", y = "Albumin concentration") +
  theme_bw()

pdf("results/figures/longitudinal_alb_comparison_sim_obs.pdf", width = 8, height = 6)
plot_comparison_distribution_sim_obs_generic(obs_data = individual_coefs, sim_data = sim_poly, 
                                             variables = c("b0", "b1", "b2"), plot_type = "density")
dev.off()

pdf("results/figures/longitudinal_alb_polynomials_both.pdf", width = 6, height = 4)
print(plot_simulated_curves)
dev.off()

set.seed(123)
m <- 1000
sim_poly_unif_spl <- as.data.frame(rvinecop(n = m, vine_coefs_alb_spl))

sim_poly_spl <- data.frame(b0 = marg_beta0_spl$pdf(sim_poly_unif_spl$b0),
                       b1 = marg_beta1_spl$pdf(sim_poly_unif_spl$b1),
                       b2 = marg_beta2_spl$pdf(sim_poly_unif_spl$b2))


pdf("results/figures/longitudinal_alb_comparison_sim_obs_spl.pdf", width = 8, height = 6)
plot_comparison_distribution_sim_obs_generic(obs_data = individual_coefs, sim_data = sim_poly_spl, 
                                             variables = c("b0", "b1", "b2"), plot_type = "density")
dev.off()



#Polynomials for multiple variables
variables_of_interest <- c("Platelets", "SCr", "Albumin", "Bilirubin", "Lymphocytes", "Neutrophils")

coef_data_raw <- as.data.frame(matrix(NA, nrow = length(unique(data_clean$ID)), ncol = 3*length(variables_of_interest)))
names(coef_data_raw) <- paste0("b", 0:2, "_", rep(variables_of_interest, each = 3))

for (variable in variables_of_interest) {
  formula_v <- as.formula(paste(variable, "~ poly(gest, 2, raw = TRUE) + (1 + poly(gest,2, raw = TRUE)|ID)"))
  
  re_lm <- lmer(data = data_clean, formula_v, REML = TRUE)
  individual_coefs <- coef(re_lm)$ID
  
  coef_data_raw[grep(paste0("_", variable), names(coef_data_raw))] <-  coef(re_lm)$ID
}

coef_data_unif <- coef_data_raw
marginals <- list()
for (variable in names(coef_data_unif)) {
  marginals[[variable]] <- estimate_spline_marginal(coef_data_unif[, variable])
  ind_na <- is.na(coef_data_unif[, variable])
  coef_data_unif[!ind_na, variable] <- marginals[[variable]]$pit(coef_data_unif[!ind_na, variable])
}



vine_coefs <- vinecop(coef_data_unif)
plot(vine_coefs, var_names = "use", tree = 1:2)
pdf("results/figures/longitudinal_parameter_copula.pdf", width = 20, height = 20)
contour(vine_coefs)
pairs_copula_data(coef_data_unif)
dev.off()

#Simulation
set.seed(123)
m <- 123
sim_unif <- as.data.frame(rvinecop(n = m, vine_coefs))

sim_data <- sim_unif
for (variable in names(sim_data)) {
  ind_na <- is.na(sim_data[, variable])
  sim_data[!ind_na, variable] <- marginals[[variable]]$pdf(sim_data[!ind_na, variable])
}

gest_times <- seq(min(data_clean$gest), max(data_clean$gest), by = 1/7)
df_sim <- sim_data %>% 
  rownames_to_column("ID") %>%
  right_join(as.data.frame(expand_grid(gest = gest_times, ID = rownames(sim_data))))

for (variable in variables_of_interest) {
  col_ind <- grep(paste0("_", variable), names(df_sim))
  df_sim[, variable] <- df_sim[, col_ind[1]] + df_sim[, col_ind[2]]*df_sim$gest + df_sim[, col_ind[3]]*df_sim$gest^2
}


df_poly <- coef_data_raw %>% 
  rownames_to_column("ID") %>%
  right_join(as.data.frame(expand_grid(gest = gest_times, ID = rownames(coef_data_raw))))

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

#simulated polynomials
plot_sim <- df_sim %>% 
  select(all_of(c("ID", "gest", variables_of_interest))) %>% 
  pivot_longer(-c(ID, gest), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw()

#estimated polynomials from observed data
plot_curves <- df_poly %>% 
  select(all_of(c("ID", "gest", variables_of_interest))) %>% 
  pivot_longer(-c(ID, gest), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  facet_wrap(~ variable, ncol = 1, scales = "free_y", strip.position = "right") + 
  theme_bw()


#side-by-side comparison of the estimated and simulated curves
df_sim %>% 
  mutate(type = "simulation") %>% 
  bind_rows(df_poly %>% mutate(type = "observed")) %>% 
  select(all_of(c("ID", "gest", variables_of_interest, "type"))) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey35") +
  geom_smooth(color = "red", se = F) +
  facet_grid(variable ~ type, scales = "free_y") + 
  theme_bw()

df_comparing_sim_org <- df_sim %>% 
  mutate(type = "Simulated") %>% 
  bind_rows(df_poly %>% mutate(type = "Observed")) %>% 
  select(all_of(c("ID", "gest", variables_of_interest, "type"))) %>% 
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

df_comparing_sim_org <- df_sim %>% 
  mutate(type = "Simulated") %>% 
  bind_rows(df_poly %>% mutate(type = "Observed")) %>% 
  select(all_of(c("ID", "gest", variables_of_interest, "type"))) %>% 
  pivot_longer(-c(ID, gest, type), names_to = "variable") %>% 
  ggplot(aes(x = gest, y = value)) +
  geom_line(aes(group = as.factor(ID)), color = "grey15", alpha = 0.5) +
  geom_smooth(color = "red", se = F) +
  labs(y = NULL, x = "Gestational age (weeks)") +
  facet_grid(type ~ variable, scales = "free_y", switch = "y") + 
  theme_bw() +
  theme(strip.placement = "outside")

pdf("results/figures/longitudinal_polynomials_flipped.pdf", width = 9, height = 6)
print(df_comparing_sim_org)
dev.off()


#correlations at different timepoints
cov(df_sim %>% filter(gest == gest_times[10]) %>% select(variables_of_interest))
cov(df_poly %>% filter(gest == gest_times[10]) %>% select(variables_of_interest))



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


#Try long format - variables: variables_of_interest, gest, ID

df_clean_long <- data_clean %>% 
  select(all_of(c("ID", "gest", variables_of_interest))) %>% 
  pivot_longer(-c(ID, gest), names_to = "variable") %>% 
  left_join(df_poly %>% 
    select(all_of(c("ID", variables_poly))) %>% 
    mutate(ID = as.integer(ID)) %>% distinct() %>% 
    pivot_longer(-ID) %>% 
    separate(name, into = c("coefficient", "variable")) %>% 
  pivot_wider(c(ID, variable), names_from = coefficient)) %>% 
  mutate(pred = b0 + b1*gest + b2*gest^2) %>% 
  mutate(res = value - pred,
         scaled_res = (value - pred)/value)

df_clean_long %>% 
  ggplot(aes(x = gest, y = scaled_res)) +
  geom_point(show.legend = F) +
  geom_smooth(color = "red") +
  geom_hline(yintercept = 0, color = "blue", linetype = 2) +
  geom_hline(yintercept = c(-2,2), color = "light blue", linetype = 2) +
  facet_grid(variable ~., scales = "free_y") +
  theme_bw()


#check residuals distributions
residuals <- df_clean_long %>% 
  select(ID, gest, variable, res) %>% 
  distinct(ID, gest, variable, .keep_all = TRUE) %>% 
  pivot_wider(c(ID, gest), names_from = "variable", values_from = "res")


marg_gest <- estimate_spline_marginal(residuals$gest)
check_fit_plot(residuals$gest, marg_gest$density)

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

