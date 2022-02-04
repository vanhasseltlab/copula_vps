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
vine_coefs <- vinecop(uniform_coefs)

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
  ggplot(aes(x = gest, y = res, color = as.factor(ID))) +
  geom_point(show.legend = F) +
  geom_line(aes(group = as.factor(ID)), show.legend = F) +
  theme_bw()
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


#correlations at different timepoints
cov(df_sim %>% filter(gest == gest_times[10]) %>% select(variables_of_interest))
cov(df_poly %>% filter(gest == gest_times[10]) %>% select(variables_of_interest))


gest_times
