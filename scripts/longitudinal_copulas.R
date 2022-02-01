#Longitudinal analysis
library(fitdistrplus)
library(tidyverse)
library(rvinecopulib)
library(actuar)
library(EnvStats)
library(mvtnorm)
library(kde1d)

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

#Remove (most) imputed values
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

hist(data_reduced$Albumin)
wide_albumin <- data_reduced %>% 
  select(ID, gest, Albumin) %>% 
  mutate(gest = round(gest)) %>% 
  group_by(ID, gest) %>% 
  summarize(Albumin = mean(Albumin)) %>% ungroup() %>%
  pivot_wider(id_cols = ID, names_from = gest, names_prefix = "t_", values_from = Albumin)
order_cols <- suppressWarnings(order(as.numeric(str_remove(names(wide_albumin), "t_"))))
wide_albumin <- wide_albumin[, order_cols]

#Fit marginal distributions
alb_spline <- estimate_spline_marginal(data_reduced$Albumin)
alb_para <- estimate_parametric_marginal(data_reduced$Albumin[!is.na(data_reduced$Albumin)], "norm")
check_fit_plot(data_reduced$Albumin, alb_spline$density)
check_fit_plot(data_reduced$Albumin, alb_para$density)
alb_spline$dist$nobs
alb_para$dist$aic
cormat <- cor(wide_albumin[, -ncol(wide_albumin)], use = "pairwise.complete.obs")

melted_cormat <- reshape2::melt(cormat)
cor_plot <- ggplot(data = melted_cormat, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() +
  scale_fill_gradient2() +
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  labs(x = NULL, y = NULL) +
  theme_bw() +
  coord_equal()

pdf("results/figures/correlation_Albumin.pdf", width = 8)
print(cor_plot)
dev.off()


#get empirical unif marginals for each time point
uniform_data <- as.data.frame(apply(wide_albumin[, !names(wide_albumin) %in% "ID"], 2, 
                                    function(x) {rank(x, na.last = "keep")/(sum(!is.na(x)) + 1)})) %>% 
  select(which(colSums(!is.na(.)) >1 ))

vine_longitudinal <- vinecop(uniform_data, family_set = "parametric")

sum_vine <- summary(vine_longitudinal) #%>% filter(family != "indep")
vine_model <- suppressMessages(
  sum_vine %>% unnest_wider(conditioned) %>%
    rename(var1 = `...1`, var2 = `...2`) %>% 
  left_join(data.frame(var_name1 = vine_longitudinal$names, 
                       var1 = 1:length(vine_longitudinal$names),
                       time1 = as.numeric(gsub("t_", "", vine_longitudinal$names)))) %>% 
    left_join(data.frame(var_name2 = vine_longitudinal$names, 
                         var2 = 1:length(vine_longitudinal$names),
                         time2 = as.numeric(gsub("t_", "", vine_longitudinal$names))))
)


pdf("results/figures/copula_Albumin.pdf", height = 2.7, width = 35)
contour(vine_longitudinal, tree = 1:2)

dev.off()
pdf("results/figures/copula_Albumin2.pdf", height = 20, width = 14)
plot(vine_longitudinal, tree = 1:2, edge_labels = "family")

dev.off()


library(igraph)

graph_data <- vine_model %>% filter(family != "indep") %>% 
  select(var_name1, var_name2, time1, time2, family, tau)
vertice_list <- unique(c(graph_data$var_name1, graph_data$var_name2))
vertice_list <- vertice_list[order(as.numeric(gsub("t_", "", vertice_list)))]

g_Albumin <- graph_from_data_frame(d = graph_data, directed = F, vertices = vertice_list)
set.seed(123)
plot(x = g_Albumin)

g_Albumin %>%
  #add_layout_(in_circle()) %>%
  add_layout_(as_tree()) %>%
  plot(vertex.size = 1)



### Using polynomials

library(lme4)
data_reduced %>% 
  filter(ID %in% c(1:3, 5:13)) %>% 
  ggplot(aes(x = gest, y = Albumin)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = F) +
  facet_wrap(~ ID) +
  theme_bw()

data_reduced %>% 
  filter(ID %in% c(1:3, 5:13)) %>% 
  ggplot(aes(x = gest, y = Albumin)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y ~ x, se = F) +
  facet_grid(ID ~ BABY, scales = "free_x") +
  theme_bw()

library(splines)

fixed_lm_poly <- lm(data = data_reduced, Albumin ~ poly(gest, 2, raw = TRUE))
fixed_lm_bs <- lm(data = data_reduced, Albumin ~ bs(gest, knots = 1, degree = 3))

fixed_lm_function <- function(x, lm_obj) {
  coefs <- coef(lm_obj)
  y <- coefs[1] + coefs[2]*x + coefs[3]*x^2
  return(y)
}

with(data_reduced, plot(Albumin ~ gest))
lines(fixed_lm_function(sort(unique(data_reduced$gest)), fixed_lm_poly) ~ sort(unique(data_reduced$gest)), col = "red")

plot(bs(data_reduced$gest, knots = 1, degree = 3))

spline_gest <- bs(data_reduced$gest, knots = 1, degree = 3)
coefs_bs <- coef(fixed_lm_bs)

gest_x <- seq(2/7, 526/7, by = 1/7)
with(data_reduced, plot(Albumin ~ gest))
lines(gest_x, predict(fixed_lm_bs, data.frame(gest = gest_x)), col = "red")
coefs_bs[1] + coefs_bs[2]*spline_gest[, 1]

#try with baby

fixed_lm <- lm(data = data_reduced, Albumin ~  gest*BABY + I(gest^2)*BABY +I(gest^3)*BABY )

fixed_lm <- lm(data = data_reduced, Albumin ~  gest*BABY + I(gest^2):BABY)

with(data_reduced, {
  plot(Albumin ~ gest, pch = 1 + BABY)
  points(gest, predict(fixed_lm, data_reduced[, c("gest", "BABY")]), col = "red", pch = 1 + BABY)
})

fixed_lm <- lm(data = data_reduced, log(Albumin) ~  gest*BABY + I(gest^2):BABY)

with(data_reduced, {
  plot(log(Albumin) ~ gest, pch = 1 + BABY)
  points(gest, predict(fixed_lm, data_reduced[, c("gest", "BABY")]), col = "red", pch = 1 + BABY)
})

re_lm <- lmer(data = data_reduced,
              Albumin ~  gest*BABY + I(gest^2)*BABY +I(gest^3)*BABY + (1 + gest*BABY + I(gest^2)*BABY +I(gest^3)*BABY|ID), REML = TRUE)

sum_model <- summary(re_lm)
individual_coefs <- coef(re_lm)$ID

#make nice names
make_nice_names <- function(x, betas) {
  x[x == "(Intercept)"] <- "b0"
  x <- str_remove(x, fixed("I("))
  x <- str_remove(x, fixed(")"))
  x <- str_replace(x, fixed("^"),  "_")
  for (i in 1:length(betas)) {
    x <- str_replace(x, betas[i],  paste0("b", i))
  }
  return(x)
}

names(individual_coefs) <- make_nice_names(names(individual_coefs), c("gest", "BABY"))
# b0 + b1*gest + b2*BABY + b1_2*gest^2 + b1_3*gest^3 + `b1:b2`*gest*BABY + `b2:b1_2`*gest^2*BABY + `b2:b1_3`*gest^3*BABY

poly_df <- individual_coefs %>% 
  rownames_to_column("ID") %>%
  mutate(ID = as.integer(ID)) %>%
  right_join(data_reduced %>% select(ID, Albumin, gest, BABY)) %>% 
  mutate(pred =  b0 + b1*gest + b2*BABY + b1_2*gest^2 + b1_3*gest^3 + `b1:b2`*gest*BABY + `b2:b1_2`*gest^2*BABY + `b2:b1_3`*gest^3*BABY) %>% 
  mutate(res = Albumin - pred)

poly_df %>% 
  filter(ID %in% 1:12) %>% 
  ggplot(aes(x = gest)) +
  geom_line(aes(y = pred), color = "red") +
  geom_point(aes(y = Albumin)) +
  #geom_point(aes(y = fake_pred), color = "blue") +
  facet_wrap(~ ID) +
  theme_bw()

###
#simpler:
re_lm <- lmer(data = data_reduced,
              Albumin ~  gest*BABY + (1 + gest*BABY|ID), REML = TRUE)

#b0 + b1*gest + b2*BABY + b1:b2*gest*BABY

sum_model <- summary(re_lm)
individual_coefs <- coef(re_lm)$ID
names(individual_coefs) <- make_nice_names(names(individual_coefs), c("gest", "BABY"))

poly_df <- individual_coefs %>% 
  rownames_to_column("ID") %>%
  mutate(ID = as.integer(ID)) %>%
  right_join(data_reduced %>% select(ID, Albumin, gest, BABY)) %>% 
  mutate(pred =  b0 + b1*gest + b2*BABY + `b1:b2`*gest*BABY) %>% 
  mutate(res = Albumin - pred)

poly_df %>% 
  filter(ID %in% 20:50) %>% 
  ggplot(aes(x = gest)) +
  geom_line(aes(y = pred), color = "red") +
  geom_point(aes(y = Albumin)) +
  facet_wrap(~ ID) +
  theme_bw()


#mixed effect model
re_lm <- lmer(data = data_reduced,
              Albumin ~ poly(gest, 2, raw = TRUE) + (1 + poly(gest,2, raw = TRUE)|ID), REML = TRUE)
sum_model <- summary(re_lm)
individual_coefs <- coef(re_lm)$ID
names(individual_coefs) <- paste0("beta", 0:2)

#plot results
set.seed(123)
noise <- rnorm(nrow(data_reduced), mean = 0, sd = summary(re_lm)$sigma)

poly_df <- individual_coefs %>% 
  rownames_to_column("ID") %>%
  mutate(ID = as.integer(ID)) %>%
  right_join(data_reduced %>% select(ID, Albumin, gest)) %>% 
  mutate(pred =  beta0 + beta1*gest + beta2*gest^2) %>% 
  mutate(res = Albumin - pred) %>% 
  mutate(fake_pred = pred + noise)


poly_df %>% 
  filter(ID %in% 1:12) %>% 
  ggplot(aes(x = gest)) +
  geom_line(aes(y = pred), color = "red") +
  geom_point(aes(y = Albumin)) +
  #geom_point(aes(y = fake_pred), color = "blue") +
  facet_wrap(~ ID) +
  theme_bw()

poly_df %>% 
  ggplot(aes(group = ID, y = res, x = ID)) +
  geom_boxplot(aes(color = ID%%2 == 0)) +
  theme_bw()

poly_df_sim %>% 
  ggplot(aes(group = ID, y = res, x = ID)) +
  geom_boxplot(aes(color = ID%%2 == 0)) +
  theme_bw()

check_diffs <- poly_df %>% group_by(ID) %>% 
  summarize(mean = mean(res), sd = sd(res),
            median = median(res))


pdf("results/figures/longitudinal_alb_polynomials.pdf", width = 20, height = 20)
print(poly_df %>% 
        ggplot(aes(x = gest)) +
        geom_line(aes(y = pred), color = "red") +
        geom_point(aes(y = Albumin)) +
        #geom_point(aes(y = fake_pred), color = "blue", alpha = 0.5, shape = 17) +
        facet_wrap(~ ID) +
        theme_bw())
dev.off()


##residuals
poly_df %>% 
  ggplot(aes(x = gest, y = res)) +
  geom_point() +
  theme_bw()
hist(poly_df$gest)
marg_gest <- estimate_spline_marginal(poly_df$gest)
check_fit_plot(poly_df$gest, marg_gest$density)


marg_res <- estimate_parametric_marginal(poly_df$res[!is.na(poly_df$res)], "norm")
#marg_res <- estimate_parametric_marginal(poly_df$res[!is.na(poly_df$res)], "t",
 #                                        param = list(estimate = c()))

marg_id <-  estimate_parametric_marginal(poly_df$ID, "unif")
check_fit_plot(poly_df$res, marg_res$density)

uniform_data_res <- data.frame(res = marg_res$pit(poly_df$res),
                               gest = marg_gest$pit(poly_df$gest))
cop_res <- bicop(uniform_data_res)
contour(cop_res)

n <- 2000
m <- 120
set.seed(123)
test_sim <- rbicop(n = n, cop_res)
test_sim_org <- data.frame(gest = marg_gest$pdf(test_sim[, "gest"]),
                           res = marg_res$pdf(test_sim[, "res"]),
                           ID = rep(1:m, length.out = n))

sim_org <- test_sim_org %>% 
  arrange(ID, gest)


#fit copula curves

hist(individual_coefs$beta0)
hist(individual_coefs$beta1)
hist(individual_coefs$beta2)
marg_beta0 <- estimate_parametric_marginal(individual_coefs$beta0, "norm")
marg_beta1 <- estimate_parametric_marginal(individual_coefs$beta1, "norm")
marg_beta2 <- estimate_parametric_marginal(individual_coefs$beta2, "norm")

uniform_coefs <- data.frame(beta0 = marg_beta0$pit(individual_coefs$beta0),
                            beta1 = marg_beta1$pit(individual_coefs$beta1),
                            beta2 = marg_beta2$pit(individual_coefs$beta2))
vine_coefs <- vinecop(uniform_coefs)
contour(vine_coefs)
pairs_copula_data(uniform_coefs)
cor(uniform_coefs, method = "kendall")

pdf("results/figures/longitudinal_parameter_copula.pdf")
contour(vine_coefs)
pairs_copula_data(uniform_coefs)
dev.off()


set.seed(123)
sim_poly_unif <- as.data.frame(rvinecop(n = m, vine_coefs))

sim_poly <- data.frame(beta0 = marg_beta0$pdf(sim_poly_unif$beta0),
                       beta1 = marg_beta1$pdf(sim_poly_unif$beta1),
                       beta2 = marg_beta2$pdf(sim_poly_unif$beta2))

#combine residuals and polynomials

sim_org


poly_df_sim <- sim_poly %>% 
  rownames_to_column("ID") %>%
  mutate(ID = as.integer(ID)) %>%
  right_join(sim_org) %>% 
  mutate(pred =  beta0 + beta1*gest + beta2*gest^2) %>% 
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
        geom_line(aes(y = pred), color = "red") +
        geom_point(aes(y = pred_res)) +
        facet_wrap(~ ID) +
        theme_bw())
dev.off()


names(poly_df)
names(poly_df_sim)

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
              mutate(ID = ID + 200,
                     type = "simulated")) %>% 
  ggplot(aes(x = gest, y = pred, color = as.factor(ID), group = ID)) +
  geom_line(show.legend = F, alpha = 0.5) +
  facet_grid(~ type) +
  theme_bw()

pdf("results/figures/longitudinal_alb_polynomials_both.pdf", width = 6, height = 4)
print(plot_simulated_curves)
dev.off()
