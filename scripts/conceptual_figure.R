#Making conceptual figure
library(rvinecopulib)
library(tidyverse)

#cop_distribution
#family: "gumbel", 0, parameters = c(3)
#marginals:   color = round(qunif(x[, 1], min = 1, max = 10)),
#             size  = qnorm(x[, 2], mean = 20, sd = 5)

# example observed data:
concept_data <- read.csv("example/simulated_data_conceptual_figure.csv", row.names = 1)

point_plot_intro <- concept_data %>% 
  ggplot(aes(x = size, y = color)) +
  geom_point(aes(color = color_hexa2, fill = color_hexa, size = size), 
             show.legend = FALSE, alpha = 0.8, shape = 16) +
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(limits = range(concept_data$size)*c(0.90, 1.05), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = range(concept_data$color)*c(0.6, 1.07), expand = expansion(mult = c(0, 0))) +
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())



cairo_pdf(file = "example/person_plot.pdf", width = 2.4, height = 2.4, onefile = T)
print(point_plot_intro)
dev.off()



#with density:
set.seed(123) #Not the right seed for the data example above
contour_sim_unif <- as.data.frame(rbicop(n = 10000, family = "gumbel", 0, parameters = c(3)))
contour_sim <- data.frame(color = round(qunif(contour_sim_unif[, 1], min = 1, max = 10)),
                          color_raw = qunif(contour_sim_unif[, 1], min = 1, max = 10),
                          size = qnorm(contour_sim_unif[, 2], mean = 20, sd = 5))


point_plot_intro_density <- contour_sim %>%
  ggplot(aes(x = size, y = color_raw)) +
  geom_density_2d(show.legend = FALSE, color = "grey", linetype = 5, h = 4) +
  scale_x_continuous(limits = range(concept_data$size)*c(0.90, 1.05), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = range(concept_data$color)*c(0.6, 1.07), expand = expansion(mult = c(0, 0))) +
  labs(y = "color") +
  scale_color_identity() +
  scale_fill_identity() +
  
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


cairo_pdf(file = "example/person_plot.pdf", width = 2.4, height = 2.4, onefile = T)
print(point_plot_intro)
print(point_plot_intro_density)
print(point_plot_intro_density + geom_point(data = concept_data, aes(color = color_hexa2, fill = color_hexa, size = size, y = color), 
                                         show.legend = FALSE, alpha = 0.8, shape = 16))
dev.off()

#new simulated data
color_translation <- data.frame(color = 1:10, 
                                color_hexa = c("#95aaff", "#8597eb", "#7585d8", "#6673c5", "#5661b2", "#47509f", "#383f8d", "#282f7b", "#172069", "#001158"),
                                color_hexa2 = colorRampPalette(c("#f46e32", "#8592BC", "#001158"))(10))


set.seed(12345) #Not the right seed for the data example!!
test_sim <- as.data.frame(rbicop(n = 50, family = "gumbel", 0, parameters = c(3)))

simulated_data <- data.frame(color = round(qunif(test_sim[, 1], min = 1, max = 10)),
                             color_raw = qunif(test_sim[, 1], min = 1, max = 10),
                             size = qnorm(test_sim[, 2], mean = 20, sd = 5)) %>% 
  left_join(color_translation)

point_plot_intro_simulated <- contour_sim %>%
  ggplot(aes(x = size, y = color_raw)) +
  geom_density_2d(show.legend = FALSE, color = "grey", linetype = 5, h = 4) +
  scale_x_continuous(limits = range(concept_data$size)*c(0.90, 1.05), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = range(concept_data$color)*c(0.6, 1.07), expand = expansion(mult = c(0, 0))) +
  labs(y = "color") +
  geom_point(data = simulated_data, aes(color = color_hexa2, fill = color_hexa, size = size), 
             show.legend = FALSE, alpha = 0.8, shape = 15) +
  scale_color_identity() +
  scale_fill_identity() +
  
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

cairo_pdf(file = "example/person_plot.pdf", width = 2.8, height = 2, onefile = T)
print(point_plot_intro)
print(point_plot_intro_density)
print(point_plot_intro_simulated)

print(point_plot_intro + labs(x = NULL, y = NULL))
print(point_plot_intro_density + labs(x = NULL, y = NULL))
print(point_plot_intro_simulated + labs(x = NULL, y = NULL))

print(point_plot_intro + labs(x = "covariate A", y = "covariate B"))
print(point_plot_intro_density + labs(x = "covariate A", y = "covariate B"))
print(point_plot_intro_simulated + labs(x = "covariate A", y = "covariate B"))
dev.off()

#Marginal
set.seed(12345) #Not the right seed for the data example!!
test_sim_marg <- as.data.frame(rbicop(n = 50, family = "indep"))

#okayish
simulated_data_marg <- data.frame(color = round(qunif(test_sim_marg[, 1], min = 1, max = 10)),
                             color_raw = qunif(test_sim_marg[, 1], min = 1, max = 10),
                             size = qnorm(test_sim_marg[, 2], mean = 20, sd = 5)) %>% 
  left_join(color_translation)

#worse
simulated_data_marg <- data.frame(color = round(qunif(test_sim_marg[, 1], min = 1, max = 10)),
                                  color_raw = qunif(test_sim_marg[, 1], min = 1, max = 10),
                                  size = qunif(test_sim_marg[, 2], min = 8.3, max = 28)) %>% 
  left_join(color_translation)

point_plot_intro_simulated_marg <- contour_sim %>%
  ggplot(aes(x = size, y = color_raw)) +
  geom_density_2d(show.legend = FALSE, color = "grey", linetype = 5, h = 4) +
  scale_x_continuous(limits = range(concept_data$size)*c(0.90, 1.05), expand = expansion(mult = c(0, 0))) +
  scale_y_continuous(limits = range(concept_data$color)*c(0.6, 1.07), expand = expansion(mult = c(0, 0))) +
  labs(y = "color") +
  geom_point(data = simulated_data_marg, aes(color = color_hexa2, fill = color_hexa, size = size), 
             show.legend = FALSE, alpha = 0.8, shape = 15) +
  scale_color_identity() +
  scale_fill_identity() +
  
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


cairo_pdf(file = "example/person_plot.pdf", width = 2.8, height = 2, onefile = T)
print(point_plot_intro + labs(x = "covariate A", y = "covariate B"))
print(point_plot_intro_density + labs(x = "covariate A", y = "covariate B"))
print(point_plot_intro_simulated + labs(x = "covariate A", y = "covariate B"))
print(point_plot_intro_simulated_marg + labs(x = "covariate A", y = "covariate B"))
print(point_plot_intro_density + geom_point(data = concept_data, aes(color = color_hexa2, fill = color_hexa, size = size, y = color), 
                                            show.legend = FALSE, alpha = 0.8, shape = 16) + labs(x = "covariate A", y = "covariate B"))
dev.off()

#### PK Plots ####
load("results/PK_24h.Rdata")
#Typical PK
plot_typical_PK <- df_pk_summary_24 %>% 
  filter(type == "observed") %>% 
  ggplot(aes(x = time/60)) +
  geom_line(aes(y = median), color = "#8592BC", size = 1.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time (hour)", y = "Vancomycin concentration (mg/L)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())

pdf(file = "presentation/figures/line_typical_pk_Grimsey.pdf", width = 2.7, height = 2.7)
print(plot_typical_PK)
dev.off()

plot_one_PK <- df_pk_summary_24 %>% 
  filter(type == "observed") %>% 
  ggplot(aes(x = time/60)) +
  geom_ribbon(aes(ymin = p_25,ymax = p_75), fill = "grey80") +
  geom_line(aes(y = median), color = "grey36", size = 1.5) +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24), limits = c(0, 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time", y = "Concentration") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

plot_data_PK <- plot_data %>% 
  filter(time %in% c(60, 120, 600, 1400)) %>% 
  filter(id %in% 1:10) %>% 
  filter(Type == "copula") %>% 
  ggplot(aes(y = conc, x = time/60)) +
  stat_summary(fun = mean, geom = "line", color = "grey49", size = 1.5, linetype = 1) +
  geom_point()  +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24), limits = c(0, 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time", y = "Concentration") +
  labs(color = "Weight (kg)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


col_pats <- read.csv("example/simulated_data_conceptual_figure.csv") %>% 
  arrange(size) %>% 
  mutate(id = 1:20)

plot_model_PK_variation <- plot_data %>% 
  filter(id %in% 1:20) %>% 
  left_join(col_pats) %>% 
  filter(Type == "copula") %>% 
  ggplot(aes(y = conc, x = time/60, color = color_hexa2)) +
  geom_ribbon(data = df_pk_summary_24 %>% filter(type == "copula"), 
              aes(ymin = p_low, ymax = p_high, x = time/60), fill = "grey80", inherit.aes = F) +
  geom_line(aes(group = id),size = 1.5)  +
  scale_color_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0)), breaks = c(seq(0, 20, 5), 24), limits = c(0, 24)) + 
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(x = "Time", y = "Concentration") +
  labs(color = "Weight (kg)") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

pdf(file = "presentation/figures/line_pks_conceptual.pdf", width = 3.57, height = 2)
print(plot_one_PK + labs(x = NULL, y = NULL))
print(plot_data_PK + labs(x = NULL, y = NULL))
print(plot_model_PK_variation + labs(x = NULL, y = NULL))

print(plot_one_PK)
print(plot_data_PK)
print(plot_model_PK_variation)
dev.off()






###############
library(VennDiagram)
list("marginal distribution", "multivariate normal distribution", "conditional distribution", "bootstrap", "copula")

dependency <- c("multivariate normal distribution", "conditional distribution", "bootstrap", "copula")
dist_based <- c("marginal distribution", "multivariate normal distribution", "copula")
new_vps <- c("marginal distribution", "multivariate normal distribution", "conditional distribution", "copula")
flexible <- c("marginal distribution", "conditional distribution", "bootstrap", "copula")
list_venn <- list(dependency = dependency, new_vps = new_vps, dist_based = dist_based, flexible = flexible)

venn_diag <- venn.diagram(list_venn, output = T, filename = NULL)

overlaps <- calculate.overlap(list_venn)
overlaps <- rev(overlaps)

posOverlap <- as.numeric (gsub ("a","", (names (overlaps))))
for (i in 1:length(overlaps)){
  pos <- posOverlap [i]
  venn_diag[[pos+8]]$label <- paste(overlaps[[i]], collapse = "\n")
}

grid.newpage()
grid.draw(venn_diag)

###
source("scripts/functions/functions.R")

data_sharing <- c("marginal distribution", "multivariate normal distribution", "copula")
correct_specification <- c("conditional distribution", "bootstrap", "copula")

list_venn <- list(`Data sharing` = data_sharing, `Correct distribution` = correct_specification)


venn_diag <- venn.diagram(list(`Data sharing` = data_sharing, `Correct distribution` = correct_specification), 
                          output = T, filename = NULL, cat.pos = c(0, 0), disable.logging = T, 
                          fill = create_colors(selected = c("light blue", "dark green")),
                          lty = 'blank')

overlaps <- list(a1 = c("- marginal distribution", "- multivariate normal \ndistribution"), 
                 a2 = c("- conditional distribution", "- bootstrap"), a3 = "- copula")

for (i in 1:length(overlaps)){
  venn_diag[[i + 4]]$label <- paste(overlaps[[i]], collapse = "\n")
}

cairo_pdf(file = "example/venn_methods.pdf", width = 5, height = 5, onefile = T)
grid.newpage()
grid.draw(venn_diag)
dev.off()

########################

set.seed(123)
test_sim <- as.data.frame(rbicop(n = 20, family = "gumbel", 0, parameters = c(3)))
test_sim %>% 
  mutate(transformed_x = qnorm(V1),
         transformed_y = qnorm(V2)) %>% 
  ggplot(aes(x = transformed_x, y = V2)) +
  geom_point() +
  theme_bw()



windowsFonts(`Segoe UI Emoji` = windowsFont("Segoe UI Emoji"))

population_plot <- simulated_data %>% left_join(color_translation) %>% 
  bind_cols(as.data.frame(expand.grid(4:1, 5:1))) %>% 
  ggplot(aes(x = Var1, y = Var2, color = color_hexa2, fill = color_hexa, size = size)) +
  geom_point(shape = "\U0001F9CD", show.legend = FALSE) +
  
  scale_color_identity() +
  scale_fill_identity() +
  scale_x_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  scale_y_continuous(expand = expansion(mult = c(0.3, 0.3))) +
  theme_void()

windowsFonts(`Noto Emoji` = windowsFont("Noto Emoji"))

cairo_pdf("example/population_people.pdf", width = 1, height = 1, family = "Segoe UI Emoji")
print(population_plot)
dev.off()


#windowsFonts(`Segoe UI Emoji` = windowsFont("Segoe UI Emoji"))
point_plot_intro <- simulated_data %>% left_join(color_translation) %>% 
  ggplot(aes(x = size, y = color, color = color_hexa2, fill = color_hexa, size = size)) +
  geom_point(show.legend = FALSE, alpha = 0.8) +
  
  scale_color_identity() +
  scale_fill_identity() +
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


cairo_pdf(file = "example/person_plot2.pdf", width = 2.4, height = 2.4)
print(point_plot_intro)
dev.off()


#emoji::emoji("standing")
#windowsFonts()
#gradient:
# rgb(149, 170, 255)
# rgb(133, 151, 235)
# rgb(117, 133, 216)
# rgb(102, 115, 197)
# rgb(86, 97, 178)
# rgb(71, 80, 159)
# rgb(56, 63, 141)
# rgb(40, 47, 123)
# rgb(23, 32, 105)
# rgb(0, 17, 88)



emoji::emoji("bacteria")


"\uFE0F"


#PK curve:

