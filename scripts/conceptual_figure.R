#Making conceptual figure
library(rvinecopulib)
library(tidyverse)

#cop_distribution
#family: "gumbel", 0, parameters = c(3)
#marginals:   color = round(qunif(x[, 1], min = 1, max = 10)),
#             size  = qnorm(x[, 2], mean = 20, sd = 5)

# example observed data:
concept_data <- read.csv("example/simulated_data_conceptual_figure.csv", row.names = 1)

point_plot_intro <- concept_data %>% left_join(color_translation) %>% 
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
  scale_y_continuous(limits = range(concept_data$color)*c(0.6, 1.07), expand = expansion(mult = c(0, 0)), name = "color") +
  
  geom_point(data = concept_data, aes(color = color_hexa2, fill = color_hexa, size = size, y = color), 
             show.legend = FALSE, alpha = 0.8, shape = 16) +
  scale_color_identity() +
  scale_fill_identity() +
  
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())


cairo_pdf(file = "example/person_plot.pdf", width = 2.4, height = 2.4, onefile = T)
print(point_plot_intro)
print(point_plot_intro_density)
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
  scale_y_continuous(limits = range(concept_data$color)*c(0.6, 1.07), expand = expansion(mult = c(0, 0)), name = "color") +
  
  geom_point(data = simulated_data, aes(color = color_hexa2, fill = color_hexa, size = size), 
             show.legend = FALSE, alpha = 0.8, shape = 15) +
  scale_color_identity() +
  scale_fill_identity() +
  
  theme_bw() +
  theme(axis.title = element_text(size = 12), panel.grid = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

cairo_pdf(file = "example/person_plot.pdf", width = 2.4, height = 2.4, onefile = T)
print(point_plot_intro)
print(point_plot_intro_density)
print(point_plot_intro_simulated)
dev.off()


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
