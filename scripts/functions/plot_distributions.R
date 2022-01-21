
#plot function for comparing simulation with observed data

plot_comparison_distribution_sim_obs <- function(sim_data, obs_data, sim_nr = 1, title = NULL, plot_type = "points") {
  #Combine truth and simulation results
  total_data <- obs_data %>% mutate(type = "observed") %>%
    bind_rows(sim_data %>% mutate(type = "simulated"))

  total_data$type <- factor(total_data$type, levels = c("simulated", "observed"))
  
  # Use only 1 simulation for plotting
  plot_data <- total_data %>%
   # mutate(simulation_nr = as.character(simulation_nr))
    filter(simulation_nr %in% sim_nr | is.na(simulation_nr))

  
  p_age_crea_sim <- ggplot(mapping = aes(x = age, y = CREA)) +
    geom_point(data = plot_data[plot_data$type == "simulated", ], alpha = 0.3, shape = 16, color = "#F46E32") +
    geom_point(data = plot_data[plot_data$type == "observed", ], alpha = 0.8, shape = 16, color = "#5063B9") +
    theme_bw()
  
  p_age_bw_sim <- ggplot(mapping = aes(x = age, y = BW)) +
    geom_point(data = plot_data[plot_data$type == "simulated", ], alpha = 0.3, shape = 16, color = "#F46E32") +
    geom_point(data = plot_data[plot_data$type == "observed", ], alpha = 0.8, shape = 16, color = "#5063B9") +
    theme_bw()
  
  p_bw_crea_sim <- ggplot(mapping = aes(x = BW, y = CREA)) +
    geom_point(data = plot_data[plot_data$type == "simulated", ], alpha = 0.3, shape = 16, color = "#F46E32") +
    geom_point(data = plot_data[plot_data$type == "observed", ], alpha = 0.8, shape = 16, color = "#5063B9") +
    theme_bw()
  
  # p_simulation_observation <- list(geom_point(shape = 16),
  #                                  scale_color_manual(values = c("#F46E32", "#5063B9")),
  #                                  scale_alpha_manual(values = c(0.1, 0.5)),
  #                                  theme_bw())
  # p_age_crea_sim <- plot_data %>% 
  #   ggplot(aes(x = age, y = CREA, color = type, alpha = type)) +
  #   p_simulation_observation
  # p_age_bw_sim <- plot_data %>% 
  #   ggplot(aes(x = age, y = BW, color = type, alpha = type)) +
  #   p_simulation_observation
  # p_bw_crea_sim <- plot_data %>% 
  #   ggplot(aes(x = BW, y = CREA, color = type, alpha = type)) +
  #   p_simulation_observation
  
  legend_plots <- ggpubr::get_legend(p_age_crea_sim)
  remove_legend <- function(plot) {
    plot + theme(legend.position = "none")
  }
  
  p_univ_simulation_observation <- list(geom_boxplot(), 
                                        scale_fill_manual(values = c("#F46E32", "#5063B9")),
                                        theme_bw())
  p_crea_sim <- plot_data %>% 
    ggplot(aes(y = CREA, x = type, fill = type)) +
    p_univ_simulation_observation
  p_age_sim <- plot_data %>% 
    ggplot(aes(y = age, x = type, fill = type)) +
    p_univ_simulation_observation
  p_bw_sim <- plot_data %>% 
    ggplot(aes(y = BW, x = type, fill = type)) +
    p_univ_simulation_observation
  
  
  
  
  #with densities
  p_age_crea_sim_dens <- ggplot(mapping = aes(x = age, y = CREA)) +
    #geom_point(data = plot_data[plot_data$type == "simulated", ], alpha = 0.1, shape = 16, color = "#F46E32") +
    #geom_point(data = plot_data[plot_data$type == "observed", ], shape = 16, color = "#5063B9") +
    geom_density2d(data = plot_data[plot_data$type == "simulated", ], color = "#F46E32") +
    geom_density2d(data = plot_data[plot_data$type == "observed", ],  alpha = 1, color = "#5063B9", linetype = 2) +
    theme_bw()
  p_age_bw_sim_dens <- ggplot(mapping = aes(x = age, y = BW)) +
    #geom_point(data = plot_data[plot_data$type == "simulated", ], alpha = 0.1, shape = 16, color = "#F46E32") +
    #geom_point(data = plot_data[plot_data$type == "observed", ], shape = 16, color = "#5063B9") +
    geom_density2d(data = plot_data[plot_data$type == "simulated", ], color = "#F46E32") +
    geom_density2d(data = plot_data[plot_data$type == "observed", ], alpha = 1, color = "#5063B9", linetype = 2) +
    theme_bw()
  
  p_bw_crea_sim_dens <- ggplot(mapping = aes(x = BW, y = CREA)) +
    #geom_point(data = plot_data[plot_data$type == "simulated", ], alpha = 0.1, shape = 16, color = "#F46E32") +
    #geom_point(data = plot_data[plot_data$type == "observed", ], shape = 16, color = "#5063B9") +
    geom_density2d(data = plot_data[plot_data$type == "simulated", ], color = "#F46E32") +
    geom_density2d(data = plot_data[plot_data$type == "observed", ], alpha = 1, color = "#5063B9", linetype = 2) +
    theme_bw()
  
  if (plot_type == "points") {
  gridExtra::grid.arrange(remove_legend(p_age_sim),
                          remove_legend(p_age_bw_sim), 
                          remove_legend(p_age_crea_sim), 
                          remove_legend(p_bw_sim),
                          remove_legend(p_bw_crea_sim),
                          remove_legend(p_crea_sim),
                          ggpubr::as_ggplot(legend_plots),
                          layout_matrix = matrix(c(1, 2, 3, NA, 4, 5, 7, NA, 6), ncol = 3), top = title)
  } else if (plot_type == "density") {
  gridExtra::grid.arrange(remove_legend(p_age_sim),
                          remove_legend(p_age_bw_sim_dens), 
                          remove_legend(p_age_crea_sim_dens), 
                          remove_legend(p_bw_sim),
                          remove_legend(p_bw_crea_sim_dens),
                          remove_legend(p_crea_sim),
                          ggpubr::as_ggplot(legend_plots),
                          layout_matrix = matrix(c(1, 2, 3, NA, 4, 5, 7, NA, 6), ncol = 3), top = title)
  } else if (plot_type == "both") {
    gridExtra::grid.arrange(remove_legend(p_age_sim),
                            remove_legend(p_age_bw_sim_dens), 
                            remove_legend(p_age_crea_sim_dens), 
                            remove_legend(p_bw_sim),
                            remove_legend(p_bw_crea_sim_dens),
                            remove_legend(p_crea_sim),
                            remove_legend(p_age_crea_sim),
                            remove_legend(p_age_bw_sim), 
                            remove_legend(p_bw_crea_sim),
                            layout_matrix = matrix(c(1, 2, 3, 8, 4, 5, 7, 9, 6), ncol = 3), top = title)
  }
  
  
}
