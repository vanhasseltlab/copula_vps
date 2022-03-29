#plot function for comparing simulation with observed data

plot_comparison_distribution_sim_obs_generic <- function(sim_data, obs_data, variables = NULL, 
                                                         sim_nr = 1, title = NULL, plot_type = "points", 
                                                         pick_color = c("#F46E32", "#5063B9")) {
  #Combine truth and simulation results
  total_data <- obs_data %>% mutate(type = "observed") %>%
    bind_rows(sim_data %>% mutate(type = "simulated"))
  
  total_data$type <- factor(total_data$type, levels = c("simulated", "observed"))
  
  if ("simulation_nr" %in% colnames(sim_data)) {
    # Use only 1 simulation for plotting
    total_data <- total_data %>%
      filter(simulation_nr %in% sim_nr | is.na(simulation_nr))
  }
  
  if (is.null(variables)) {
    variables <- setdiff(colnames(sim_data), "simulation_nr")
  }
  

  point_plots <- density_plots <- list()
  combination_variables <- t(combn(variables, 2))
  for (i in 1:nrow(combination_variables)) {
    combination <- paste(combination_variables[i, ], collapse = "_")
    point_plots[[combination]] <- ggplot(mapping = aes_string(x = combination_variables[i, 1], y = combination_variables[i, 2])) +
      geom_point(data = total_data[total_data$type == "simulated", ], alpha = 0.3, shape = 16, color = pick_color[1]) +
      geom_point(data = total_data[total_data$type == "observed", ], alpha = 0.8, shape = 16, color = pick_color[2]) +
      theme_bw() + 
      theme(legend.position = "none")
    
    density_plots[[combination]] <- ggplot(mapping = aes_string(x = combination_variables[i, 1], y = combination_variables[i, 2])) +
      geom_density2d(data = total_data[total_data$type == "simulated", ], color = pick_color[1]) +
      geom_density2d(data = total_data[total_data$type == "observed", ],  alpha = 1, color = pick_color[2], linetype = 2) +
      theme_bw() +
      theme(legend.position = "none")
  }
  
  univariate_plots <- list()
  for (i in variables) {
    univariate_plots[[i]] <- total_data %>% 
      ggplot(aes_string(y = i, x = "type", fill = "type")) +
      geom_boxplot() +
      scale_fill_manual(values = c(pick_color[1], pick_color[2])) +
      theme_bw() +
      theme(legend.position = "none")
  }
  nr_cross_plots <- choose(length(variables), 2)
  layout_mat <- matrix(NA, nrow = length(variables), ncol = length(variables))
  layout_mat[lower.tri(layout_mat)] <- 1:nr_cross_plots
  diag(layout_mat) <- nr_cross_plots + (1:length(variables))
 
  
  
  if (plot_type == "points") {
    list_plots <- c(point_plots, univariate_plots)
  } else if (plot_type == "density") {
    list_plots <- c(density_plots, univariate_plots)
  } else if (plot_type == "both") {
    layout_mat[upper.tri(layout_mat)] <- nr_cross_plots + length(variables) + (1:nr_cross_plots)
    list_plots <- c(density_plots, univariate_plots, point_plots)
  }
  
  gridExtra::grid.arrange(grobs = list_plots, layout_matrix = layout_mat, top = title)
}
