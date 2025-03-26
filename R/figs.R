my_ggsave <- function(filename, plot, units = "mm", height = NA, width = NA, dpi = NA, ...) {
  
  # # Save the plot as .tiff
  # ggsave(
  #   filename = paste0(filename, ".tiff"),
  #   plot = plot,
  #   height = height,
  #   width = width,
  #   units = units,
  #   dpi = dpi,
  #   ...
  # )
  
  # Save the plot as .pdf
  ggsave(
    filename = paste0(filename, ".pdf"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    ...
  )
  
    # Save the plot as .png
  ggsave(
    filename = paste0(filename, ".png"),
    plot = plot,
    height = height,
    width = width,
    units = units,
    dpi = dpi,
    ...
  )
  # Return the paths to the saved files
  # return(paste0(filename, c(".tiff", ".pdf", ".png")))
  return(paste0(filename, c(".pdf", ".png")))

}


my_theme <- function(){
  theme_bw() %+replace%
  theme(
    plot.margin = margin(t = 5, r = 10, b = 5, l = 5, unit = "pt"),
    legend.key.size = unit(0.3, "cm"),
    legend.spacing.y = unit(0.15, "cm"),
    legend.text.align = 0,
    legend.text = element_text(size = 8),
    legend.title = element_text(size = 9),
    legend.background = element_blank()
  )
}



#' Plot Latitude and Longitude Data on a World Map
#' Generates a map displaying locations based on latitude and longitude data. Points are overlaid on a world map to visualize their distribution.
#' @param tallo_wd_df0 A data frame containing latitude and longitude columns.
#' @return A ggplot object with the world map and plotted latitude and longitude points.
#' @export

data_map <- function(tallo_wd_df0) {
  lat_long <- tallo_wd_df0 |> 
    dplyr::select(latitude, longitude) |> 
    na.omit()

  world_map <- map_data("world")

  p <- ggplot() +
    geom_polygon(data = world_map, aes(x = long, y = lat, group = group), 
                 fill = "gray70", color = NA) +
    geom_point(data = lat_long, aes(x = longitude, y = latitude), 
               color = "#d7191c", alpha = 0.5, size = 1.5) +
    labs(x = "Longitude", y = "Latitude") +
    my_theme()

  return(p)
}



#' Generate a Custom Curve Plot
#' Creates a plot with custom curves representing different mathematical functions: Power-law, generalized Michaelis-Menten (gMM), and Weibull. Curves are normalized for better comparison.
#' @return A ggplot object with custom curves plotted.
#' @export
generate_custom_curve_plot <- function() {
  # Define the functions with appropriate names
  gmm <- function(x, a2, b2, k1) {
    (a2 * x^b2) / (k1 + x^b2)
  }
  
  weibull <- function(x, a3, b3, k2) {
    a3 * (1 - exp(-b3 * x^k2))
  }
  
  power_law <- function(x, a1, b1) {
    a1 * x^b1
  }
  
  # Custom colors
  color_power_law <- "#FF5733"  # Light Orange color
  color_gmm <- "#1b9e77" 
  color_weibull <- "#0f92e9"    # Light Blue color

  # Create a data frame with the x values and the y values for each curve
  x_values <- seq(0, 5, length.out = 1000)  # Increased points for smoother curves
  y_power_law <- power_law(x_values, a1 = 1, b1 = 0.5)
  y_weibull <- weibull(x_values, a3 = 1.5, b3 = 0.5, k2 = 2)
  y_gmm <- gmm(x_values, a2 = 2, b2 = 1.5, k1 = 1)
  
  # Normalize the y-values to make the upper limits similar
  y_max <- max(c(y_power_law, y_weibull, y_gmm))
  y_power_law <- y_power_law / max(y_power_law) * y_max
  y_weibull <- y_weibull / max(y_weibull) * y_max
  y_gmm <- y_gmm / max(y_gmm) * y_max
  
  # Create the data frame
  data <- data.frame(
    x = rep(x_values, 3),
    y = c(y_power_law, y_weibull, y_gmm),
    curve_type = rep(c("Power-law", "Weibull", "gMM"), each = length(x_values))
  )
  
  # Plotting the curves with default grid lines
  p <- ggplot(data, aes(x = x, y = y, color = curve_type)) +
    geom_line(size = 0.4) +
    scale_color_manual(values = c("Power-law" = color_power_law, "Weibull" = color_weibull, "gMM" = color_gmm)) +
    labs(x = "DBH", y = "Tree Height/ Crown Radius", color = "Model") +  # Updated axis labels and legend title
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),           # Remove X-axis text
      axis.text.y = element_blank()            # Remove Y-axis text
    ) +
    my_theme()

  return(p)
}


#====================
# Allometry plot
#====================
generate_allo_plot <- function(
  tallo_reduced_lr_df_ang_h, tallo_reduced_lr_df_ang_cr,
  tallo_reduced_lr_df_gym_h, tallo_reduced_lr_df_gym_cr,
  fit_lr_nou_summary_pl_ang_h, fit_lr_nou_summary_pl_gym_h,
  fit_nlr_nou_summary_gmm_ang_h, fit_nlr_nou_summary_gmm_gym_h,
  fit_nlr_nou_summary_weibull_ang_h, fit_nlr_nou_summary_weibull_gym_h,
  fit_lr_nou_summary_pl_ang_cr, fit_lr_nou_summary_pl_gym_cr,
  fit_nlr_nou_summary_gmm_ang_cr, fit_nlr_nou_summary_gmm_gym_cr,
  fit_nlr_nou_summary_weibull_ang_cr, fit_nlr_nou_summary_weibull_gym_cr,
  log_scale = TRUE
) {

  # Define custom colors for the curve types
  color_power_law <- "#FF5733"
  # color_gmm <- "#31a354"
  color_gmm <- "#1b9e77" 
  color_weibull <- "#0f92e9"

  # Separate the datasets for Angiosperm and Gymnosperm
  ang_data_h <- tallo_reduced_lr_df_ang_h |> filter(division == "Angiosperm")
  gym_data_h <- tallo_reduced_lr_df_gym_h |> filter(division == "Gymnosperm")
  
  ang_data_cr <- tallo_reduced_lr_df_ang_cr |> filter(division == "Angiosperm")
  gym_data_cr <- tallo_reduced_lr_df_gym_cr |> filter(division == "Gymnosperm")

  # Define function to compute model predictions
  compute_fitted_values <- function(model_summary, x_seq, model_type) {
    gamma <- model_summary |> filter(grepl("gamma", variable))
    
    if (model_type == "pl") {
      log_a <- gamma |> filter(variable == "gamma[1]") |> pull(q50)
      a <- exp(log_a)  # Only exponentiate log_a
      b <- gamma |> filter(variable == "gamma[2]") |> pull(q50)
      log_y_pred <- log(a) + b * log(x_seq)
      fitted_values <- tibble(y = exp(log_y_pred), x = x_seq)
      
    } else if (model_type == "gmm") {
      log_a <- gamma |> filter(variable == "gamma[1]") |> pull(q50)
      a <- exp(log_a)
      b <- gamma |> filter(variable == "gamma[2]") |> pull(q50)
      k <- gamma |> filter(variable == "gamma[3]") |> pull(q50)
      
      log_y_pred <- log(a) + b * log(x_seq) - log(k + x_seq^b)
      fitted_values <- tibble(y = exp(log_y_pred), x = x_seq)
      
    } else if (model_type == "weibull") {
      log_a <- gamma |> filter(variable == "gamma[1]") |> pull(q50)
      a <- exp(log_a)
      b <- gamma |> filter(variable == "gamma[2]") |> pull(q50)
      k <- gamma |> filter(variable == "gamma[3]") |> pull(q50)
      
      log_y_pred <- log(a) + log(1 - exp(-b * (x_seq ^ k)))
      fitted_values <- tibble(y = exp(log_y_pred), x = x_seq)
      
    } else {
      stop("Invalid model type.")
    }
    
    return(fitted_values)
  }

  # Create prediction sequences for DBH
  x_seq_ang_h <- seq(min(ang_data_h$dbh), max(ang_data_h$dbh), length.out = 100)
  x_seq_gym_h <- seq(min(gym_data_h$dbh), max(gym_data_h$dbh), length.out = 100)
  
  x_seq_ang_cr <- seq(min(ang_data_cr$dbh), max(ang_data_cr$dbh), length.out = 100)
  x_seq_gym_cr <- seq(min(gym_data_cr$dbh), max(gym_data_cr$dbh), length.out = 100)
  
  # Compute fitted values for H-DBH models
  df_pl_ang_h <- compute_fitted_values(fit_lr_nou_summary_pl_ang_h, x_seq_ang_h, "pl")
  df_pl_gym_h <- compute_fitted_values(fit_lr_nou_summary_pl_gym_h, x_seq_gym_h, "pl")
  
  df_gmm_ang_h <- compute_fitted_values(fit_nlr_nou_summary_gmm_ang_h, x_seq_ang_h, "gmm")
  df_gmm_gym_h <- compute_fitted_values(fit_nlr_nou_summary_gmm_gym_h, x_seq_gym_h, "gmm")
  
  df_wb_ang_h <- compute_fitted_values(fit_nlr_nou_summary_weibull_ang_h, x_seq_ang_h, "weibull")
  df_wb_gym_h <- compute_fitted_values(fit_nlr_nou_summary_weibull_gym_h, x_seq_gym_h, "weibull")

  # Compute fitted values for CR-DBH models
  df_pl_ang_cr <- compute_fitted_values(fit_lr_nou_summary_pl_ang_cr, x_seq_ang_cr, "pl")
  df_pl_gym_cr <- compute_fitted_values(fit_lr_nou_summary_pl_gym_cr, x_seq_gym_cr, "pl")
  
  df_gmm_ang_cr <- compute_fitted_values(fit_nlr_nou_summary_gmm_ang_cr, x_seq_ang_cr, "gmm")
  df_gmm_gym_cr <- compute_fitted_values(fit_nlr_nou_summary_gmm_gym_cr, x_seq_gym_cr, "gmm")
  
  df_wb_ang_cr <- compute_fitted_values(fit_nlr_nou_summary_weibull_ang_cr, x_seq_ang_cr, "weibull")
  df_wb_gym_cr <- compute_fitted_values(fit_nlr_nou_summary_weibull_gym_cr, x_seq_gym_cr, "weibull")

  # Plotting function for density-based comparison plots
  plot_density <- function(data, df_pl, df_gmm, df_wb, x_var, y_var, x_label, y_label, log_scale) {
    if (log_scale) {
      dens <- with(data, kde2d(log10(data[[x_var]]), log10(data[[y_var]]), n = 300))
      ix <- findInterval(log10(data[[x_var]]), dens$x)
      iy <- findInterval(log10(data[[y_var]]), dens$y)
    } else {
      dens <- with(data, kde2d(data[[x_var]], data[[y_var]], n = 300))
      ix <- findInterval(data[[x_var]], dens$x)
      iy <- findInterval(data[[y_var]], dens$y)
    }
    
    data$density <- dens$z[cbind(ix, iy)]
    
    p <- ggplot(data, aes(x = !!sym(x_var), y = !!sym(y_var))) +
      geom_point(aes(color = density), size = 0.1, alpha = 0.3) +
      scale_color_gradientn(colors = c("black", "purple", "orange", "yellow"), name = "Density",
      guide = guide_colourbar(barheight = unit(1.2, "cm"))) +
      geom_line(data = df_pl, aes(x = x, y = y), color = color_power_law, size = 0.5, show.legend = FALSE) +
      geom_line(data = df_gmm, aes(x = x, y = y), color = color_gmm, size = 0.5, show.legend = FALSE) +
      geom_line(data = df_wb, aes(x = x, y = y), color = color_weibull, size = 0.5, show.legend = FALSE) +
      labs(x = x_label, y = y_label) +
      my_theme() +
      theme(
        # legend.position = c(0.9, 0.2),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6)
      )
    if (log_scale) {
      p <- p + scale_x_log10() + scale_y_log10() +
        theme(legend.position = c(0.9, 0.2))  # Default for log-scale
    } else {
      p <- p + theme(legend.position = c(0.85, 0.2))  # Move the legend left for non-log-scale plot
    }
    

    if (log_scale) {
      p <- p + scale_x_log10() + scale_y_log10()
    }
    
    return(p)
  }

  # Generate individual plots with log_scale
  p_ang_h <- plot_density(ang_data_h, df_pl_ang_h, df_gmm_ang_h, df_wb_ang_h, "dbh", "h", "DBH (cm)", "Height (m)", log_scale)
  p_ang_cr <- plot_density(ang_data_cr, df_pl_ang_cr, df_gmm_ang_cr, df_wb_ang_cr, "dbh", "cr", "DBH (cm)", "Crown Radius (m)", log_scale)
  p_gym_h <- plot_density(gym_data_h, df_pl_gym_h, df_gmm_gym_h, df_wb_gym_h, "dbh", "h", "DBH (cm)", "Height (m)", log_scale)
  p_gym_cr <- plot_density(gym_data_cr, df_pl_gym_cr, df_gmm_gym_cr, df_wb_gym_cr, "dbh", "cr", "DBH (cm)", "Crown Radius (m)", log_scale)

  # Combine the plots into a grid
  combined_plot <- (p_ang_h | p_gym_h | p_ang_cr | p_gym_cr) +
    plot_layout(ncol = 2, nrow = 2) +
    plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)")))

  # Extract the Curve Type legend separately
  legend_curve_type <- get_legend(
    ggplot() + 
      geom_line(aes(x = 1:10, y = 1:10, linetype = "Power-law"), color = color_power_law, size = 0.5) +
      geom_line(aes(x = 1:10, y = 1:10, linetype = "gMM"), color = color_gmm, size = 0.5) +
      geom_line(aes(x = 1:10, y = 1:10, linetype = "Weibull"), color = color_weibull, size = 0.5) +
      scale_linetype_manual(values = c("Power-law" = "solid", "gMM" = "solid", "Weibull" = "solid"), 
                            name = "Model") +
      my_theme()
  )

  # Use cowplot to add the Curve Type legend to the right side
  p <- plot_grid(combined_plot, legend_curve_type, rel_widths = c(0.6, 0.1))

  # Add labels for Angiosperm and Gymnosperm sections
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
    draw_label("Angiosperms", x = 0.26, y = 0.93, fontface = 'bold', size = 10) +  
    draw_label("Gymnosperms", x = 0.685, y = 0.93, fontface = 'bold', size = 10)

  return(p)
}

#=============================
# DBH-CR and/or H allometry
#=============================



#====================
# WD vs. H-DBH
#====================
generate_wd_para_h <- function(tallo_reduced_nlr_df_ang_h, tallo_reduced_nlr_df_gym_h, fit_wd_ang_h_draws_weibull_wd, fit_wd_gym_h_draws_weibull_wd, stan_data_nlr_ang_h, stan_data_nlr_gym_h) {

# Filter Angiosperms and Gymnosperms
  ang_h_df <- tallo_reduced_nlr_df_ang_h |> filter(division == "Angiosperm")
  gym_h_df <- tallo_reduced_nlr_df_gym_h |> filter(division == "Gymnosperm")

  # Group by species and calculate mean wood density for Angiosperms
  wd_df_ang <- ang_h_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Group by species and calculate mean wood density for Gymnosperms
  wd_df_gym <- gym_h_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Create the sequence for wood density for Angiosperms and Gymnosperms separately
  wd_seq_ang <- seq(min(wd_df_ang$wd_s), max(wd_df_ang$wd_s), length.out = 100)
  wd_seq_gym <- seq(min(wd_df_gym$wd_s), max(wd_df_gym$wd_s), length.out = 100)

  wd_seq_df_ang <- data.frame(
    wd_s = wd_seq_ang,
    wd = approx(wd_df_ang$wd_s, wd_df_ang$wd, xout = wd_seq_ang)$y
  )

  wd_seq_df_gym <- data.frame(
    wd_s = wd_seq_gym,
    wd = approx(wd_df_gym$wd_s, wd_df_gym$wd, xout = wd_seq_gym)$y
  )

  # Function to prepare and plot parameters (log_a, b, k) for H-DBH of Ang and Gym
  prepare_and_plot_parameter <- function(stan_data, fit_draws, wd_df, param_name, color_point, color_ribbon, line_color, line_type = "solid", show_ribbon = TRUE) {
    # Prepare wood density data
    wd_df <- wd_df |>
      filter(!is.na(wd_s)) |>
      group_by(sp) |>
      summarize(wd = mean(wd, na.rm = TRUE)) |>
      mutate(wd_s = scale(wd) |> as.numeric(), sp_id = 1:n())

    # Extract gamma and beta parameters from fit_draws
    gamma_params <- fit_draws |>
      dplyr::select(starts_with("gamma"))

    # Reshape the beta parameters for analysis
    beta_params <- fit_draws |>
      dplyr::select(starts_with("beta[")) |>
      pivot_longer(
        cols = everything(),
        names_to = c("species", "parameter"),
        names_pattern = "beta\\[(\\d+),(\\d+)\\]",
        names_transform = list(species = as.integer, parameter = as.integer),
        values_to = "value"
      ) |>
      mutate(parameter = case_when(
        parameter == 1 ~ "log_a",
        parameter == 2 ~ "b",
        parameter == 3 ~ "k"
      ))

    # Join wood density data with beta parameters
    plot_data <- beta_params |>
      left_join(wd_df, by = c("species" = "sp_id"))

    # Summarize the posterior distributions for plotting
    summary_stats <- plot_data |>
      group_by(species, parameter, wd, wd_s) |>
      summarize(
        median = median(value, na.rm = TRUE),
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    # Create a sequence for wd_s and merge with corresponding wd values
    wd_seq <- seq(min(wd_df$wd_s), max(wd_df$wd_s), length.out = 100)
    wd_seq_df <- data.frame(
      wd_s = wd_seq,
      wd = approx(wd_df$wd_s, wd_df$wd, xout = wd_seq)$y
    )

    # Extract relevant gamma parameters for the chosen model
    gamma_param <- gamma_params |>
      dplyr::select(paste0("gamma[1,", which(c("log_a", "b", "k") == param_name), "]"), 
             paste0("gamma[2,", which(c("log_a", "b", "k") == param_name), "]"))

    # Generate lines and fit data for the selected parameter
    transform_value <- if (param_name == "log_a") {
      function(x) x
    } else {
      function(x) x
    }
    
    lines_data <- do.call(rbind, lapply(1:nrow(gamma_param), function(i) {
      data.frame(
        wd_s = wd_seq_df$wd_s,
        wd = wd_seq_df$wd,
        value = transform_value(gamma_param[[1]][i] + gamma_param[[2]][i] * wd_seq_df$wd_s),
        draw = i
      )
    }))
    
    # Compute the fit line and 95% credible intervals
    median_intercept <- median(gamma_param[[1]], na.rm = TRUE)
    median_slope <- median(gamma_param[[2]], na.rm = TRUE)

    fit_line <- data.frame(
      wd_s = wd_seq_df$wd_s,
      wd = wd_seq_df$wd,
      value = transform_value(median_intercept + median_slope * wd_seq_df$wd_s)
    )
    
    ci_ribbon <- lines_data |>
      group_by(wd) |>
      summarize(
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    fit_with_ci <- left_join(fit_line, ci_ribbon, by = "wd")

    # Apply transformation for log_a in summary_stats
    if (param_name == "log_a") {
      summary_stats <- summary_stats |>
        mutate(median = median, lower = lower, upper = upper)
    }
    
    # Set y-axis label based on parameter
    y_label <- case_when(
      param_name == "log_a" ~ "Asymptote parameter, a",
      param_name == "b" ~ "Scale parameter, b",
      param_name == "k" ~ "Shape parameter, k"
    )

    # Create plot
    plot <- ggplot() +
      geom_errorbar(data = summary_stats |> filter(parameter == param_name), 
                    aes(x = wd, ymin = lower, ymax = upper), 
                    width = 0, linewidth = 0.2, color = color_ribbon) +
      geom_point(data = summary_stats |> filter(parameter == param_name), 
                 aes(x = wd, y = median), color = color_point, size = 0.5) +
      geom_line(data = fit_line, aes(x = wd, y = value), 
                color = line_color, linewidth = 0.6, linetype = line_type) +
      theme_minimal() +
      labs(x = expression("Wood Density (g cm"^"-3"~")"), y = y_label) +
      my_theme()

    if (show_ribbon) {
      plot <- plot + geom_ribbon(data = fit_with_ci, aes(x = wd, ymin = lower, ymax = upper), alpha = 0.5, fill = color_ribbon)
    }

    return(plot)
  }

  # Generate individual plots for both H-DBH of Ang and Gym
  p1 <- prepare_and_plot_parameter(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE)
  p2 <- prepare_and_plot_parameter(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "b", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)
  p3 <- prepare_and_plot_parameter(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "k", "#6baed6", "#8fcbee", "#2171b5", show_ribbon = TRUE)
  p4 <- prepare_and_plot_parameter(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE)
  p5 <- prepare_and_plot_parameter(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "b", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)
  p6 <- prepare_and_plot_parameter(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "k", "#6baed6", "#8fcbee", "#2171b5", "dashed", show_ribbon = FALSE)

  # Combine all plots
  p <- (p1 + labs(tag = "(a)", x = expression("Wood Density (g cm"^"-3"~")")) |
                    p4 + labs(tag = "(b)", x = expression("Wood Density (g cm"^"-3"~")")) |
                    p2 + labs(tag = "(c)", x = expression("Wood Density (g cm"^"-3"~")")) |
                    p5 + labs(tag = "(d)", x = expression("Wood Density (g cm"^"-3"~")")) |
                    p3 + labs(tag = "(e)", x = expression("Wood Density (g cm"^"-3"~")")) |
                    p6 + labs(tag = "(f)", x = expression("Wood Density (g cm"^"-3"~")"))) +
    plot_layout(ncol = 2, nrow = 3) +
    plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))) &
    my_theme()

  # Add labels
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.9) +  
    draw_label("WD vs. H-DBH: Angiosperm", x = 0.31, y = 0.91, fontface = 'bold', size = 10) +  
    draw_label("WD vs. H-DBH: Gymnosperm", x = 0.785, y = 0.91, fontface = 'bold', size = 10)

  return(p)
}


#====================
# WD vs. CR-DBH
#====================
generate_wd_para_cr <- function(
  tallo_reduced_lr_df_ang_cr, 
  tallo_reduced_nlr_df_gym_cr, 
  fit_wd_ang_cr_draws_pl_wd, 
  fit_wd_gym_cr_draws_gmm_wd, 
  stan_data_lr_ang_cr, 
  stan_data_nlr_gym_cr
) {
  # Helper function to prepare and plot parameters
  prepare_and_plot_parameter <- function(
    stan_data, 
    fit_draws, 
    wd_df, 
    wd_seq_df, 
    param_name, 
    color_point, 
    color_ribbon, 
    line_color, 
    line_type = "solid", 
    show_ribbon = TRUE, 
    is_angiosperm = TRUE
  ) {
    # Extract gamma parameters from fit_draws
    gamma_param <- fit_draws |> dplyr::select(starts_with("gamma"))

    if (is_angiosperm) {
      if (param_name == "log_a") {
        gamma_intercept <- gamma_param[["gamma[1,1]"]]
        gamma_slope <- gamma_param[["gamma[1,2]"]]
      } else if (param_name == "b") {
        gamma_intercept <- gamma_param[["gamma[2,1]"]]
        gamma_slope <- gamma_param[["gamma[2,2]"]]
      }
    } else {
      if (param_name == "log_a") {
        gamma_intercept <- gamma_param[["gamma[1,1]"]]
        gamma_slope <- gamma_param[["gamma[2,1]"]]
      } else if (param_name == "b") {
        gamma_intercept <- gamma_param[["gamma[1,2]"]]
        gamma_slope <- gamma_param[["gamma[2,2]"]]
      } else if (param_name == "k") {
        gamma_intercept <- gamma_param[["gamma[1,3]"]]
        gamma_slope <- gamma_param[["gamma[2,3]"]]
      }
    }

    # Reshape the beta parameters for analysis
    if (is_angiosperm) {
      # Angiosperm beta parameters (log_a, b only)
      beta_params <- fit_draws |>
        dplyr::select(starts_with("beta[")) |>
        pivot_longer(
          cols = everything(),
          names_to = c("parameter", "species"),
          names_pattern = "beta\\[(\\d+),(\\d+)\\]",
          names_transform = list(parameter = as.integer, species = as.integer),
          values_to = "value"
        ) |>
        mutate(parameter = case_when(
          parameter == 1 ~ "log_a",
          parameter == 2 ~ "b"
        ))
    } else {
      # Gymnosperm beta parameters (log_a, b, k)
      beta_params <- fit_draws |>
        dplyr::select(starts_with("beta[")) |>
        pivot_longer(
          cols = everything(),
          names_to = c("species", "parameter"),
          names_pattern = "beta\\[(\\d+),(\\d+)\\]",
          names_transform = list(species = as.integer, parameter = as.integer),
          values_to = "value"
        ) |>
        mutate(parameter = case_when(
          parameter == 1 ~ "log_a",
          parameter == 2 ~ "b",
          parameter == 3 ~ "k"
        ))
    }

    # Filter the beta parameters for the specific parameter (log_a, b, or k)
    beta_filtered <- beta_params |> filter(parameter == param_name)

    # Combine species-specific beta with wd_df
    plot_data <- beta_filtered |> left_join(wd_df, by = c("species" = "sp_id"))

    # Summarize the posterior distributions for plotting
    summary_stats <- plot_data |>
      group_by(species, parameter, wd, wd_s) |>
      summarize(
        median = median(value, na.rm = TRUE),
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    # Prepare line and CI
    lines_data <- do.call(rbind, lapply(1:length(gamma_intercept), function(i) {
      data.frame(
        wd_s = wd_seq_df$wd_s,
        wd = wd_seq_df$wd,
        value = if (param_name == "log_a") {
          gamma_intercept[i] + gamma_slope[i] * wd_seq_df$wd_s
        } else {
          gamma_intercept[i] + gamma_slope[i] * wd_seq_df$wd_s
        },
        draw = i
      )
    }))

    fit_line <- data.frame(
      wd_s = wd_seq_df$wd_s,
      wd = wd_seq_df$wd,
      value = median(gamma_intercept) + median(gamma_slope) * wd_seq_df$wd_s
    )

    ci_ribbon <- lines_data |>
      group_by(wd) |>
      summarize(
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    fit_with_ci <- left_join(fit_line, ci_ribbon, by = "wd")
   
    if (param_name == "log_a") {
      summary_stats <- summary_stats |>
        mutate(median = median, lower = lower, upper = upper)
    }

    # Define y-labels for Angiosperms and Gymnosperms
    y_label <- if (is_angiosperm) {
      ifelse(param_name == "log_a", "Scaling factor, a", "Exponent parameter, b")
    } else {
      case_when(
        param_name == "log_a" ~ "Asymptote parameter, a",
        param_name == "b" ~ "Exponent parameter, b",
        param_name == "k" ~ "Scaling constant, k"
      )
    }
    
    # Create plot
    plot <- ggplot() +
      geom_errorbar(data = summary_stats, aes(x = wd, ymin = lower, ymax = upper), 
                    width = 0, linewidth = 0.2, color = color_ribbon) +
      geom_point(data = summary_stats, aes(x = wd, y = median), color = color_point, size = 0.5) +
      geom_line(data = fit_line, aes(x = wd, y = value), color = line_color, linewidth = 0.6, linetype = line_type) +
      theme_minimal() +
      my_theme() +
      labs(x = "Wood Density (g cm⁻³)", y = y_label)

    # Conditionally show ribbon if line_type is not "dashed"
    if (show_ribbon && line_type != "dashed") {
      plot <- plot + geom_ribbon(data = fit_with_ci, aes(x = wd, ymin = lower, ymax = upper), alpha = 0.5, fill = color_ribbon)
    }

    return(plot)
  }

  # 1. ANGIOSPERMS SECTION

  # Filter Angiosperms
  ang_cr_df <- tallo_reduced_lr_df_ang_cr |> filter(division == "Angiosperm")

  # Group by species and calculate mean wood density for Angiosperms
  wd_df_ang <- ang_cr_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Create the sequence for wood density for Angiosperms
  wd_seq_ang <- seq(min(wd_df_ang$wd_s), max(wd_df_ang$wd_s), length.out = 100)

  wd_seq_df_ang <- data.frame(
    wd_s = wd_seq_ang,
    wd = approx(wd_df_ang$wd_s, wd_df_ang$wd, xout = wd_seq_ang)$y
  )

  # Generate plots for Angiosperms
  p1_ang_log_a <- prepare_and_plot_parameter(stan_data_lr_ang_cr, fit_wd_ang_cr_draws_pl_wd, wd_df_ang, wd_seq_df_ang, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", is_angiosperm = TRUE)
  p2_ang_b <- prepare_and_plot_parameter(stan_data_lr_ang_cr, fit_wd_ang_cr_draws_pl_wd, wd_df_ang, wd_seq_df_ang, "b", "#74c476", "#bae4b3", "#238b45", is_angiosperm = TRUE)

  # 2. GYMNOSPERMS SECTION

  # Filter Gymnosperms
  gym_cr_df <- tallo_reduced_nlr_df_gym_cr |> filter(division == "Gymnosperm")

  # Group by species and calculate mean wood density for Gymnosperms
  wd_df_gym <- gym_cr_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Create the sequence for wood density for Gymnosperms
  wd_seq_gym <- seq(min(wd_df_gym$wd_s), max(wd_df_gym$wd_s), length.out = 100)

  wd_seq_df_gym <- data.frame(
    wd_s = wd_seq_gym,
    wd = approx(wd_df_gym$wd_s, wd_df_gym$wd, xout = wd_seq_gym)$y
  )

  # Generate plots for Gymnosperm
  p1_gym_log_a <- prepare_and_plot_parameter(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym, wd_seq_df_gym, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE, is_angiosperm = FALSE)
  p2_gym_b <- prepare_and_plot_parameter(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym, wd_seq_df_gym, "b", "#74c476", "#bae4b3", "#238b45", show_ribbon = TRUE, is_angiosperm = FALSE)
  p3_gym_k <- prepare_and_plot_parameter(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym, wd_seq_df_gym, "k", "#6baed6", "#8fcbee", "#2171b5", "dashed", show_ribbon = FALSE, is_angiosperm = FALSE)

  # Create a blank placeholder plot
  p_placeholder <- ggplot() + theme_void() + theme(plot.margin = margin(0, 0, 0, 0))

  # Combine all plots
  p <- (p1_ang_log_a + labs(tag = "(a)", x = "Wood Density (g cm"^"-3"~")") |
        p1_gym_log_a + labs(tag = "(c)", x = "Wood Density (g cm"^"-3"~")") |
        p2_ang_b + labs(tag = "(b)", x = "Wood Density (g cm"^"-3"~")") |
        p2_gym_b + labs(tag = "(d)", x = "Wood Density (g cm"^"-3"~")") |
        plot_spacer() |  # Leave empty space without a tag
        p3_gym_k + labs(tag = "(e)", x = "Wood Density (g cm"^"-3"~")")) +
    plot_layout(ncol = 2, nrow = 3, heights = c(1, 1, 1)) + 
    plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)"))) & 
    my_theme()

  # Add labels
  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 0.9) +  
    draw_label("WD vs. CR-DBH: Angiosperms", x = 0.295, y = 0.91, fontface = 'bold', size = 10) +  
    draw_label("WD vs. CR-DBH: Gymnosperms", x = 0.78, y = 0.91, fontface = 'bold', size = 10)
  # p <- ggdraw() +
  #   draw_plot(p, x = 0, y = 0, width = 1, height = 0.9) +  
  #   draw_label("Angiosperm", x = 0.295, y = 0.91, fontface = 'bold', size = 10) +  
  #   draw_label("Gymnosperm", x = 0.78, y = 0.91, fontface = 'bold', size = 10)

  return(p)
}



# #============================
# # WOOD DENSITY-PARA COMBINED
# #============================
generate_wd_para_com <- function(
  tallo_reduced_nlr_df_ang_h, tallo_reduced_nlr_df_gym_h, 
  fit_wd_ang_h_draws_weibull_wd, fit_wd_gym_h_draws_weibull_wd, 
  stan_data_nlr_ang_h, stan_data_nlr_gym_h,
  tallo_reduced_lr_df_ang_cr, tallo_reduced_nlr_df_gym_cr, 
  fit_wd_ang_cr_draws_pl_wd, fit_wd_gym_cr_draws_gmm_wd, 
  stan_data_lr_ang_cr, stan_data_nlr_gym_cr
) {
  #====================
  # WD vs. H-DBH Section
  #====================
  
  # Filter Angiosperms and Gymnosperms
  ang_h_df <- tallo_reduced_nlr_df_ang_h |> filter(division == "Angiosperm")
  gym_h_df <- tallo_reduced_nlr_df_gym_h |> filter(division == "Gymnosperm")

  # Group by species and calculate mean wood density for Angiosperms
  wd_df_ang <- ang_h_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Group by species and calculate mean wood density for Gymnosperms
  wd_df_gym <- gym_h_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Create the sequence for wood density for Angiosperms and Gymnosperms separately
  wd_seq_ang <- seq(min(wd_df_ang$wd_s), max(wd_df_ang$wd_s), length.out = 100)
  wd_seq_gym <- seq(min(wd_df_gym$wd_s), max(wd_df_gym$wd_s), length.out = 100)

  wd_seq_df_ang <- data.frame(
    wd_s = wd_seq_ang,
    wd = approx(wd_df_ang$wd_s, wd_df_ang$wd, xout = wd_seq_ang)$y
  )

  wd_seq_df_gym <- data.frame(
    wd_s = wd_seq_gym,
    wd = approx(wd_df_gym$wd_s, wd_df_gym$wd, xout = wd_seq_gym)$y
  )

  # Function to prepare and plot parameters (log_a, b, k) for H-DBH of Ang and Gym
  prepare_and_plot_parameter_h <- function(
    stan_data,
    fit_draws,
    wd_df,
    param_name,
    color_point,
    color_ribbon,
    line_color,
    line_type = "solid",
    show_ribbon = TRUE
  ) {
    # Prepare wood density data
    wd_df <- wd_df |>
      filter(!is.na(wd_s)) |>
      group_by(sp) |>
      summarize(wd = mean(wd, na.rm = TRUE)) |>
      mutate(wd_s = scale(wd) |> as.numeric(), sp_id = 1:n())

    # Extract gamma and beta parameters from fit_draws
    gamma_params <- fit_draws |>
      dplyr::select(starts_with("gamma"))

    # Reshape the beta parameters for analysis
    beta_params <- fit_draws |>
      dplyr::select(starts_with("beta[")) |>
      pivot_longer(
        cols = everything(),
        names_to = c("species", "parameter"),
        names_pattern = "beta\\[(\\d+),(\\d+)\\]",
        names_transform = list(species = as.integer, parameter = as.integer),
        values_to = "value"
      ) |>
      mutate(parameter = case_when(
        parameter == 1 ~ "log_a",
        parameter == 2 ~ "b",
        parameter == 3 ~ "k"
      ))

    # Join wood density data with beta parameters
    plot_data <- beta_params |>
      left_join(wd_df, by = c("species" = "sp_id"))

    # Summarize the posterior distributions for plotting
    summary_stats <- plot_data |>
      group_by(species, parameter, wd, wd_s) |>
      summarize(
        median = median(value, na.rm = TRUE),
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    # Create a sequence for wd_s and merge with corresponding wd values
    wd_seq <- seq(min(wd_df$wd_s), max(wd_df$wd_s), length.out = 100)
    wd_seq_df <- data.frame(
      wd_s = wd_seq,
      wd = approx(wd_df$wd_s, wd_df$wd, xout = wd_seq)$y
    )

    # Extract relevant gamma parameters for the chosen model
    gamma_param <- gamma_params |>
      dplyr::select(paste0("gamma[1,", which(c("log_a", "b", "k") == param_name), "]"), 
             paste0("gamma[2,", which(c("log_a", "b", "k") == param_name), "]"))

    # Generate lines and fit data for the selected parameter
    transform_value <- if (param_name == "log_a") {
      function(x) x
    } else {
      function(x) x
    }
    
    lines_data <- do.call(rbind, lapply(1:nrow(gamma_param), function(i) {
      data.frame(
        wd_s = wd_seq_df$wd_s,
        wd = wd_seq_df$wd,
        value = transform_value(gamma_param[[1]][i] + gamma_param[[2]][i] * wd_seq_df$wd_s),
        draw = i
      )
    }))
    
    # Compute the fit line and 95% credible intervals
    median_intercept <- median(gamma_param[[1]], na.rm = TRUE)
    median_slope <- median(gamma_param[[2]], na.rm = TRUE)

    fit_line <- data.frame(
      wd_s = wd_seq_df$wd_s,
      wd = wd_seq_df$wd,
      value = transform_value(median_intercept + median_slope * wd_seq_df$wd_s)
    )
    
    ci_ribbon <- lines_data |>
      group_by(wd) |>
      summarize(
        lower = quantile(value, 0.025, na.rm = TRUE),
        upper = quantile(value, 0.975, na.rm = TRUE),
        .groups = "drop"
      )

    fit_with_ci <- left_join(fit_line, ci_ribbon, by = "wd")

    # Apply transformation for log_a in summary_stats
    if (param_name == "log_a") {
      summary_stats <- summary_stats |>
        mutate(median = median, lower = lower, upper = upper)
    }
    
    # Set y-axis label based on parameter
    y_label <- case_when(
      param_name == "log_a" ~ "Asymptote parameter, log a",
      param_name == "b" ~ "Scale parameter, b",
      param_name == "k" ~ "Shape parameter, k"
    )

    # Create plot
    plot <- ggplot() +
      geom_errorbar(data = summary_stats |> filter(parameter == param_name), 
                    aes(x = wd, ymin = lower, ymax = upper), 
                    width = 0, linewidth = 0.2, color = color_ribbon) +
      geom_point(data = summary_stats |> filter(parameter == param_name), 
                 aes(x = wd, y = median), color = color_point, size = 0.5) +
      geom_line(data = fit_line, aes(x = wd, y = value), 
                color = line_color, linewidth = 0.6, linetype = line_type) +
      labs(x = expression("Wood Density (g cm"^"-3"~")"), y = y_label) +
      # theme_minimal()
      my_theme()

    if (show_ribbon) {
      plot <- plot + geom_ribbon(data = fit_with_ci, aes(x = wd, ymin = lower, ymax = upper), alpha = 0.5, fill = color_ribbon)
    }

    return(plot)
  }

  # Generate individual plots for both H-DBH of Ang and Gym
  # p1 <- prepare_and_plot_parameter_h(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE)
  # p2 <- prepare_and_plot_parameter_h(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "b", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)
  # p3 <- prepare_and_plot_parameter_h(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "k", "#6baed6", "#8fcbee", "#2171b5", show_ribbon = TRUE)
  # p4 <- prepare_and_plot_parameter_h(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE)
  # p5 <- prepare_and_plot_parameter_h(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "b", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)
  # p6 <- prepare_and_plot_parameter_h(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "k", "#6baed6", "#8fcbee", "#2171b5", "dashed", show_ribbon = FALSE)

  p1 <- prepare_and_plot_parameter_h(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE)
  p2 <- prepare_and_plot_parameter_h(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "b", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE)
  p3 <- prepare_and_plot_parameter_h(stan_data_nlr_ang_h, fit_wd_ang_h_draws_weibull_wd, wd_df_ang, "k", "#9e9ac8", "#cbc9e2", "#6a51a3", show_ribbon = TRUE)
  p4 <- prepare_and_plot_parameter_h(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "log_a", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)
  p5 <- prepare_and_plot_parameter_h(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "b", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)
  p6 <- prepare_and_plot_parameter_h(stan_data_nlr_gym_h, fit_wd_gym_h_draws_weibull_wd, wd_df_gym, "k", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE)

  #====================
  # WD vs. CR-DBH Section
  #====================
  
  # Filter Angiosperms and Gymnosperms for CR-DBH
  ang_cr_df <- tallo_reduced_lr_df_ang_cr |> filter(division == "Angiosperm")
  gym_cr_df <- tallo_reduced_nlr_df_gym_cr |> filter(division == "Gymnosperm")

  # Group by species and calculate mean wood density for Angiosperms
  wd_df_ang_cr <- ang_cr_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Group by species and calculate mean wood density for Gymnosperms
  wd_df_gym_cr <- gym_cr_df |>
    group_by(sp) |>
    summarize(wd = mean(wd, na.rm = TRUE)) |>
    mutate(wd_s = scale(wd) |> as.numeric(), sp_id = row_number())

  # Create the sequence for wood density for Angiosperms and Gymnosperms separately
  wd_seq_ang_cr <- seq(min(wd_df_ang_cr$wd_s), max(wd_df_ang_cr$wd_s), length.out = 100)
  wd_seq_gym_cr <- seq(min(wd_df_gym_cr$wd_s), max(wd_df_gym_cr$wd_s), length.out = 100)

  wd_seq_df_ang_cr <- data.frame(
    wd_s = wd_seq_ang_cr,
    wd = approx(wd_df_ang_cr$wd_s, wd_df_ang_cr$wd, xout = wd_seq_ang_cr)$y
  )

  wd_seq_df_gym_cr <- data.frame(
    wd_s = wd_seq_gym_cr,
    wd = approx(wd_df_gym_cr$wd_s, wd_df_gym_cr$wd, xout = wd_seq_gym_cr)$y
  )

  prepare_and_plot_parameter_cr <- function(
  stan_data, 
  fit_draws, 
  wd_df, 
  wd_seq_df, 
  param_name, 
  color_point, 
  color_ribbon, 
  line_color, 
  line_type = "solid", 
  show_ribbon = TRUE, 
  is_angiosperm = TRUE
) {
  # Extract gamma parameters from fit_draws
  gamma_param <- fit_draws |> dplyr::select(starts_with("gamma"))

  # Initialize gamma_intercept and gamma_slope
  gamma_intercept <- NULL
  gamma_slope <- NULL
  
  # Handle gamma parameters for Angiosperms
  if (is_angiosperm) {
    if (param_name == "log_a") {
      gamma_intercept <- gamma_param[["gamma[1,1]"]]
      gamma_slope <- gamma_param[["gamma[1,2]"]]
    } else if (param_name == "b") {
      gamma_intercept <- gamma_param[["gamma[2,1]"]]
      gamma_slope <- gamma_param[["gamma[2,2]"]]
    } else {
      stop("Unrecognized parameter for Angiosperm: ", param_name)
    }
  } else {
    # Handle gamma parameters for Gymnosperms
    if (param_name == "log_a") {
      gamma_intercept <- gamma_param[["gamma[1,1]"]]
      gamma_slope <- gamma_param[["gamma[2,1]"]]
    } else if (param_name == "b") {
      gamma_intercept <- gamma_param[["gamma[1,2]"]]
      gamma_slope <- gamma_param[["gamma[2,2]"]]
    } else if (param_name == "k") {
      gamma_intercept <- gamma_param[["gamma[1,3]"]]
      gamma_slope <- gamma_param[["gamma[2,3]"]]
    } else {
      stop("Unrecognized parameter for Gymnosperm: ", param_name)
    }
  }

  # Reshape the beta parameters for analysis
  if (is_angiosperm) {
    # Angiosperm beta parameters (log_a, b only)
    beta_params <- fit_draws |>
      dplyr::select(starts_with("beta[")) |>
      pivot_longer(
        cols = everything(),
        names_to = c("parameter", "species"),
        names_pattern = "beta\\[(\\d+),(\\d+)\\]",
        names_transform = list(parameter = as.integer, species = as.integer),
        values_to = "value"
      ) |>
      mutate(parameter = case_when(
        parameter == 1 ~ "log_a",
        parameter == 2 ~ "b"
      ))
  } else {
    # Gymnosperm beta parameters (log_a, b, k)
    beta_params <- fit_draws |>
      dplyr::select(starts_with("beta[")) |>
      pivot_longer(
        cols = everything(),
        names_to = c("species", "parameter"),
        names_pattern = "beta\\[(\\d+),(\\d+)\\]",
        names_transform = list(species = as.integer, parameter = as.integer),
        values_to = "value"
      ) |>
      mutate(parameter = case_when(
        parameter == 1 ~ "log_a",
        parameter == 2 ~ "b",
        parameter == 3 ~ "k"
      ))
  }

  # Filter the beta parameters for the specific parameter (log_a, b, or k)
  beta_filtered <- beta_params |> filter(parameter == param_name)

  # Combine species-specific beta with wd_df
  plot_data <- beta_filtered |> left_join(wd_df, by = c("species" = "sp_id"))

  # Summarize the posterior distributions for plotting
  summary_stats <- plot_data |>
    group_by(species, parameter, wd, wd_s) |>
    summarize(
      median = median(value, na.rm = TRUE),
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  # Prepare line and CI
  lines_data <- do.call(rbind, lapply(1:length(gamma_intercept), function(i) {
    data.frame(
      wd_s = wd_seq_df$wd_s,
      wd = wd_seq_df$wd,
      value = gamma_intercept[i] + gamma_slope[i] * wd_seq_df$wd_s,
      draw = i
    )
  }))

  fit_line <- data.frame(
    wd_s = wd_seq_df$wd_s,
    wd = wd_seq_df$wd,
    value = median(gamma_intercept) + median(gamma_slope) * wd_seq_df$wd_s
  )

  ci_ribbon <- lines_data |>
    group_by(wd) |>
    summarize(
      lower = quantile(value, 0.025, na.rm = TRUE),
      upper = quantile(value, 0.975, na.rm = TRUE),
      .groups = "drop"
    )

  fit_with_ci <- left_join(fit_line, ci_ribbon, by = "wd")
  
  # Define y-labels for Angiosperms and Gymnosperms
  y_label <- if (is_angiosperm) {
    ifelse(param_name == "log_a", "Scaling factor, log a", "Exponent parameter, b")
  } else {
    case_when(
      param_name == "log_a" ~ "Asymptote parameter, log a",
      param_name == "b" ~ "Exponent parameter, b",
      param_name == "k" ~ "Scaling constant, k"
    )
  }
  
  # Create plot
  plot <- ggplot() +
    geom_errorbar(data = summary_stats, aes(x = wd, ymin = lower, ymax = upper), 
                  width = 0, linewidth = 0.2, color = color_ribbon) +
    geom_point(data = summary_stats, aes(x = wd, y = median), color = color_point, size = 0.5) +
    geom_line(data = fit_line, aes(x = wd, y = value), color = line_color, linewidth = 0.6, linetype = line_type) +
    labs(x = expression("Wood Density (g cm"^"-3"~")"), y = y_label) +
    # theme_minimal()
    my_theme()

  # Conditionally show ribbon if line_type is not "dashed"
  if (show_ribbon && line_type != "dashed") {
    plot <- plot + geom_ribbon(data = fit_with_ci, aes(x = wd, ymin = lower, ymax = upper), alpha = 0.5, fill = color_ribbon)
  }

  return(plot)
}

  # Generate individual plots for CR-DBH of Angiosperms and Gymnosperms
  # p7 <- prepare_and_plot_parameter_cr(stan_data_lr_ang_cr, fit_wd_ang_cr_draws_pl_wd, wd_df_ang_cr, wd_seq_df_ang_cr, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", is_angiosperm = TRUE)
  # p8 <- prepare_and_plot_parameter_cr(stan_data_lr_ang_cr, fit_wd_ang_cr_draws_pl_wd, wd_df_ang_cr, wd_seq_df_ang_cr, "b", "#74c476", "#bae4b3", "#238b45", is_angiosperm = TRUE)
  # p9 <- prepare_and_plot_parameter_cr(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym_cr, wd_seq_df_gym_cr, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", "dashed", show_ribbon = FALSE, is_angiosperm = FALSE)
  # p10 <- prepare_and_plot_parameter_cr(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym_cr, wd_seq_df_gym_cr, "b", "#74c476", "#bae4b3", "#238b45", show_ribbon = TRUE, is_angiosperm = FALSE)
  # p11 <- prepare_and_plot_parameter_cr(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym_cr, wd_seq_df_gym_cr, "k", "#6baed6", "#8fcbee", "#2171b5", "dashed", show_ribbon = FALSE, is_angiosperm = FALSE)
  
  p7 <- prepare_and_plot_parameter_cr(stan_data_lr_ang_cr, fit_wd_ang_cr_draws_pl_wd, wd_df_ang_cr, wd_seq_df_ang_cr, "log_a", "#9e9ac8", "#cbc9e2", "#6a51a3", is_angiosperm = TRUE)
  p8 <- prepare_and_plot_parameter_cr(stan_data_lr_ang_cr, fit_wd_ang_cr_draws_pl_wd, wd_df_ang_cr, wd_seq_df_ang_cr, "b", "#9e9ac8", "#cbc9e2", "#6a51a3", is_angiosperm = TRUE)
  p9 <- prepare_and_plot_parameter_cr(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym_cr, wd_seq_df_gym_cr, "log_a", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE, is_angiosperm = FALSE)
  p10 <- prepare_and_plot_parameter_cr(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym_cr, wd_seq_df_gym_cr, "b", "#74c476", "#bae4b3", "#238b45", show_ribbon = TRUE, is_angiosperm = FALSE)
  p11 <- prepare_and_plot_parameter_cr(stan_data_nlr_gym_cr, fit_wd_gym_cr_draws_gmm_wd, wd_df_gym_cr, wd_seq_df_gym_cr, "k", "#74c476", "#bae4b3", "#238b45", "dashed", show_ribbon = FALSE, is_angiosperm = FALSE)
  
  # Combine all panels
  p <- (p1 + labs(tag = "(a)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p2 + labs(tag = "(b)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p3 + labs(tag = "(c)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p4 + labs(tag = "(d)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p5 + labs(tag = "(e)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p6 + labs(tag = "(f)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p7 + labs(tag = "(g)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p8 + labs(tag = "(h)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      plot_spacer() |
      p9 + labs(tag = "(i)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p10 + labs(tag = "(j)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9)) |
      p11 + labs(tag = "(k)", x = expression("Wood Density (g cm"^"-3"~")")) + 
         theme(axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9))
     ) +
     plot_layout(ncol = 3, nrow = 4) +
     plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                                         "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))) & 
     my_theme() + 
     theme(
       plot.tag = element_text(size = 8),
       text = element_text(size = 5.5),
       axis.title = element_text(size = 8),
       axis.text = element_text(size = 7),
       plot.margin = margin(t = 3, r = 5, b = 3, l = 9, unit = "pt")
     )

  p <- ggdraw() +
    draw_plot(p, x = 0, y = 0, width = 1, height = 1) +
    draw_label("WD vs. H-DBH", x = 0.02, y = 0.76, angle = 90, fontface = 'bold', size = 10) + 
    draw_label("WD vs. CR-DBH", x = 0.02, y = 0.27, angle = 90, fontface = 'bold', size = 10) 

  return(p)
}



#=======================
# WOOD DENSITY's EFFECT
#=======================
# generate_wd_ef <- function(
#   tallo_reduced_nlr_df_ang_h,
#   fit_wd_ang_h_draws_weibull_wd,
#   tallo_reduced_lr_df_ang_cr,
#   fit_wd_ang_cr_draws_pl_wd,
#   tallo_reduced_nlr_df_gym_cr,
#   fit_wd_gym_cr_draws_gmm_wd
# ) {

#   dbh_seq <- seq(1, 150, length.out = 100)

#   # Prepare data for Angiosperm height
#   wd_h_ang_dat <- tallo_reduced_nlr_df_ang_h |>
#     filter(division == "Angiosperm") |>
#     group_by(sp) |>
#     summarize(wd = mean(wd)) |>
#     mutate(wd_s = scale(wd) |> as.numeric())

#   light_h_ang_wd_s <- quantile(wd_h_ang_dat$wd_s, 0.25, na.rm = TRUE)
#   dense_h_ang_wd_s <- quantile(wd_h_ang_dat$wd_s, 0.75, na.rm = TRUE)

#   gamma_params_h <- fit_wd_ang_h_draws_weibull_wd |>
#     dplyr::select(starts_with("gamma_int"), starts_with("gamma_slope"))

#   median_a <- median(exp(gamma_params_h$`gamma_int[1]`), na.rm = TRUE)
#   median_b <- median(gamma_params_h$`gamma_int[2]`, na.rm = TRUE)

#   gamma_params_h <- gamma_params_h |>
#     mutate(
#       k_light = `gamma_int[3]` + `gamma_slope[3]` * light_h_ang_wd_s,
#       k_dense = `gamma_int[3]` + `gamma_slope[3]` * dense_h_ang_wd_s
#     )

#   pred_h <- gamma_params_h |>
#     rowwise() |>
#     mutate(
#       h_light = list(tibble(
#         dbh = dbh_seq,
#         h_light = median_a * (1 - exp(-median_b * dbh^k_light))
#       )),
#       h_dense = list(tibble(
#         dbh = dbh_seq,
#         h_dense = median_a * (1 - exp(-median_b * dbh^k_dense))
#       ))
#     ) |>
#     dplyr::select(h_light, h_dense) |>
#     unnest(cols = c(h_light, h_dense), names_sep = "_") |>
#     rename(
#       dbh = h_light_dbh,
#       h_light = h_light_h_light,
#       h_dense_dbh = h_dense_dbh,
#       h_dense = h_dense_h_dense
#     )

#   h_ang_light <- pred_h |>
#     group_by(dbh) |>
#     summarize(
#       h_median_light = median(h_light),
#       h_lower_light = quantile(h_light, 0.025),
#       h_upper_light = quantile(h_light, 0.975),
#       .groups = "drop"
#     )

#   h_ang_dense <- pred_h |>
#     group_by(dbh) |>
#     summarize(
#       h_median_dense = median(h_dense),
#       h_lower_dense = quantile(h_dense, 0.025),
#       h_upper_dense = quantile(h_dense, 0.975),
#       .groups = "drop"
#     )

#   h_ang_dat <- left_join(h_ang_light, h_ang_dense, by = "dbh")

#   p1 <- ggplot(h_ang_dat, aes(x = dbh)) +
#     geom_ribbon(aes(ymin = h_lower_light, ymax = h_upper_light, fill = "Light WD species"), alpha = 0.2) +
#     geom_line(aes(y = h_median_light, color = "Light WD species")) +
#     geom_ribbon(aes(ymin = h_lower_dense, ymax = h_upper_dense, fill = "Dense WD species"), alpha = 0.2) +
#     geom_line(aes(y = h_median_dense, color = "Dense WD species")) +
#     scale_color_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
#     scale_fill_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
#     labs(x = "DBH (cm)", y = "Height (m)") +
#     my_theme()

#   # Prepare data for Angiosperm crown radius
#   wd_cr_ang_dat <- tallo_reduced_lr_df_ang_cr |>
#     filter(division == "Angiosperm") |>
#     group_by(sp) |>
#     summarize(wd = mean(wd)) |>
#     mutate(wd_s = scale(wd) |> as.numeric())

#   light_cr_ang_wd_s <- quantile(wd_cr_ang_dat$wd_s, 0.25, na.rm = TRUE)
#   dense_cr_ang_wd_s <- quantile(wd_cr_ang_dat$wd_s, 0.75, na.rm = TRUE)

#   gamma_params_cr_ang <- fit_wd_ang_cr_draws_pl_wd |>
#     dplyr::select(starts_with("gamma"))

#   gamma_params_cr_ang <- gamma_params_cr_ang |>
#     mutate(
#       a_light = exp(`gamma[1,1]` + `gamma[1,2]` * light_cr_ang_wd_s),
#       a_dense = exp(`gamma[1,1]` + `gamma[1,2]` * dense_cr_ang_wd_s),
#       b_light = `gamma[2,1]` + `gamma[2,2]` * light_cr_ang_wd_s,
#       b_dense = `gamma[2,1]` + `gamma[2,2]` * dense_cr_ang_wd_s
#     )


#   pred_ang_cr <- gamma_params_cr_ang |>
#     rowwise() |>
#     mutate(
#       cr_light = list(tibble(
#         dbh = dbh_seq,
#         cr_light = a_light * dbh^b_light
#       )),
#       cr_dense = list(tibble(
#         dbh = dbh_seq,
#         cr_dense = a_dense * dbh^b_dense
#       ))
#     ) |>
#     dplyr::select(cr_light, cr_dense) |>
#     unnest(cols = c(cr_light, cr_dense), names_sep = "_") |>
#     rename(
#       dbh = cr_light_dbh,
#       cr_light = cr_light_cr_light,
#       cr_dense_dbh = cr_dense_dbh,
#       cr_dense = cr_dense_cr_dense
#     )

#   cr_ang_light <- pred_ang_cr |>
#     group_by(dbh) |>
#     summarize(
#       cr_median_light = median(cr_light),
#       cr_lower_light = quantile(cr_light, 0.025),
#       cr_upper_light = quantile(cr_light, 0.975),
#       .groups = "drop"
#     )

#   cr_ang_dense <- pred_ang_cr |>
#     group_by(dbh) |>
#     summarize(
#       cr_median_dense = median(cr_dense),
#       cr_lower_dense = quantile(cr_dense, 0.025),
#       cr_upper_dense = quantile(cr_dense, 0.975),
#       .groups = "drop"
#     )

#   cr_ang_dat <- left_join(cr_ang_light, cr_ang_dense, by = "dbh")

#   p2 <- ggplot(cr_ang_dat, aes(x = dbh)) +
#     geom_ribbon(aes(ymin = cr_lower_light, ymax = cr_upper_light, fill = "Light WD species"), alpha = 0.2) +
#     geom_line(aes(y =  cr_median_light, color = "Light WD species")) +
#     geom_ribbon(aes(ymin = cr_lower_dense, ymax = cr_upper_dense, fill = "Dense WD species"), alpha = 0.2) +
#     geom_line(aes(y = cr_median_dense, color = "Dense WD species")) +
#     scale_color_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
#     scale_fill_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
#     labs(x = "DBH (cm)", y = "Crown Radius (m)") +
#     my_theme()


#   # Prepare data for Gymnosperm crown radius
#   wd_cr_gym_dat <- tallo_reduced_nlr_df_gym_cr |>
#     filter(division == "Gymnosperm") |>
#     group_by(sp) |>
#     summarize(wd = mean(wd)) |>
#     mutate(wd_s = scale(wd) |> as.numeric())

#   light_cr_gym_wd_s <- quantile(wd_cr_gym_dat$wd_s, 0.25, na.rm = TRUE)
#   dense_cr_gym_wd_s <- quantile(wd_cr_gym_dat$wd_s, 0.75, na.rm = TRUE)

#   gamma_params_cr_gym <- fit_wd_gym_cr_draws_gmm_wd |>
#     dplyr::select(starts_with("gamma_int"), starts_with("gamma_slope"))

#   median_cr_gym_a <- median(exp(gamma_params_cr_gym$`gamma_int[1]`), na.rm = TRUE)
#   median_cr_gym_k <- median(gamma_params_cr_gym$`gamma_int[3]`, na.rm = TRUE)

#   gamma_params_cr_gym <- gamma_params_cr_gym |>
#     mutate(
#       b_light = `gamma_int[2]` + `gamma_slope[2]` * light_cr_gym_wd_s,
#       b_dense = `gamma_int[2]` + `gamma_slope[2]` * dense_cr_gym_wd_s
#     )

#   pred_gym_cr <- gamma_params_cr_gym |>
#     rowwise() |>
#     mutate(
#       cr_light = list(tibble(
#         dbh = dbh_seq,
#         cr_light = median_cr_gym_a * dbh^b_light / (median_cr_gym_k + dbh^b_light)
#       )),
#       cr_dense = list(tibble(
#         dbh = dbh_seq,
#         cr_dense = median_cr_gym_a * dbh^b_dense / (median_cr_gym_k + dbh^b_dense)
#       ))
#     ) |>
#     dplyr::select(cr_light, cr_dense) |>
#     unnest(cols = c(cr_light, cr_dense), names_sep = "_") |>
#     rename(
#       dbh = cr_light_dbh,
#       cr_light = cr_light_cr_light,
#       cr_dense_dbh = cr_dense_dbh,
#       cr_dense = cr_dense_cr_dense
#     )

#   cr_gym_light <- pred_gym_cr |>
#     group_by(dbh) |>
#     summarize(
#       cr_median_light = median(cr_light),
#       cr_lower_light = quantile(cr_light, 0.025),
#       cr_upper_light = quantile(cr_light, 0.975),
#       .groups = "drop"
#     )

#   cr_gym_dense <- pred_gym_cr |>
#     group_by(dbh) |>
#     summarize(
#       cr_median_dense = median(cr_dense),
#       cr_lower_dense = quantile(cr_dense, 0.025),
#       cr_upper_dense = quantile(cr_dense, 0.975),
#       .groups = "drop"
#     )

#   cr_gym_dat <- left_join(cr_gym_light, cr_gym_dense, by = "dbh")

#   p3 <- ggplot(cr_gym_dat, aes(x = dbh)) +
#     geom_ribbon(aes(ymin = cr_lower_light, ymax = cr_upper_light, fill = "Light WD species"), alpha = 0.2) +
#     geom_line(aes(y =  cr_median_light, color = "Light WD species")) +
#     geom_ribbon(aes(ymin = cr_lower_dense, ymax = cr_upper_dense, fill = "Dense WD species"), alpha = 0.2) +
#     geom_line(aes(y = cr_median_dense, color = "Dense WD species")) +
#     scale_color_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
#     scale_fill_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
#     labs(x = "DBH (cm)", y = "Crown Radius (m)") +
#     my_theme()


#   p <- p1 + p2 + p3 + plot_layout(nrow = 1, guides = "collect") +
#        plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
#        theme(
#         plot.tag = element_text(size = 8),
#         axis.title = element_text(size = 9),
#         text = element_text(size = 8),
#         legend.title = element_blank(),
#         plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
#       )
  
#   p <- ggdraw() +
#   draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +  
#   draw_label("Angiosperm", x = 0.172, y = 0.94, fontface = 'bold', size = 10) +  
#   draw_label("Angiosperm", x = 0.429, y = 0.94, fontface = 'bold', size = 10) +  
#   draw_label("Gymnosperm", x = 0.689, y = 0.94, fontface = 'bold', size = 10)

#   return(p)
# }

generate_wd_ef <- function(
  tallo_reduced_nlr_df_ang_h,
  fit_wd_ang_h_draws_weibull_wd,
  tallo_reduced_lr_df_ang_cr,
  fit_wd_ang_cr_draws_pl_wd,
  tallo_reduced_nlr_df_gym_cr,
  fit_wd_gym_cr_draws_gmm_wd
) {
  dbh_seq_h <- seq(1, 100, length.out = 100)
  dbh_seq_cr_ang <- seq(1, 100, length.out = 100)
  dbh_seq_cr_gym <- seq(1, 100, length.out = 100)

  process_data <- function(data, division, fit_draws, light_quantile, dense_quantile, model_type) {
    wd_dat <- data |>
      filter(division == !!division) |>
      group_by(sp) |>
      summarize(wd = mean(wd), .groups = "drop") |>
      mutate(wd_s = scale(wd) |> as.numeric())

    light_wd_s <- quantile(wd_dat$wd_s, light_quantile, na.rm = TRUE)
    dense_wd_s <- quantile(wd_dat$wd_s, dense_quantile, na.rm = TRUE)

    gamma_params <- fit_draws |>
      dplyr::select(starts_with("gamma"))

    if (model_type == "weibull") {
      median_a <- median(exp(gamma_params$`gamma_int[1]`), na.rm = TRUE)
      median_b <- median(gamma_params$`gamma_int[2]`, na.rm = TRUE)

      gamma_params <- gamma_params |>
        mutate(
          k_light = `gamma_int[3]` + `gamma_slope[3]` * light_wd_s,
          k_dense = `gamma_int[3]` + `gamma_slope[3]` * dense_wd_s
        )

      pred <- gamma_params |>
        rowwise() |>
        mutate(
          light = list(tibble(
            dbh = dbh_seq_h,
            value = median_a * (1 - exp(-median_b * dbh^k_light))
          )),
          dense = list(tibble(
            dbh = dbh_seq_h,
            value = median_a * (1 - exp(-median_b * dbh^k_dense))
          ))
        )
    } else if (model_type == "power_law") {
      gamma_params <- gamma_params |>
        mutate(
          a_light = exp(`gamma[1,1]` + `gamma[1,2]` * light_wd_s),
          a_dense = exp(`gamma[1,1]` + `gamma[1,2]` * dense_wd_s),
          b_light = `gamma[2,1]` + `gamma[2,2]` * light_wd_s,
          b_dense = `gamma[2,1]` + `gamma[2,2]` * dense_wd_s
        )

      pred <- gamma_params |>
        rowwise() |>
        mutate(
          light = list(tibble(
            dbh = dbh_seq_cr_ang,
            value = a_light * dbh^b_light
          )),
          dense = list(tibble(
            dbh = dbh_seq_cr_ang,
            value = a_dense * dbh^b_dense
          ))
        )
    } else if (model_type == "gmm") {
      median_a <- median(exp(gamma_params$`gamma_int[1]`), na.rm = TRUE)
      median_k <- median(gamma_params$`gamma_int[3]`, na.rm = TRUE)

      gamma_params <- gamma_params |>
        mutate(
          b_light = `gamma_int[2]` + `gamma_slope[2]` * light_wd_s,
          b_dense = `gamma_int[2]` + `gamma_slope[2]` * dense_wd_s
        )

      pred <- gamma_params |>
        rowwise() |>
        mutate(
          light = list(tibble(
            dbh = dbh_seq_cr_gym,
            value = median_a * dbh^b_light / (median_k + dbh^b_light)
          )),
          dense = list(tibble(
            dbh = dbh_seq_cr_gym,
            value = median_a * dbh^b_dense / (median_k + dbh^b_dense)
          ))
        )
    } else {
      stop("Unsupported model type")
    }

    pred |>
      dplyr::select(light, dense) |>
      unnest(cols = c(light, dense), names_sep = "_") |>
      rename(
        dbh = light_dbh,
        light = light_value,
        dense = dense_value
      ) |>
      group_by(dbh) |>
      summarize(
        median_light = median(light),
        lower_light = quantile(light, 0.025),
        upper_light = quantile(light, 0.975),
        median_dense = median(dense),
        lower_dense = quantile(dense, 0.025),
        upper_dense = quantile(dense, 0.975),
        .groups = "drop"
      )
  }

  plot_data <- function(data, x_label, y_label) {
    ggplot(data, aes(x = dbh)) +
      geom_ribbon(aes(ymin = lower_light, ymax = upper_light, fill = "Light WD species"), alpha = 0.2) +
      geom_line(aes(y = median_light, color = "Light WD species")) +
      geom_ribbon(aes(ymin = lower_dense, ymax = upper_dense, fill = "Dense WD species"), alpha = 0.2) +
      geom_line(aes(y = median_dense, color = "Dense WD species")) +
      scale_color_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
      scale_fill_manual(name = "Wood Type", values = c("Light WD species" = "#f48d10", "Dense WD species" = "#7f590e")) +
      labs(x = x_label, y = y_label) +
      my_theme()
  }

  h_ang_dat <- process_data(
    tallo_reduced_nlr_df_ang_h,
    "Angiosperm",
    fit_wd_ang_h_draws_weibull_wd,
    0.1,
    0.9,
    "weibull"
  )
  p1 <- plot_data(h_ang_dat, "DBH (cm)", "Height (m)")

  cr_ang_dat <- process_data(
    tallo_reduced_lr_df_ang_cr,
    "Angiosperm",
    fit_wd_ang_cr_draws_pl_wd,
    0.1,
    0.9,
    "power_law"
  )
  p2 <- plot_data(cr_ang_dat, "DBH (cm)", "Crown Radius (m)")

  cr_gym_dat <- process_data(
    tallo_reduced_nlr_df_gym_cr,
    "Gymnosperm",
    fit_wd_gym_cr_draws_gmm_wd,
    0.1,
    0.9,
    "gmm"
  )
  p3 <- plot_data(cr_gym_dat, "DBH (cm)", "Crown Radius (m)")

  p2_y_limits <- ggplot_build(p2)$layout$panel_scales_y[[1]]$range$range
  p3_y_limits <- ggplot_build(p3)$layout$panel_scales_y[[1]]$range$range

  # Set the same limits for both p2 and p3
  shared_y_limits <- range(p2_y_limits, p3_y_limits)

  # Recreate p2 and p3 with shared y-axis limits
  p2 <- p2 + ylim(shared_y_limits)
  p3 <- p3 + ylim(shared_y_limits)

  p <- p1 + p2 + p3 + plot_layout(nrow = 1, guides = "collect") +
       plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)"))) &
       theme(
        plot.tag = element_text(size = 8),
        axis.title = element_text(size = 9),
        text = element_text(size = 8),
        legend.title = element_blank(),
        plot.margin = margin(t = 2, r = 2, b = 2, l = 2, unit = "pt")
      )
  
  p <- ggdraw() +
  draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +  
  draw_label("Angiosperm", x = 0.172, y = 0.94, fontface = 'bold', size = 10) +  
  draw_label("Angiosperm", x = 0.429, y = 0.94, fontface = 'bold', size = 10) +  
  draw_label("Gymnosperm", x = 0.689, y = 0.94, fontface = 'bold', size = 10)

  return(p)
}

#======================
# AGB ESTIMATION
#======================

#' AGB Estimation Function
#'
#' @description This function generates Above-Ground Biomass (AGB) predictions using different models 
#' (Weibull and Power-Law) for both Angiosperms and Gymnosperms. It computes AGB estimates based on
#' species-level posterior parameters and wood density data, and evaluates model performance through
#' bias and RMSE metrics. The function also generates diagnostic plots comparing observed vs. predicted
#' AGB for multiple models and allows for optional export of results in YAML format.
#'
#' @param tallo_wd_df0 Data frame containing the original wood density and tree measurements.
#' @param sp_posterior_agb_df Data frame of posterior distributions for species-specific parameters 
#' (Weibull and Power-Law models).
#' @param export_yaml Logical flag to indicate whether the AGB metrics (RMSE and Bias) should be
#' exported to a YAML file. Default is FALSE.
#' @param yaml_file File path where the YAML metrics will be exported if `export_yaml` is TRUE.
#' 
#' @details The function first splits the data into Angiosperms and Gymnosperms based on the 
#' `Division` column. It then joins the species-level posterior distributions (Weibull and Power-Law) 
#' with wood density data, computes AGB using different model formulations, and evaluates model performance 
#' by calculating RMSE and Bias for each model. The function generates diagnostic plots for each model, 
#' comparing predicted vs. observed AGB on a log-log scale.
#' 

generate_agb_estimation_com <- function(tallo_wd_df0, sp_posterior_agb_df, export_yaml = FALSE, yaml_file = NULL) {
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }
  # Data preparation
  tallo_wd_df <- tallo_wd_df0 |>
    filter(!is.na(wd)) |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_wd = log(wd)) |>
    filter(!is.na(log_dbh), !is.na(log_cr), !is.na(log_h), !is.na(log_wd)) |>
    dplyr::select(tree_id, division, sp, dbh, h, cr, wd, log_dbh, log_h, log_cr, log_wd) |>
    group_by(sp) |>
    filter(n() >= 20) |>
    ungroup()
  
  # Split the tallo_wd_df data by Division
  tallo_wd_ang <- tallo_wd_df |>
    filter(division == "Angiosperm")

  tallo_wd_gym <- tallo_wd_df |>
    filter(division == "Gymnosperm")

  # Splitting the posterior data into Angiosperms and Gymnosperms
  sp_posterior_ang <- sp_posterior_agb_df |> filter(Division == "Angiosperm")
  sp_posterior_gym <- sp_posterior_agb_df |> filter(Division == "Gymnosperm")

  # Function to filter common species across all relevant datasets
  filter_common_species <- function(tallo_wd_df, sp_posterior_agb_df) {
    h_wb_df <- sp_posterior_agb_df |>
    filter(Dependent_variable == "Tree Height", !is.na(a), !is.na(b),!is.na(k)) |>
      dplyr::select(sp, a, b, k) |>
      distinct()

    h_pl_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", is.na(k)) |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH", !is.na(a), !is.na(b), !is.na(c)) |>
      dplyr::select(sp, a, b, c) |>
      distinct()

    dbh1_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH1") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh2_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH2") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh3_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH3") |>
      dplyr::select(sp, a, b) |>
      distinct()

    common_species <- Reduce(intersect, list(
      tallo_wd_df$sp,
      h_wb_df$sp,
      h_pl_df$sp,
      dbh_df$sp,
      dbh1_df$sp,
      dbh2_df$sp,
      dbh3_df$sp
    ))

    list(
      filtered_tallo_wd_df = tallo_wd_df |> filter(sp %in% common_species),
      filtered_h_wb_df = h_wb_df |> filter(sp %in% common_species),
      filtered_h_pl_df = h_pl_df |> filter(sp %in% common_species),
      filtered_dbh_df = dbh_df |> filter(sp %in% common_species),
      filtered_dbh1_df = dbh1_df |> filter(sp %in% common_species),
      filtered_dbh2_df = dbh2_df |> filter(sp %in% common_species),
      filtered_dbh3_df = dbh3_df |> filter(sp %in% common_species),
      num_species = length(common_species),
      num_trees = tallo_wd_df |> filter(sp %in% common_species) |> nrow()

    )
  }

  # Filter common species for Angiosperms and Gymnosperms
  filtered_ang <- filter_common_species(tallo_wd_ang, sp_posterior_ang)
  filtered_gym <- filter_common_species(tallo_wd_gym, sp_posterior_gym)

 
  # Function to join posterior data and estimate AGB
  calculate_agb <- function(filtered_data) {
    tallo_wd_df <- filtered_data$filtered_tallo_wd_df
    agb_df <- tallo_wd_df |>
      left_join(
        filtered_data$filtered_h_wb_df |> rename(a_h_wb = a, b_h_wb = b, k_h_wb = k),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_h_pl_df |> rename(a_h_pl = a, b_h_pl = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh_df |> rename(a_dbh = a, b_dbh = b, c_dbh = c),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh1_df |> rename(a_dbh1 = a, b_dbh1 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh2_df |> rename(a_dbh2 = a, b_dbh2 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh3_df |> rename(a_dbh3 = a, b_dbh3 = b),
        by = "sp"
      )

    # Calculate the log_AGB based on various models
    agb_df <- agb_df |>
      mutate(
        a_h_wb = as.numeric(gsub(" \\(.*\\)", "", a_h_wb)),
        b_h_wb = as.numeric(gsub(" \\(.*\\)", "", b_h_wb)),
        k_h_wb = as.numeric(gsub(" \\(.*\\)", "", k_h_wb)),
        a_h_pl = as.numeric(gsub(" \\(.*\\)", "", a_h_pl)),
        b_h_pl = as.numeric(gsub(" \\(.*\\)", "", b_h_pl)),
        a_dbh = as.numeric(gsub(" \\(.*\\)", "", a_dbh)),
        b_dbh = as.numeric(gsub(" \\(.*\\)", "", b_dbh)),
        c_dbh = as.numeric(gsub(" \\(.*\\)", "", c_dbh)),
        a_dbh1 = as.numeric(gsub(" \\(.*\\)", "", a_dbh1)),
        b_dbh1 = as.numeric(gsub(" \\(.*\\)", "", b_dbh1)),
        a_dbh2 = as.numeric(gsub(" \\(.*\\)", "", a_dbh2)),
        b_dbh2 = as.numeric(gsub(" \\(.*\\)", "", b_dbh2)),
        a_dbh3 = as.numeric(gsub(" \\(.*\\)", "", a_dbh3)),
        b_dbh3 = as.numeric(gsub(" \\(.*\\)", "", b_dbh3))
      ) |>
      mutate(
        log_AGB_bl = log(0.0559) + log_wd + 2 * log_dbh + log_h,
        
        height_wb = a_h_wb * (1 - exp(-b_h_wb * dbh^k_h_wb)),    # Weibull-based height estimation
        log_AGB_wb_h = log(0.0559) + log_wd + 2 * log_dbh + log(height_wb),  # Biomass from Weibull height
        
        height_pl = a_h_pl * dbh^b_h_pl,  # Power-law height estimation
        log_AGB_pl_h = log(0.0559) + log_wd + 2 * log_dbh + log(height_pl),  # Biomass from Power-law height
        
        dbh_pl = a_dbh * (cr^b_dbh) * (h^c_dbh),  # Power-law DBH estimation using CR and H
        log_AGB_pl_dbh = log(0.0559) + log_wd + 2 * log(dbh_pl) + log_h,  # Biomass from Power-law DBH
        
        dbh_pl1 = a_dbh1 * (cr * h)^b_dbh1,  # Power-law DBH estimation using CR*H
        log_AGB_pl_dbh1 = log(0.0559) + log_wd + 2 * log(dbh_pl1) + log_h,  # Biomass from Power-law DBH1
        
        dbh_pl2 = a_dbh2 * (cr^b_dbh2),  # Power-law DBH estimation using CR only
        log_AGB_pl_dbh2 = log(0.0559) + log_wd + 2 * log(dbh_pl2) + log_h,  # Biomass from Power-law DBH2
        
        dbh_pl3 = a_dbh3 * (h^b_dbh3),  # Power-law DBH estimation using H only
        log_AGB_pl_dbh3 = log(0.0559) + log_wd + 2 * log(dbh_pl3) + log_h  # Biomass from Power-law DBH3
      )

    return(agb_df)
  }

  # Calculate AGB for Angiosperms and Gymnosperms
  agb_ang <- calculate_agb(filtered_ang)
  agb_gym <- calculate_agb(filtered_gym)

  # Combine metrics from Angiosperms and Gymnosperms
  combine_metrics <- function(agb_df) {
    agb_df |>
      summarize(
        rmse_wb_h = sqrt(mean((log_AGB_bl - log_AGB_wb_h)^2, na.rm = TRUE)),
        bias_wb_h = mean(log_AGB_wb_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_h = sqrt(mean((log_AGB_bl - log_AGB_pl_h)^2, na.rm = TRUE)),
        bias_pl_h = mean(log_AGB_pl_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh)^2, na.rm = TRUE)),
        bias_pl_dbh = mean(log_AGB_pl_dbh - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh1 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh1)^2, na.rm = TRUE)),
        bias_pl_dbh1 = mean(log_AGB_pl_dbh1 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh2 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh2)^2, na.rm = TRUE)),
        bias_pl_dbh2 = mean(log_AGB_pl_dbh2 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh3 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh3)^2, na.rm = TRUE)),
        bias_pl_dbh3 = mean(log_AGB_pl_dbh3 - log_AGB_bl, na.rm = TRUE)
      )
  }

  metrics_ang <- combine_metrics(agb_ang)
  metrics_gym <- combine_metrics(agb_gym)
  # Add species and tree counts to metrics
  metrics_ang$num_species <- filtered_ang$num_species
  metrics_ang$num_trees <- filtered_ang$num_trees
  metrics_gym$num_species <- filtered_gym$num_species
  metrics_gym$num_trees <- filtered_gym$num_trees

  metrics_ang$num_species <- format_number(metrics_ang$num_species)
  metrics_ang$num_trees <- format_number(metrics_ang$num_trees)
  metrics_gym$num_species <- format_number(metrics_gym$num_species)
  metrics_gym$num_trees <- format_number(metrics_gym$num_trees)

  # Define the plotting function for Angiosperms and Gymnosperms
  plot_agb <- function(agb_df, metrics_df, model_colors, division_label) {
      create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
        ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
          geom_point(alpha = 0.5, color = "grey", size = 0.3) +
          geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
          scale_x_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
          scale_y_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
            labs(tag = tag_label, 
                x = x_axis_label, 
                y = switch(model_name,
                            "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
                            "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
                            "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " + ", H)]~"(kg)"),
                            "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
                            "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
                            "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

          annotate("text", x = 0.1, y = 1e5, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 1.5) +
          annotate("text", x = 0.1, y = 3e4, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 1.5) +
          my_theme()
      }

      model_colors <- c(
        "h_wb" = "#D55E00",
        "h_pl" = "#0072B2",
        "dbh_pl" = "#E69F00",
        "dbh1_pl" = "#F0E442",
        "dbh2_pl" = "#009E73",
        "dbh3_pl" = "#56B4E9"
      )

      # Create all 6 subplots for one division (Angiosperm or Gymnosperm)
      p1 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_wb_h, 
                        metrics_df$rmse_wb_h, metrics_df$bias_wb_h, "h_wb", "(a)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p2 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_h, 
                        metrics_df$rmse_pl_h, metrics_df$bias_pl_h, "h_pl", "(b)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p3 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh, 
                        metrics_df$rmse_pl_dbh, metrics_df$bias_pl_dbh, "dbh_pl", "(c)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p4 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh1, 
                        metrics_df$rmse_pl_dbh1, metrics_df$bias_pl_dbh1, "dbh1_pl", "(d)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p5 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh2, 
                        metrics_df$rmse_pl_dbh2, metrics_df$bias_pl_dbh2, "dbh2_pl", "(e)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p6 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh3, 
                        metrics_df$rmse_pl_dbh3, metrics_df$bias_pl_dbh3, "dbh3_pl", "(f)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      return(list(p1, p2, p3, p4, p5, p6))
    }

    # Generate individual plots for Angiosperms and Gymnosperms
    ang_plots <- plot_agb(agb_ang, metrics_ang, model_colors, "Angiosperm")
    gym_plots <- plot_agb(agb_gym, metrics_gym, model_colors, "Gymnosperm")

    # Combine the plots into a 4x3 grid: 4 columns and 3 rows
    p <- (ang_plots[[1]] | ang_plots[[2]] |
          gym_plots[[1]] | gym_plots[[2]] |
          ang_plots[[3]] | ang_plots[[4]] |
          gym_plots[[3]] | gym_plots[[4]] |
          ang_plots[[5]] | ang_plots[[6]] |
          gym_plots[[5]] | gym_plots[[6]]) +
      plot_layout(ncol = 4, nrow = 3) +
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                                          "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))) &
      my_theme() +
      theme(
        text = element_text(size = 5.5),
        axis.title = element_text(size = 6.5),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.margin = margin(t = 3, r = 3, b = 3, l = 2, unit = "pt")
      )
    
    p <- ggdraw() +
      draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label("Angiosperm", x = 0.28, y = 0.965, fontface = 'bold', size = 10) +  
      draw_label("Gymnosperm", x = 0.78, y = 0.965, fontface = 'bold', size = 10)

  if (export_yaml) {
    yaml::write_yaml(list(metrics_ang = as.list(metrics_ang), metrics_gym = as.list(metrics_gym)), yaml_file)
  }

    return(list(p = p, metrics_ang = metrics_ang, metrics_gym = metrics_gym))
}


#=================================
# AGB ESTIMATION ADDITIONAL Eq1
#================================
generate_agb_estimation_com1 <- function(tallo_wd_df0, sp_posterior_agb_df, export_yaml = FALSE, yaml_file = NULL) {
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }
  # Data preparation
  tallo_wd_df <- tallo_wd_df0 |>
    filter(!is.na(wd)) |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_wd = log(wd)) |>
    filter(!is.na(log_dbh), !is.na(log_cr), !is.na(log_h), !is.na(log_wd)) |>
    dplyr::select(tree_id, division, sp, dbh, h, cr, wd, log_dbh, log_h, log_cr, log_wd) |>
    group_by(sp) |>
    filter(n() >= 20) |>
    ungroup()
  
  # Split the tallo_wd_df data by Division
  tallo_wd_ang <- tallo_wd_df |>
    filter(division == "Angiosperm")

  tallo_wd_gym <- tallo_wd_df |>
    filter(division == "Gymnosperm")

  # Splitting the posterior data into Angiosperms and Gymnosperms
  sp_posterior_ang <- sp_posterior_agb_df |> filter(Division == "Angiosperm")
  sp_posterior_gym <- sp_posterior_agb_df |> filter(Division == "Gymnosperm")

  # Function to filter common species across all relevant datasets
  filter_common_species <- function(tallo_wd_df, sp_posterior_agb_df) {
    h_wb_df <- sp_posterior_agb_df |>
    filter(Dependent_variable == "Tree Height", !is.na(a), !is.na(b),!is.na(k)) |>
      dplyr::select(sp, a, b, k) |>
      distinct()

    h_pl_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", is.na(k)) |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH", !is.na(a), !is.na(b), !is.na(c)) |>
      dplyr::select(sp, a, b, c) |>
      distinct()

    dbh1_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH1") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh2_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH2") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh3_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH3") |>
      dplyr::select(sp, a, b) |>
      distinct()

    common_species <- Reduce(intersect, list(
      tallo_wd_df$sp,
      h_wb_df$sp,
      h_pl_df$sp,
      dbh_df$sp,
      dbh1_df$sp,
      dbh2_df$sp,
      dbh3_df$sp
    ))

    list(
      filtered_tallo_wd_df = tallo_wd_df |> filter(sp %in% common_species),
      filtered_h_wb_df = h_wb_df |> filter(sp %in% common_species),
      filtered_h_pl_df = h_pl_df |> filter(sp %in% common_species),
      filtered_dbh_df = dbh_df |> filter(sp %in% common_species),
      filtered_dbh1_df = dbh1_df |> filter(sp %in% common_species),
      filtered_dbh2_df = dbh2_df |> filter(sp %in% common_species),
      filtered_dbh3_df = dbh3_df |> filter(sp %in% common_species),
      num_species = length(common_species),
      num_trees = tallo_wd_df |> filter(sp %in% common_species) |> nrow()

    )
  }

  # Filter common species for Angiosperms and Gymnosperms
  filtered_ang <- filter_common_species(tallo_wd_ang, sp_posterior_ang)
  filtered_gym <- filter_common_species(tallo_wd_gym, sp_posterior_gym)

 
calculate_agb <- function(filtered_data) {
    tallo_wd_df <- filtered_data$filtered_tallo_wd_df
    agb_df <- tallo_wd_df |>
      left_join(
        filtered_data$filtered_h_wb_df |> rename(a_h_wb = a, b_h_wb = b, k_h_wb = k),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_h_pl_df |> rename(a_h_pl = a, b_h_pl = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh_df |> rename(a_dbh = a, b_dbh = b, c_dbh = c),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh1_df |> rename(a_dbh1 = a, b_dbh1 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh2_df |> rename(a_dbh2 = a, b_dbh2 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh3_df |> rename(a_dbh3 = a, b_dbh3 = b),
        by = "sp"
      )

    # Calculate the log_AGB based on various models
    agb_df <- agb_df |>
      mutate(
        a_h_wb = as.numeric(gsub(" \\(.*\\)", "", a_h_wb)),
        b_h_wb = as.numeric(gsub(" \\(.*\\)", "", b_h_wb)),
        k_h_wb = as.numeric(gsub(" \\(.*\\)", "", k_h_wb)),
        a_h_pl = as.numeric(gsub(" \\(.*\\)", "", a_h_pl)),
        b_h_pl = as.numeric(gsub(" \\(.*\\)", "", b_h_pl)),
        a_dbh = as.numeric(gsub(" \\(.*\\)", "", a_dbh)),
        b_dbh = as.numeric(gsub(" \\(.*\\)", "", b_dbh)),
        c_dbh = as.numeric(gsub(" \\(.*\\)", "", c_dbh)),
        a_dbh1 = as.numeric(gsub(" \\(.*\\)", "", a_dbh1)),
        b_dbh1 = as.numeric(gsub(" \\(.*\\)", "", b_dbh1)),
        a_dbh2 = as.numeric(gsub(" \\(.*\\)", "", a_dbh2)),
        b_dbh2 = as.numeric(gsub(" \\(.*\\)", "", b_dbh2)),
        a_dbh3 = as.numeric(gsub(" \\(.*\\)", "", a_dbh3)),
        b_dbh3 = as.numeric(gsub(" \\(.*\\)", "", b_dbh3))
      ) |>
      mutate(
        # Updated baseline AGB calculation with new equation
        log_AGB_bl = log(0.0673) + 0.976 * (log_wd + 2 * log_dbh + log_h),
        
        height_wb = a_h_wb * (1 - exp(-b_h_wb * dbh^k_h_wb)),    # Weibull-based height estimation
        log_AGB_wb_h = log(0.0673) + 0.976 * (log_wd + 2 * log_dbh + log(height_wb)),  # Biomass from Weibull height
        
        height_pl = a_h_pl * dbh^b_h_pl,  # Power-law height estimation
        log_AGB_pl_h = log(0.0673) + 0.976 * (log_wd + 2 * log_dbh + log(height_pl)),  # Biomass from Power-law height
        
        dbh_pl = a_dbh * (cr^b_dbh) * (h^c_dbh),  # Power-law DBH estimation using CR and H
        log_AGB_pl_dbh = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl) + log_h),  # Biomass from Power-law DBH
        
        dbh_pl1 = a_dbh1 * (cr * h)^b_dbh1,  # Power-law DBH estimation using CR*H
        log_AGB_pl_dbh1 = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl1) + log_h),  # Biomass from Power-law DBH1
        
        dbh_pl2 = a_dbh2 * (cr^b_dbh2),  # Power-law DBH estimation using CR only
        log_AGB_pl_dbh2 = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl2) + log_h),  # Biomass from Power-law DBH2
        
        dbh_pl3 = a_dbh3 * (h^b_dbh3),  # Power-law DBH estimation using H only
        log_AGB_pl_dbh3 = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl3) + log_h)  # Biomass from Power-law DBH3
      )

    return(agb_df)
  }

  # Calculate AGB for Angiosperms and Gymnosperms
  agb_ang <- calculate_agb(filtered_ang)
  agb_gym <- calculate_agb(filtered_gym)

  # Combine metrics from Angiosperms and Gymnosperms
  combine_metrics <- function(agb_df) {
    agb_df |>
      summarize(
        rmse_wb_h = sqrt(mean((log_AGB_bl - log_AGB_wb_h)^2, na.rm = TRUE)),
        bias_wb_h = mean(log_AGB_wb_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_h = sqrt(mean((log_AGB_bl - log_AGB_pl_h)^2, na.rm = TRUE)),
        bias_pl_h = mean(log_AGB_pl_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh)^2, na.rm = TRUE)),
        bias_pl_dbh = mean(log_AGB_pl_dbh - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh1 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh1)^2, na.rm = TRUE)),
        bias_pl_dbh1 = mean(log_AGB_pl_dbh1 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh2 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh2)^2, na.rm = TRUE)),
        bias_pl_dbh2 = mean(log_AGB_pl_dbh2 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh3 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh3)^2, na.rm = TRUE)),
        bias_pl_dbh3 = mean(log_AGB_pl_dbh3 - log_AGB_bl, na.rm = TRUE)
      )
  }

  metrics_ang1 <- combine_metrics(agb_ang)
  metrics_gym1 <- combine_metrics(agb_gym)
  # Add species and tree counts to metrics
  metrics_ang1$num_species <- filtered_ang$num_species
  metrics_ang1$num_trees <- filtered_ang$num_trees
  metrics_gym1$num_species <- filtered_gym$num_species
  metrics_gym1$num_trees <- filtered_gym$num_trees

  metrics_ang1$num_species <- format_number(metrics_ang1$num_species)
  metrics_ang1$num_trees <- format_number(metrics_ang1$num_trees)
  metrics_gym1$num_species <- format_number(metrics_gym1$num_species)
  metrics_gym1$num_trees <- format_number(metrics_gym1$num_trees)

  # Define the plotting function for Angiosperms and Gymnosperms
  plot_agb <- function(agb_df, metrics_df, model_colors, division_label) {
      create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
        ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
          geom_point(alpha = 0.5, color = "grey", size = 0.3) +
          geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
          scale_x_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
          scale_y_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
            labs(tag = tag_label, 
                x = x_axis_label, 
                y = switch(model_name,
                            "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
                            "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
                            "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " + ", H)]~"(kg)"),
                            "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
                            "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
                            "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

          annotate("text", x = 0.1, y = 1e5, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 1.5) +
          annotate("text", x = 0.1, y = 3e4, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 1.5) +
          my_theme()
      }

      model_colors <- c(
        "h_wb" = "#D55E00",
        "h_pl" = "#0072B2",
        "dbh_pl" = "#E69F00",
        "dbh1_pl" = "#F0E442",
        "dbh2_pl" = "#009E73",
        "dbh3_pl" = "#56B4E9"
      )

      # Create all 6 subplots for one division (Angiosperm or Gymnosperm)
      p1 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_wb_h, 
                        metrics_df$rmse_wb_h, metrics_df$bias_wb_h, "h_wb", "(a)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p2 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_h, 
                        metrics_df$rmse_pl_h, metrics_df$bias_pl_h, "h_pl", "(b)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p3 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh, 
                        metrics_df$rmse_pl_dbh, metrics_df$bias_pl_dbh, "dbh_pl", "(c)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p4 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh1, 
                        metrics_df$rmse_pl_dbh1, metrics_df$bias_pl_dbh1, "dbh1_pl", "(d)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p5 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh2, 
                        metrics_df$rmse_pl_dbh2, metrics_df$bias_pl_dbh2, "dbh2_pl", "(e)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p6 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh3, 
                        metrics_df$rmse_pl_dbh3, metrics_df$bias_pl_dbh3, "dbh3_pl", "(f)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      return(list(p1, p2, p3, p4, p5, p6))
    }

    # Generate individual plots for Angiosperms and Gymnosperms
    ang_plots <- plot_agb(agb_ang, metrics_ang1, model_colors, "Angiosperm")
    gym_plots <- plot_agb(agb_gym, metrics_gym1, model_colors, "Gymnosperm")

    # Combine the plots into a 4x3 grid: 4 columns and 3 rows
    p <- (ang_plots[[1]] | ang_plots[[2]] |
          gym_plots[[1]] | gym_plots[[2]] |
          ang_plots[[3]] | ang_plots[[4]] |
          gym_plots[[3]] | gym_plots[[4]] |
          ang_plots[[5]] | ang_plots[[6]] |
          gym_plots[[5]] | gym_plots[[6]]) +
      plot_layout(ncol = 4, nrow = 3) +
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                                          "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))) &
      my_theme() +
      theme(
        text = element_text(size = 5.5),
        axis.title = element_text(size = 6.5),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.margin = margin(t = 3, r = 3, b = 3, l = 2, unit = "pt")
      )
    
    p <- ggdraw() +
      draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label("Angiosperm", x = 0.28, y = 0.965, fontface = 'bold', size = 10) +  
      draw_label("Gymnosperm", x = 0.78, y = 0.965, fontface = 'bold', size = 10)

  if (export_yaml) {
    yaml::write_yaml(list(metrics_ang1 = as.list(metrics_ang1), metrics_gym1 = as.list(metrics_gym1)), yaml_file)
  }

    return(list(p = p, metrics_ang1 = metrics_ang1, metrics_gym1 = metrics_gym1))
}


#===========================================
# AGB ESTIMATION ADDITIONAL Eq1 (ANGIOSPERM)
#===========================================

generate_agb_estimation <- function(tallo_wd_df0, sp_posterior_agb_df, export_yaml = FALSE, yaml_file = NULL) {
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }
  # Data preparation
  tallo_wd_df <- tallo_wd_df0 |>
    filter(!is.na(wd)) |>
    filter(dbh >= 5) |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_wd = log(wd)) |>
    filter(!is.na(log_dbh), !is.na(log_cr), !is.na(log_h), !is.na(log_wd)) |>
    dplyr::select(tree_id, division, sp, dbh, h, cr, wd, log_dbh, log_h, log_cr, log_wd) |>
    group_by(sp) |>
    filter(n() >= 20) |>
    ungroup()
  
  tallo_wd_ang <- tallo_wd_df |>
    filter(division == "Angiosperm")

  # Splitting the posterior data into Angiosperms and Gymnosperms
  sp_posterior_ang <- sp_posterior_agb_df |> filter(Division == "Angiosperm")

  # Function to filter common species across all relevant datasets
  filter_common_species <- function(tallo_wd_df, sp_posterior_agb_df) {
    h_wb_df <- sp_posterior_agb_df |>
    filter(Dependent_variable == "Tree Height", !is.na(a), !is.na(b),!is.na(k)) |>
      dplyr::select(sp, a, b, k) |>
      distinct()

    h_pl_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", is.na(k)) |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH", !is.na(a), !is.na(b), !is.na(c)) |>
      dplyr::select(sp, a, b, c) |>
      distinct()

    dbh1_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH1") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh2_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH2") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh3_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH3") |>
      dplyr::select(sp, a, b) |>
      distinct()

    common_species <- Reduce(intersect, list(
      tallo_wd_df$sp,
      h_wb_df$sp,
      h_pl_df$sp,
      dbh_df$sp,
      dbh1_df$sp,
      dbh2_df$sp,
      dbh3_df$sp
    ))

    list(
      filtered_tallo_wd_df = tallo_wd_df |> filter(sp %in% common_species),
      filtered_h_wb_df = h_wb_df |> filter(sp %in% common_species),
      filtered_h_pl_df = h_pl_df |> filter(sp %in% common_species),
      filtered_dbh_df = dbh_df |> filter(sp %in% common_species),
      filtered_dbh1_df = dbh1_df |> filter(sp %in% common_species),
      filtered_dbh2_df = dbh2_df |> filter(sp %in% common_species),
      filtered_dbh3_df = dbh3_df |> filter(sp %in% common_species),
      num_species = length(common_species),
      num_trees = tallo_wd_df |> filter(sp %in% common_species) |> nrow()

    )
  }

  # Filter common species for Angiosperms
  filtered_ang <- filter_common_species(tallo_wd_ang, sp_posterior_ang)

 
calculate_agb <- function(filtered_data) {
    tallo_wd_df <- filtered_data$filtered_tallo_wd_df
    agb_df <- tallo_wd_df |>
      left_join(
        filtered_data$filtered_h_wb_df |> rename(a_h_wb = a, b_h_wb = b, k_h_wb = k),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_h_pl_df |> rename(a_h_pl = a, b_h_pl = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh_df |> rename(a_dbh = a, b_dbh = b, c_dbh = c),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh1_df |> rename(a_dbh1 = a, b_dbh1 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh2_df |> rename(a_dbh2 = a, b_dbh2 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh3_df |> rename(a_dbh3 = a, b_dbh3 = b),
        by = "sp"
      )

    # Calculate the log_AGB based on various models
    agb_df <- agb_df |>
      mutate(
        a_h_wb = as.numeric(gsub(" \\(.*\\)", "", a_h_wb)),
        b_h_wb = as.numeric(gsub(" \\(.*\\)", "", b_h_wb)),
        k_h_wb = as.numeric(gsub(" \\(.*\\)", "", k_h_wb)),
        a_h_pl = as.numeric(gsub(" \\(.*\\)", "", a_h_pl)),
        b_h_pl = as.numeric(gsub(" \\(.*\\)", "", b_h_pl)),
        a_dbh = as.numeric(gsub(" \\(.*\\)", "", a_dbh)),
        b_dbh = as.numeric(gsub(" \\(.*\\)", "", b_dbh)),
        c_dbh = as.numeric(gsub(" \\(.*\\)", "", c_dbh)),
        a_dbh1 = as.numeric(gsub(" \\(.*\\)", "", a_dbh1)),
        b_dbh1 = as.numeric(gsub(" \\(.*\\)", "", b_dbh1)),
        a_dbh2 = as.numeric(gsub(" \\(.*\\)", "", a_dbh2)),
        b_dbh2 = as.numeric(gsub(" \\(.*\\)", "", b_dbh2)),
        a_dbh3 = as.numeric(gsub(" \\(.*\\)", "", a_dbh3)),
        b_dbh3 = as.numeric(gsub(" \\(.*\\)", "", b_dbh3))
      ) |>
      mutate(
        # Updated baseline AGB calculation with new equation
        log_AGB_bl = log(0.0673) + 0.976 * (log_wd + 2 * log_dbh + log_h),
        
        height_wb = a_h_wb * (1 - exp(-b_h_wb * dbh^k_h_wb)),    # Weibull-based height estimation
        log_AGB_wb_h = log(0.0673) + 0.976 * (log_wd + 2 * log_dbh + log(height_wb)),  # Biomass from Weibull height
        
        height_pl = a_h_pl * dbh^b_h_pl,  # Power-law height estimation
        log_AGB_pl_h = log(0.0673) + 0.976 * (log_wd + 2 * log_dbh + log(height_pl)),  # Biomass from Power-law height
        
        dbh_pl = a_dbh * (cr^b_dbh) * (h^c_dbh),  # Power-law DBH estimation using CR and H
        log_AGB_pl_dbh = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl) + log_h),  # Biomass from Power-law DBH
        
        dbh_pl1 = a_dbh1 * (cr * h)^b_dbh1,  # Power-law DBH estimation using CR*H
        log_AGB_pl_dbh1 = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl1) + log_h),  # Biomass from Power-law DBH1
        
        dbh_pl2 = a_dbh2 * (cr^b_dbh2),  # Power-law DBH estimation using CR only
        log_AGB_pl_dbh2 = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl2) + log_h),  # Biomass from Power-law DBH2
        
        dbh_pl3 = a_dbh3 * (h^b_dbh3),  # Power-law DBH estimation using H only
        log_AGB_pl_dbh3 = log(0.0673) + 0.976 * (log_wd + 2 * log(dbh_pl3) + log_h)  # Biomass from Power-law DBH3
      )

    return(agb_df)
  }

  # Calculate AGB for Angiosperms
  agb_ang <- calculate_agb(filtered_ang)

  # Combine metrics from Angiosperms
  combine_metrics <- function(agb_df) {
    agb_df |>
      summarize(
        rmse_wb_h = sqrt(mean((log_AGB_bl - log_AGB_wb_h)^2, na.rm = TRUE)),
        bias_wb_h = mean(log_AGB_wb_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_h = sqrt(mean((log_AGB_bl - log_AGB_pl_h)^2, na.rm = TRUE)),
        bias_pl_h = mean(log_AGB_pl_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh)^2, na.rm = TRUE)),
        bias_pl_dbh = mean(log_AGB_pl_dbh - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh1 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh1)^2, na.rm = TRUE)),
        bias_pl_dbh1 = mean(log_AGB_pl_dbh1 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh2 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh2)^2, na.rm = TRUE)),
        bias_pl_dbh2 = mean(log_AGB_pl_dbh2 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh3 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh3)^2, na.rm = TRUE)),
        bias_pl_dbh3 = mean(log_AGB_pl_dbh3 - log_AGB_bl, na.rm = TRUE)
      )
  }

  metrics_agb_ang <- combine_metrics(agb_ang)

  # Add species and tree counts to metrics
  metrics_agb_ang$num_species <- filtered_ang$num_species
  metrics_agb_ang$num_trees <- filtered_ang$num_trees

  metrics_agb_ang$num_species <- format_number(metrics_agb_ang$num_species)
  metrics_agb_ang$num_trees <- format_number(metrics_agb_ang$num_trees)

  # Define the plotting function for Angiosperms and Gymnosperms
  plot_agb <- function(agb_df, metrics_df, model_colors, division_label) {
      # create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
      #   ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
      #     geom_point(alpha = 0.5, color = "grey", size = 0.3) +
      #     geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
      #     geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
      #     scale_x_log10(
      #       breaks = c(0.1, 10, 1000, 1e5), 
      #       # limits = c(NA, 1e5),
      #       labels = scales::trans_format("log10", scales::math_format(10^.x))
      #     ) +
      #     scale_y_log10(
      #       breaks = c(0.1, 10, 1000, 1e5), 
      #       # limits = c(NA, 1e5),
      #       labels = scales::trans_format("log10", scales::math_format(10^.x))
      #     ) +
      #       labs(tag = tag_label, 
      #           x = x_axis_label, 
      #           y = switch(model_name,
      #                       "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
      #                       "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
      #                       "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " , ", H)]~"(kg)"),
      #                       "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
      #                       "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
      #                       "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

      #     annotate("text", x = 0.1, y = 1e5, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 2) +
      #     annotate("text", x = 0.1, y = 3e4, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 2) +
      #     my_theme()
      # }
    create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
      # Calculate the range of predicted values to adjust the y position dynamically
      max_predicted <- max(exp(predicted), na.rm = TRUE)  # Using exp to scale back from log

      # Adjust positions based on the max value of predicted AGB
      y_position_rmse <- max_predicted * 0.9  # Place RMSE slightly below the top
      y_position_bias <- max_predicted * 0.4  # Place Bias slightly lower

      ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
        geom_point(alpha = 0.5, color = "grey", size = 0.3) +
        geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
        scale_x_log10(
          breaks = c(0.1, 10, 1000, 1e5), 
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        scale_y_log10(
          breaks = c(0.1, 10, 1000, 1e5), 
          labels = scales::trans_format("log10", scales::math_format(10^.x))
        ) +
        labs(tag = tag_label, 
            x = x_axis_label, 
            y = switch(model_name,
                        "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
                        "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
                        "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " , ", H)]~"(kg)"),
                        "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
                        "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
                        "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

        # Annotate RMSE and Bias with dynamic positions
        # annotate("text", x = 0.1, y = y_position_rmse, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 2) +
        # annotate("text", x = 0.1, y = y_position_bias, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 2) +
        annotate("text", x = 0.1, y = y_position_rmse, 
                label = paste("RMSE =", sprintf("%.3g", rmse), "Mg"), 
                hjust = 0, size = 2) +
        annotate("text", x = 0.1, y = y_position_bias, 
                label = paste("Bias =", sprintf("%.3g", bias * 100), "%"), 
                hjust = 0, size = 2) +

        my_theme()
    }

      model_colors <- c(
        "h_wb" = "#D55E00",
        "h_pl" = "#0072B2",
        "dbh_pl" = "#E69F00",
        "dbh1_pl" = "#F0E442",
        "dbh2_pl" = "#009E73",
        "dbh3_pl" = "#56B4E9"
      )

      # Create all 6 subplots for one division (Angiosperm or Gymnosperm)
      p1 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_wb_h, 
                        metrics_df$rmse_wb_h, metrics_df$bias_wb_h, "h_wb", "(a)",
                        x_axis_label = expression(AGB[ref]~"(kg)"))

      p2 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_h, 
                        metrics_df$rmse_pl_h, metrics_df$bias_pl_h, "h_pl", "(b)",
                        x_axis_label = expression(AGB[ref]~"(kg)"))

      p3 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh, 
                        metrics_df$rmse_pl_dbh, metrics_df$bias_pl_dbh, "dbh_pl", "(c)",
                        x_axis_label = expression(AGB[ref]~"(kg)"))

      p4 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh1, 
                        metrics_df$rmse_pl_dbh1, metrics_df$bias_pl_dbh1, "dbh1_pl", "(d)",
                        x_axis_label = expression(AGB[ref]~"(kg)"))

      p5 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh2, 
                        metrics_df$rmse_pl_dbh2, metrics_df$bias_pl_dbh2, "dbh2_pl", "(e)",
                        x_axis_label = expression(AGB[ref]~"(kg)"))

      p6 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh3, 
                        metrics_df$rmse_pl_dbh3, metrics_df$bias_pl_dbh3, "dbh3_pl", "(f)",
                        x_axis_label = expression(AGB[ref]~"(kg)"))

      return(list(p1, p2, p3, p4, p5, p6))
    }

    # Generate individual plots for Angiosperms
    ang_plots <- plot_agb(agb_ang, metrics_agb_ang, model_colors, "Angiosperm")

    # Combine the plots into a 4x3 grid: 4 columns and 3 rows
    p <- (ang_plots[[1]] | ang_plots[[2]] |
          ang_plots[[3]] | ang_plots[[4]] |
          ang_plots[[6]] | ang_plots[[5]]) +
      plot_layout(ncol = 3, nrow = 2) +
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)"))) &
      my_theme() +
      theme(
        text = element_text(size = 6.5),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 7),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 6),
        plot.margin = margin(t = 3, r = 3, b = 3, l = 2.5, unit = "pt")
      )
    
    p <- ggdraw() +
      draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label("Angiosperm", x = 0.545, y = 0.965, fontface = 'bold', size = 10)

  if (export_yaml) {
    yaml::write_yaml(list(metrics_agb_ang = as.list(metrics_agb_ang)), yaml_file)
  }

    return(list(p = p, metrics_agb_ang = metrics_agb_ang))
}

#=================================
# AGB ESTIMATION ADDITIONAL Eq2
#================================
generate_agb_estimation_com2 <- function(tallo_wd_df0, sp_posterior_agb_df, export_yaml = FALSE, yaml_file = NULL) {
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }
  # Data preparation
  tallo_wd_df <- tallo_wd_df0 |>
    filter(!is.na(wd)) |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_wd = log(wd)) |>
    filter(!is.na(log_dbh), !is.na(log_cr), !is.na(log_h), !is.na(log_wd)) |>
    dplyr::select(tree_id, division, sp, dbh, h, cr, wd, log_dbh, log_h, log_cr, log_wd) |>
    group_by(sp) |>
    filter(n() >= 20) |>
    ungroup()

  # Split the tallo_wd_df data by Division
  tallo_wd_ang <- tallo_wd_df |>
    filter(division == "Angiosperm")

  tallo_wd_gym <- tallo_wd_df |>
    filter(division == "Gymnosperm")

  # Function to filter common species across all relevant datasets
  filter_common_species <- function(tallo_wd_df, sp_posterior_agb_df) {
    h_wb_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", !is.na(a), !is.na(b),!is.na(k)) |>
      dplyr::select(sp, a, b, k) |>
      distinct()

    h_pl_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", is.na(k)) |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH", !is.na(a), !is.na(b), !is.na(c)) |>
      dplyr::select(sp, a, b, c) |>
      distinct()

    dbh1_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH1") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh2_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH2") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh3_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH3") |>
      dplyr::select(sp, a, b) |>
      distinct()

    common_species <- Reduce(intersect, list(
      tallo_wd_df$sp,
      h_wb_df$sp,
      h_pl_df$sp,
      dbh_df$sp,
      dbh1_df$sp,
      dbh2_df$sp,
      dbh3_df$sp
    ))

    list(
      filtered_tallo_wd_df = tallo_wd_df |> filter(sp %in% common_species),
      filtered_h_wb_df = h_wb_df |> filter(sp %in% common_species),
      filtered_h_pl_df = h_pl_df |> filter(sp %in% common_species),
      filtered_dbh_df = dbh_df |> filter(sp %in% common_species),
      filtered_dbh1_df = dbh1_df |> filter(sp %in% common_species),
      filtered_dbh2_df = dbh2_df |> filter(sp %in% common_species),
      filtered_dbh3_df = dbh3_df |> filter(sp %in% common_species),
      num_species = length(common_species),
      num_trees = tallo_wd_df |> filter(sp %in% common_species) |> nrow()
    )
  }

  # Filter common species for Angiosperms and Gymnosperms
  filtered_ang <- filter_common_species(tallo_wd_ang, sp_posterior_agb_df)
  filtered_gym <- filter_common_species(tallo_wd_gym, sp_posterior_agb_df)
  
  # Function to join posterior data and estimate AGB
  calculate_agb <- function(filtered_data) {
    tallo_wd_df <- filtered_data$filtered_tallo_wd_df
    agb_df <- tallo_wd_df |>
      left_join(
        filtered_data$filtered_h_wb_df |> rename(a_h_wb = a, b_h_wb = b, k_h_wb = k),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_h_pl_df |> rename(a_h_pl = a, b_h_pl = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh_df |> rename(a_dbh = a, b_dbh = b, c_dbh = c),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh1_df |> rename(a_dbh1 = a, b_dbh1 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh2_df |> rename(a_dbh2 = a, b_dbh2 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh3_df |> rename(a_dbh3 = a, b_dbh3 = b),
        by = "sp"
      )

    # Calculate the log_AGB based on various models
    agb_df <- agb_df |>
      mutate(
        a_h_wb = as.numeric(gsub(" \\(.*\\)", "", a_h_wb)),
        b_h_wb = as.numeric(gsub(" \\(.*\\)", "", b_h_wb)),
        k_h_wb = as.numeric(gsub(" \\(.*\\)", "", k_h_wb)),
        a_h_pl = as.numeric(gsub(" \\(.*\\)", "", a_h_pl)),
        b_h_pl = as.numeric(gsub(" \\(.*\\)", "", b_h_pl)),
        a_dbh = as.numeric(gsub(" \\(.*\\)", "", a_dbh)),
        b_dbh = as.numeric(gsub(" \\(.*\\)", "", b_dbh)),
        c_dbh = as.numeric(gsub(" \\(.*\\)", "", c_dbh)),
        a_dbh1 = as.numeric(gsub(" \\(.*\\)", "", a_dbh1)),
        b_dbh1 = as.numeric(gsub(" \\(.*\\)", "", b_dbh1)),
        a_dbh2 = as.numeric(gsub(" \\(.*\\)", "", a_dbh2)),
        b_dbh2 = as.numeric(gsub(" \\(.*\\)", "", b_dbh2)),
        a_dbh3 = as.numeric(gsub(" \\(.*\\)", "", a_dbh3)),
        b_dbh3 = as.numeric(gsub(" \\(.*\\)", "", b_dbh3))
      ) |>
      mutate(
        log_AGB_bl = log(0.046) + 0.883 * (2 * log_dbh + log_h),

        height_wb = a_h_wb * (1 - exp(-b_h_wb * dbh^k_h_wb)),
        log_AGB_wb_h = log(0.046) + 0.883 * (2 * log_dbh + log(height_wb)),

        height_pl = a_h_pl * dbh^b_h_pl,
        log_AGB_pl_h = log(0.046) + 0.883 * (2 * log_dbh + log(height_pl)),

        dbh_pl = a_dbh * (cr^b_dbh) * (h^c_dbh),
        log_AGB_pl_dbh = log(0.046) + 0.883 * (2 * log(dbh_pl) + log_h),

        dbh_pl1 = a_dbh1 * (cr * h)^b_dbh1,
        log_AGB_pl_dbh1 = log(0.046) + 0.883 * (2 * log(dbh_pl1) + log_h),

        dbh_pl2 = a_dbh2 * (cr^b_dbh2),
        log_AGB_pl_dbh2 = log(0.046) + 0.883 * (2 * log(dbh_pl2) + log_h),

        dbh_pl3 = a_dbh3 * (h^b_dbh3),
        log_AGB_pl_dbh3 = log(0.046) + 0.883 * (2 * log(dbh_pl3) + log_h)
      )

    return(agb_df)
  }

  agb_ang <- calculate_agb(filtered_ang)
  agb_gym <- calculate_agb(filtered_gym)

  combine_metrics <- function(agb_df) {
    agb_df |>
      summarize(
        rmse_wb_h = sqrt(mean((log_AGB_bl - log_AGB_wb_h)^2, na.rm = TRUE)),
        bias_wb_h = mean(log_AGB_wb_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_h = sqrt(mean((log_AGB_bl - log_AGB_pl_h)^2, na.rm = TRUE)),
        bias_pl_h = mean(log_AGB_pl_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh)^2, na.rm = TRUE)),
        bias_pl_dbh = mean(log_AGB_pl_dbh - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh1 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh1)^2, na.rm = TRUE)),
        bias_pl_dbh1 = mean(log_AGB_pl_dbh1 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh2 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh2)^2, na.rm = TRUE)),
        bias_pl_dbh2 = mean(log_AGB_pl_dbh2 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh3 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh3)^2, na.rm = TRUE)),
        bias_pl_dbh3 = mean(log_AGB_pl_dbh3 - log_AGB_bl, na.rm = TRUE)
      )
  }

  metrics_ang2 <- combine_metrics(agb_ang)
  metrics_gym2 <- combine_metrics(agb_gym)
 # Add species and tree counts to metrics
  metrics_ang2$num_species <- filtered_ang$num_species
  metrics_ang2$num_trees <- filtered_ang$num_trees
  metrics_gym2$num_species <- filtered_gym$num_species
  metrics_gym2$num_trees <- filtered_gym$num_trees

  metrics_ang2$num_species <- format_number(metrics_ang2$num_species)
  metrics_ang2$num_trees <- format_number(metrics_ang2$num_trees)
  metrics_gym2$num_species <- format_number(metrics_gym2$num_species)
  metrics_gym2$num_trees <- format_number(metrics_gym2$num_trees)

  # Define the plotting function for Angiosperms and Gymnosperms
  plot_agb <- function(agb_df, metrics_df, model_colors, division_label) {
      create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
        ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
          geom_point(alpha = 0.5, color = "grey", size = 0.3) +
          geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
          scale_x_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
          scale_y_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
            labs(tag = tag_label, 
                x = x_axis_label, 
                y = switch(model_name,
                            "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
                            "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
                            "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, "  ", H)]~"(kg)"),
                            "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
                            "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
                            "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

          annotate("text", x = 0.1, y = 1e5, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 1.5) +
          annotate("text", x = 0.1, y = 3e4, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 1.5) +
          my_theme()
      }

      model_colors <- c(
        "h_wb" = "#D55E00",
        "h_pl" = "#0072B2",
        "dbh_pl" = "#E69F00",
        "dbh1_pl" = "#F0E442",
        "dbh2_pl" = "#009E73",
        "dbh3_pl" = "#56B4E9"
      )

      # Create all 6 subplots for one division (Angiosperm or Gymnosperm)
      p1 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_wb_h, 
                        metrics_df$rmse_wb_h, metrics_df$bias_wb_h, "h_wb", "(a)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p2 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_h, 
                        metrics_df$rmse_pl_h, metrics_df$bias_pl_h, "h_pl", "(b)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p3 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh, 
                        metrics_df$rmse_pl_dbh, metrics_df$bias_pl_dbh, "dbh_pl", "(c)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p4 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh1, 
                        metrics_df$rmse_pl_dbh1, metrics_df$bias_pl_dbh1, "dbh1_pl", "(d)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p5 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh2, 
                        metrics_df$rmse_pl_dbh2, metrics_df$bias_pl_dbh2, "dbh2_pl", "(e)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p6 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh3, 
                        metrics_df$rmse_pl_dbh3, metrics_df$bias_pl_dbh3, "dbh3_pl", "(f)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      return(list(p1, p2, p3, p4, p5, p6))
    }

    # Generate individual plots for Angiosperms and Gymnosperms
    ang_plots <- plot_agb(agb_ang, metrics_ang2, model_colors, "Angiosperm")
    gym_plots <- plot_agb(agb_gym, metrics_gym2, model_colors, "Gymnosperm")

    # Combine the plots into a 4x3 grid: 4 columns and 3 rows
    p <- (ang_plots[[1]] | ang_plots[[2]] |
          gym_plots[[1]] | gym_plots[[2]] |
          ang_plots[[3]] | ang_plots[[4]] |
          gym_plots[[3]] | gym_plots[[4]] |
          ang_plots[[5]] | ang_plots[[6]] |
          gym_plots[[5]] | gym_plots[[6]]) +
      plot_layout(ncol = 4, nrow = 3) +
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                                          "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))) &
      my_theme() +
      theme(
        text = element_text(size = 5.5),
        axis.title = element_text(size = 6.5),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.margin = margin(t = 3, r = 3, b = 3, l = 2, unit = "pt")
      )
    
    p <- ggdraw() +
      draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label("Angiosperm", x = 0.28, y = 0.965, fontface = 'bold', size = 10) +  
      draw_label("Gymnosperm", x = 0.78, y = 0.965, fontface = 'bold', size = 10)

  if (export_yaml) {
    yaml::write_yaml(list(metrics_ang2 = as.list(metrics_ang2), metrics_gym2 = as.list(metrics_gym2)), yaml_file)
  }

  return(list(p = p, metrics_ang2 = metrics_ang2, metrics_gym2 = metrics_gym2))
}


#=================================
# AGB ESTIMATION ADDITIONAL Eq3
#================================
generate_agb_estimation_com3 <- function(tallo_wd_df0, sp_posterior_agb_df, export_yaml = FALSE, yaml_file = NULL) {
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }

  # Data preparation
  tallo_wd_df <- tallo_wd_df0 |>
    filter(!is.na(wd)) |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_wd = log(wd)) |>
    filter(!is.na(log_dbh), !is.na(log_cr), !is.na(log_h), !is.na(log_wd)) |>
    dplyr::select(tree_id, division, sp, dbh, h, cr, wd, log_dbh, log_h, log_cr, log_wd) |>
    group_by(sp) |>
    filter(n() >= 20) |>
    ungroup()

  # Split the tallo_wd_df data by Division
  tallo_wd_ang <- tallo_wd_df |>
    filter(division == "Angiosperm")

  tallo_wd_gym <- tallo_wd_df |>
    filter(division == "Gymnosperm")

  # Function to filter common species across all relevant datasets
  filter_common_species <- function(tallo_wd_df, sp_posterior_agb_df) {
    h_wb_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", !is.na(a), !is.na(b),!is.na(k)) |>
      dplyr::select(sp, a, b, k) |>
      distinct()

    h_pl_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", is.na(k)) |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH", !is.na(a), !is.na(b), !is.na(c)) |>
      dplyr::select(sp, a, b, c) |>
      distinct()

    dbh1_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH1") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh2_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH2") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh3_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH3") |>
      dplyr::select(sp, a, b) |>
      distinct()

    common_species <- Reduce(intersect, list(
      tallo_wd_df$sp,
      h_wb_df$sp,
      h_pl_df$sp,
      dbh_df$sp,
      dbh1_df$sp,
      dbh2_df$sp,
      dbh3_df$sp
    ))

    list(
      filtered_tallo_wd_df = tallo_wd_df |> filter(sp %in% common_species),
      filtered_h_wb_df = h_wb_df |> filter(sp %in% common_species),
      filtered_h_pl_df = h_pl_df |> filter(sp %in% common_species),
      filtered_dbh_df = dbh_df |> filter(sp %in% common_species),
      filtered_dbh1_df = dbh1_df |> filter(sp %in% common_species),
      filtered_dbh2_df = dbh2_df |> filter(sp %in% common_species),
      filtered_dbh3_df = dbh3_df |> filter(sp %in% common_species),
      num_species = length(common_species),
      num_trees = tallo_wd_df |> filter(sp %in% common_species) |> nrow()
    )
  }

  # Filter common species for Angiosperms and Gymnosperms
  filtered_ang <- filter_common_species(tallo_wd_ang, sp_posterior_agb_df)
  filtered_gym <- filter_common_species(tallo_wd_gym, sp_posterior_agb_df)

  # Function to join posterior data and estimate AGB
  calculate_agb <- function(filtered_data) {
    tallo_wd_df <- filtered_data$filtered_tallo_wd_df
    agb_df <- tallo_wd_df |>
      left_join(
        filtered_data$filtered_h_wb_df |> rename(a_h_wb = a, b_h_wb = b, k_h_wb = k),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_h_pl_df |> rename(a_h_pl = a, b_h_pl = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh_df |> rename(a_dbh = a, b_dbh = b, c_dbh = c),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh1_df |> rename(a_dbh1 = a, b_dbh1 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh2_df |> rename(a_dbh2 = a, b_dbh2 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh3_df |> rename(a_dbh3 = a, b_dbh3 = b),
        by = "sp"
      )

    # Calculate the log_AGB based on various models
    agb_df <- agb_df |>
      mutate(
        a_h_wb = as.numeric(gsub(" \\(.*\\)", "", a_h_wb)),
        b_h_wb = as.numeric(gsub(" \\(.*\\)", "", b_h_wb)),
        k_h_wb = as.numeric(gsub(" \\(.*\\)", "", k_h_wb)),
        a_h_pl = as.numeric(gsub(" \\(.*\\)", "", a_h_pl)),
        b_h_pl = as.numeric(gsub(" \\(.*\\)", "", b_h_pl)),
        a_dbh = as.numeric(gsub(" \\(.*\\)", "", a_dbh)),
        b_dbh = as.numeric(gsub(" \\(.*\\)", "", b_dbh)),
        c_dbh = as.numeric(gsub(" \\(.*\\)", "", c_dbh)),
        a_dbh1 = as.numeric(gsub(" \\(.*\\)", "", a_dbh1)),
        b_dbh1 = as.numeric(gsub(" \\(.*\\)", "", b_dbh1)),
        a_dbh2 = as.numeric(gsub(" \\(.*\\)", "", a_dbh2)),
        b_dbh2 = as.numeric(gsub(" \\(.*\\)", "", b_dbh2)),
        a_dbh3 = as.numeric(gsub(" \\(.*\\)", "", a_dbh3)),
        b_dbh3 = as.numeric(gsub(" \\(.*\\)", "", b_dbh3))
      ) |>
      mutate(
        log_AGB_bl = log(0.039) + 2.43 * log_dbh + 0.114 * log_h,

        height_wb = a_h_wb * (1 - exp(-b_h_wb * dbh^k_h_wb)),
        log_AGB_wb_h = log(0.039) + 2.43 * log_dbh + 0.114 * log(height_wb),

        height_pl = a_h_pl * dbh^b_h_pl,
        log_AGB_pl_h = log(0.039) + 2.43 * log_dbh + 0.114 * log(height_pl),

        dbh_pl = a_dbh * (cr^b_dbh) * (h^c_dbh),
        log_AGB_pl_dbh = log(0.039) + 2.43 * log(dbh_pl) + 0.114 * log_h,

        dbh_pl1 = a_dbh1 * (cr * h)^b_dbh1,
        log_AGB_pl_dbh1 = log(0.039) + 2.43 * log(dbh_pl1) + 0.114 * log_h,

        dbh_pl2 = a_dbh2 * (cr^b_dbh2),
        log_AGB_pl_dbh2 = log(0.039) + 2.43 * log(dbh_pl2) + 0.114 * log_h,

        dbh_pl3 = a_dbh3 * (h^b_dbh3),
        log_AGB_pl_dbh3 = log(0.039) + 2.43 * log(dbh_pl3) + 0.114 * log_h
      )

    return(agb_df)
  }

  agb_ang <- calculate_agb(filtered_ang)
  agb_gym <- calculate_agb(filtered_gym)

  combine_metrics <- function(agb_df) {
    agb_df |>
      summarize(
        rmse_wb_h = sqrt(mean((log_AGB_bl - log_AGB_wb_h)^2, na.rm = TRUE)),
        bias_wb_h = mean(log_AGB_wb_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_h = sqrt(mean((log_AGB_bl - log_AGB_pl_h)^2, na.rm = TRUE)),
        bias_pl_h = mean(log_AGB_pl_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh)^2, na.rm = TRUE)),
        bias_pl_dbh = mean(log_AGB_pl_dbh - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh1 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh1)^2, na.rm = TRUE)),
        bias_pl_dbh1 = mean(log_AGB_pl_dbh1 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh2 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh2)^2, na.rm = TRUE)),
        bias_pl_dbh2 = mean(log_AGB_pl_dbh2 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh3 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh3)^2, na.rm = TRUE)),
        bias_pl_dbh3 = mean(log_AGB_pl_dbh3 - log_AGB_bl, na.rm = TRUE)
      )
  }

  metrics_ang3 <- combine_metrics(agb_ang)
  metrics_gym3 <- combine_metrics(agb_gym)

  # Add species and tree counts to metrics
  metrics_ang3$num_species <- filtered_ang$num_species
  metrics_ang3$num_trees <- filtered_ang$num_trees
  metrics_gym3$num_species <- filtered_gym$num_species
  metrics_gym3$num_trees <- filtered_gym$num_trees

  metrics_ang3$num_species <- format_number(metrics_ang3$num_species)
  metrics_ang3$num_trees <- format_number(metrics_ang3$num_trees)
  metrics_gym3$num_species <- format_number(metrics_gym3$num_species)
  metrics_gym3$num_trees <- format_number(metrics_gym3$num_trees)

  # Define the plotting function for Angiosperms and Gymnosperms
  plot_agb <- function(agb_df, metrics_df, model_colors, division_label) {
      create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
        ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
          geom_point(alpha = 0.5, color = "grey", size = 0.3) +
          geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
          scale_x_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
          scale_y_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
            labs(tag = tag_label, 
                x = x_axis_label, 
                y = switch(model_name,
                            "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
                            "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
                            "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " + ", H)]~"(kg)"),
                            "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
                            "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
                            "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

          annotate("text", x = 0.1, y = 1e5, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 1.5) +
          annotate("text", x = 0.1, y = 3e4, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 1.5) +
          my_theme()
      }

      model_colors <- c(
        "h_wb" = "#D55E00",
        "h_pl" = "#0072B2",
        "dbh_pl" = "#E69F00",
        "dbh1_pl" = "#F0E442",
        "dbh2_pl" = "#009E73",
        "dbh3_pl" = "#56B4E9"
      )

      # Create all 6 subplots for one division (Angiosperm or Gymnosperm)
      p1 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_wb_h, 
                        metrics_df$rmse_wb_h, metrics_df$bias_wb_h, "h_wb", "(a)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p2 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_h, 
                        metrics_df$rmse_pl_h, metrics_df$bias_pl_h, "h_pl", "(b)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p3 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh, 
                        metrics_df$rmse_pl_dbh, metrics_df$bias_pl_dbh, "dbh_pl", "(c)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p4 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh1, 
                        metrics_df$rmse_pl_dbh1, metrics_df$bias_pl_dbh1, "dbh1_pl", "(d)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p5 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh2, 
                        metrics_df$rmse_pl_dbh2, metrics_df$bias_pl_dbh2, "dbh2_pl", "(e)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p6 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh3, 
                        metrics_df$rmse_pl_dbh3, metrics_df$bias_pl_dbh3, "dbh3_pl", "(f)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      return(list(p1, p2, p3, p4, p5, p6))
    }

    # Generate individual plots for Angiosperms and Gymnosperms
    ang_plots <- plot_agb(agb_ang, metrics_ang3, model_colors, "Angiosperm")
    gym_plots <- plot_agb(agb_gym, metrics_gym3, model_colors, "Gymnosperm")

    # Combine the plots into a 4x3 grid: 4 columns and 3 rows
    p <- (ang_plots[[1]] | ang_plots[[2]] |
          gym_plots[[1]] | gym_plots[[2]] |
          ang_plots[[3]] | ang_plots[[4]] |
          gym_plots[[3]] | gym_plots[[4]] |
          ang_plots[[5]] | ang_plots[[6]] |
          gym_plots[[5]] | gym_plots[[6]]) +
      plot_layout(ncol = 4, nrow = 3) +
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                                          "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))) &
      my_theme() +
      theme(
        text = element_text(size = 5.5),
        axis.title = element_text(size = 6.5),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.margin = margin(t = 3, r = 3, b = 3, l = 2, unit = "pt")
      )
    
    p <- ggdraw() +
      draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label("Angiosperm", x = 0.28, y = 0.965, fontface = 'bold', size = 10) +  
      draw_label("Gymnosperm", x = 0.78, y = 0.965, fontface = 'bold', size = 10)

  if (export_yaml) {
    yaml::write_yaml(list(metrics_ang3 = as.list(metrics_ang3), metrics_gym3 = as.list(metrics_gym3)), yaml_file)
  }

  return(list(p = p, metrics_ang3 = metrics_ang3, metrics_gym3 = metrics_gym3))
}

#=================================
# AGB ESTIMATION ADDITIONAL Eq4
#================================
generate_agb_estimation_com4 <- function(tallo_wd_df0, sp_posterior_agb_df, export_yaml = FALSE, yaml_file = NULL) {
  format_number <- function(x) {
    format(x, big.mark = ",", scientific = FALSE)
  }

  # Data preparation
  tallo_wd_df <- tallo_wd_df0 |>
    filter(!is.na(wd)) |>
    mutate(log_dbh = log(dbh),
           log_h = log(h),
           log_cr = log(cr),
           log_wd = log(wd)) |>
    filter(!is.na(log_dbh), !is.na(log_cr), !is.na(log_h), !is.na(log_wd)) |>
    dplyr::select(tree_id, division, sp, dbh, h, cr, wd, log_dbh, log_h, log_cr, log_wd) |>
    group_by(sp) |>
    filter(n() >= 20) |>
    ungroup()

  # Split the tallo_wd_df data by Division
  tallo_wd_ang <- tallo_wd_df |>
    filter(division == "Angiosperm")

  tallo_wd_gym <- tallo_wd_df |>
    filter(division == "Gymnosperm")

  # Function to filter common species across all relevant datasets
  filter_common_species <- function(tallo_wd_df, sp_posterior_agb_df) {
    h_wb_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", !is.na(a), !is.na(b),!is.na(k)) |>
      dplyr::select(sp, a, b, k) |>
      distinct()

    h_pl_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "Tree Height", is.na(k)) |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH", !is.na(a), !is.na(b), !is.na(c)) |>
      dplyr::select(sp, a, b, c) |>
      distinct()

    dbh1_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH1") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh2_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH2") |>
      dplyr::select(sp, a, b) |>
      distinct()

    dbh3_df <- sp_posterior_agb_df |>
      filter(Dependent_variable == "DBH3") |>
      dplyr::select(sp, a, b) |>
      distinct()

    common_species <- Reduce(intersect, list(
      tallo_wd_df$sp,
      h_wb_df$sp,
      h_pl_df$sp,
      dbh_df$sp,
      dbh1_df$sp,
      dbh2_df$sp,
      dbh3_df$sp
    ))

    list(
      filtered_tallo_wd_df = tallo_wd_df |> filter(sp %in% common_species),
      filtered_h_wb_df = h_wb_df |> filter(sp %in% common_species),
      filtered_h_pl_df = h_pl_df |> filter(sp %in% common_species),
      filtered_dbh_df = dbh_df |> filter(sp %in% common_species),
      filtered_dbh1_df = dbh1_df |> filter(sp %in% common_species),
      filtered_dbh2_df = dbh2_df |> filter(sp %in% common_species),
      filtered_dbh3_df = dbh3_df |> filter(sp %in% common_species),
      num_species = length(common_species),
      num_trees = tallo_wd_df |> filter(sp %in% common_species) |> nrow()
    )
  }

  # Filter common species for Angiosperms and Gymnosperms
  filtered_ang <- filter_common_species(tallo_wd_ang, sp_posterior_agb_df)
  filtered_gym <- filter_common_species(tallo_wd_gym, sp_posterior_agb_df)

  # Function to join posterior data and estimate AGB
  calculate_agb <- function(filtered_data) {
    tallo_wd_df <- filtered_data$filtered_tallo_wd_df
    agb_df <- tallo_wd_df |>
      left_join(
        filtered_data$filtered_h_wb_df |> rename(a_h_wb = a, b_h_wb = b, k_h_wb = k),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_h_pl_df |> rename(a_h_pl = a, b_h_pl = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh_df |> rename(a_dbh = a, b_dbh = b, c_dbh = c),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh1_df |> rename(a_dbh1 = a, b_dbh1 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh2_df |> rename(a_dbh2 = a, b_dbh2 = b),
        by = "sp"
      ) |>
      left_join(
        filtered_data$filtered_dbh3_df |> rename(a_dbh3 = a, b_dbh3 = b),
        by = "sp"
      )

    # Calculate the log_AGB based on various models
    agb_df <- agb_df |>
      mutate(
        a_h_wb = as.numeric(gsub(" \\(.*\\)", "", a_h_wb)),
        b_h_wb = as.numeric(gsub(" \\(.*\\)", "", b_h_wb)),
        k_h_wb = as.numeric(gsub(" \\(.*\\)", "", k_h_wb)),
        a_h_pl = as.numeric(gsub(" \\(.*\\)", "", a_h_pl)),
        b_h_pl = as.numeric(gsub(" \\(.*\\)", "", b_h_pl)),
        a_dbh = as.numeric(gsub(" \\(.*\\)", "", a_dbh)),
        b_dbh = as.numeric(gsub(" \\(.*\\)", "", b_dbh)),
        c_dbh = as.numeric(gsub(" \\(.*\\)", "", c_dbh)),
        a_dbh1 = as.numeric(gsub(" \\(.*\\)", "", a_dbh1)),
        b_dbh1 = as.numeric(gsub(" \\(.*\\)", "", b_dbh1)),
        a_dbh2 = as.numeric(gsub(" \\(.*\\)", "", a_dbh2)),
        b_dbh2 = as.numeric(gsub(" \\(.*\\)", "", b_dbh2)),
        a_dbh3 = as.numeric(gsub(" \\(.*\\)", "", a_dbh3)),
        b_dbh3 = as.numeric(gsub(" \\(.*\\)", "", b_dbh3))
      ) |>
      mutate(
        log_AGB_bl = -3.048 + 2.111 * log_dbh + 0.552 * log_h,

        height_wb = a_h_wb * (1 - exp(-b_h_wb * dbh^k_h_wb)),
        log_AGB_wb_h = -3.048 + 2.111 * log_dbh + 0.552 * log(height_wb),

        height_pl = a_h_pl * dbh^b_h_pl,
        log_AGB_pl_h = -3.048 + 2.111 * log_dbh + 0.552 * log(height_pl),

        dbh_pl = a_dbh * (cr^b_dbh) * (h^c_dbh),
        log_AGB_pl_dbh = -3.048 + 2.111 * log(dbh_pl) + 0.552 * log_h,

        dbh_pl1 = a_dbh1 * (cr * h)^b_dbh1,
        log_AGB_pl_dbh1 = -3.048 + 2.111 * log(dbh_pl1) + 0.552 * log_h,

        dbh_pl2 = a_dbh2 * (cr^b_dbh2),
        log_AGB_pl_dbh2 = -3.048 + 2.111 * log(dbh_pl2) + 0.552 * log_h,

        dbh_pl3 = a_dbh3 * (h^b_dbh3),
        log_AGB_pl_dbh3 = -3.048 + 2.111 * log(dbh_pl3) + 0.552 * log_h
      )

    return(agb_df)
  }

  agb_ang <- calculate_agb(filtered_ang)
  agb_gym <- calculate_agb(filtered_gym)

  combine_metrics <- function(agb_df) {
    agb_df |>
      summarize(
        rmse_wb_h = sqrt(mean((log_AGB_bl - log_AGB_wb_h)^2, na.rm = TRUE)),
        bias_wb_h = mean(log_AGB_wb_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_h = sqrt(mean((log_AGB_bl - log_AGB_pl_h)^2, na.rm = TRUE)),
        bias_pl_h = mean(log_AGB_pl_h - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh)^2, na.rm = TRUE)),
        bias_pl_dbh = mean(log_AGB_pl_dbh - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh1 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh1)^2, na.rm = TRUE)),
        bias_pl_dbh1 = mean(log_AGB_pl_dbh1 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh2 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh2)^2, na.rm = TRUE)),
        bias_pl_dbh2 = mean(log_AGB_pl_dbh2 - log_AGB_bl, na.rm = TRUE),

        rmse_pl_dbh3 = sqrt(mean((log_AGB_bl - log_AGB_pl_dbh3)^2, na.rm = TRUE)),
        bias_pl_dbh3 = mean(log_AGB_pl_dbh3 - log_AGB_bl, na.rm = TRUE)
      )
  }

  metrics_ang4 <- combine_metrics(agb_ang)
  metrics_gym4 <- combine_metrics(agb_gym)

  # Add species and tree counts to metrics
  metrics_ang4$num_species <- filtered_ang$num_species
  metrics_ang4$num_trees <- filtered_ang$num_trees
  metrics_gym4$num_species <- filtered_gym$num_species
  metrics_gym4$num_trees <- filtered_gym$num_trees

  metrics_ang4$num_species <- format_number(metrics_ang4$num_species)
  metrics_ang4$num_trees <- format_number(metrics_ang4$num_trees)
  metrics_gym4$num_species <- format_number(metrics_gym4$num_species)
  metrics_gym4$num_trees <- format_number(metrics_gym4$num_trees)

  # Define the plotting function for Angiosperms and Gymnosperms
  plot_agb <- function(agb_df, metrics_df, model_colors, division_label) {
      create_plot <- function(data, baseline, predicted, rmse, bias, model_name, tag_label, x_axis_label = NULL) {
        ggplot(data, aes(x = exp(baseline), y = exp(predicted))) +
          geom_point(alpha = 0.5, color = "grey", size = 0.3) +
          geom_smooth(method = "lm", se = FALSE, color = model_colors[model_name], size = 0.5) +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black", size = 0.5) +
          scale_x_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
          scale_y_log10(
            breaks = c(0.1, 10, 1000, 1e5),  # Adjust breaks to use 10^6
            limits = c(NA, 1e5),  # Adjust limits
            labels = scales::trans_format("log10", scales::math_format(10^.x))  # Use scientific notation for labels
          ) +
            labs(tag = tag_label, 
                x = x_axis_label, 
                y = switch(model_name,
                            "h_wb" = expression(Predicted~AGB[paste(H, " - ", WB)]~"(kg)"),
                            "h_pl" = expression(Predicted~AGB[paste(H, " - ", PL)]~"(kg)"),
                            "dbh_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " + ", H)]~"(kg)"),
                            "dbh1_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR, " × ", H)]~"(kg)"),
                            "dbh2_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", CR)]~"(kg)"),
                            "dbh3_pl" = expression(Predicted~AGB[paste(DBH, " ~ ", H)]~"(kg)"))) +

          annotate("text", x = 0.1, y = 1e5, label = paste("RMSE =", format(round(rmse, 3), nsmall = 3), "Mg"), hjust = 0, size = 1.5) +
          annotate("text", x = 0.1, y = 3e4, label = paste("Bias =", format(round(bias * 100, 3), nsmall = 3), "%"), hjust = 0, size = 1.5) +
          my_theme()
      }

      model_colors <- c(
        "h_wb" = "#D55E00",
        "h_pl" = "#0072B2",
        "dbh_pl" = "#E69F00",
        "dbh1_pl" = "#F0E442",
        "dbh2_pl" = "#009E73",
        "dbh3_pl" = "#56B4E9"
      )

      # Create all 6 subplots for one division (Angiosperm or Gymnosperm)
      p1 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_wb_h, 
                        metrics_df$rmse_wb_h, metrics_df$bias_wb_h, "h_wb", "(a)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p2 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_h, 
                        metrics_df$rmse_pl_h, metrics_df$bias_pl_h, "h_pl", "(b)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p3 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh, 
                        metrics_df$rmse_pl_dbh, metrics_df$bias_pl_dbh, "dbh_pl", "(c)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p4 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh1, 
                        metrics_df$rmse_pl_dbh1, metrics_df$bias_pl_dbh1, "dbh1_pl", "(d)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p5 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh2, 
                        metrics_df$rmse_pl_dbh2, metrics_df$bias_pl_dbh2, "dbh2_pl", "(e)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      p6 <- create_plot(agb_df, agb_df$log_AGB_bl, agb_df$log_AGB_pl_dbh3, 
                        metrics_df$rmse_pl_dbh3, metrics_df$bias_pl_dbh3, "dbh3_pl", "(f)",
                        x_axis_label = expression(AGB[obs]~"(kg)"))

      return(list(p1, p2, p3, p4, p5, p6))
    }

    # Generate individual plots for Angiosperms and Gymnosperms
    ang_plots <- plot_agb(agb_ang, metrics_ang4, model_colors, "Angiosperm")
    gym_plots <- plot_agb(agb_gym, metrics_gym4, model_colors, "Gymnosperm")

    # Combine the plots into a 4x3 grid: 4 columns and 3 rows
    p <- (ang_plots[[1]] | ang_plots[[2]] |
          gym_plots[[1]] | gym_plots[[2]] |
          ang_plots[[3]] | ang_plots[[4]] |
          gym_plots[[3]] | gym_plots[[4]] |
          ang_plots[[5]] | ang_plots[[6]] |
          gym_plots[[5]] | gym_plots[[6]]) +
      plot_layout(ncol = 4, nrow = 3) +
      plot_annotation(tag_levels = list(c("(a)", "(b)", "(c)", "(d)", "(e)", "(f)",
                                          "(g)", "(h)", "(i)", "(j)", "(k)", "(l)"))) &
      my_theme() +
      theme(
        text = element_text(size = 5.5),
        axis.title = element_text(size = 6.5),
        axis.text = element_text(size = 5),
        legend.text = element_text(size = 5),
        legend.title = element_text(size = 5),
        plot.margin = margin(t = 3, r = 3, b = 3, l = 2, unit = "pt")
      )
    
    p <- ggdraw() +
      draw_plot(p, x = 0, y = 0, width = 1, height = 0.95) +
      draw_label("Angiosperm", x = 0.28, y = 0.965, fontface = 'bold', size = 10) +  
      draw_label("Gymnosperm", x = 0.78, y = 0.965, fontface = 'bold', size = 10)

  if (export_yaml) {
    yaml::write_yaml(list(metrics_ang4 = as.list(metrics_ang4), metrics_gym4 = as.list(metrics_gym4)), yaml_file)
  }

  return(list(p = p, metrics_ang4 = metrics_ang4, metrics_gym4 = metrics_gym4))
}



#================
#SUPPLEMENT
#================

#' @title Extract Wood Density References from TRY Version 6 Dataset
#'
process_try_data <- function(file_path) {
  # Load the data
  try_data <- read_delim(file_path, delim = "\t")
  
  # Clean and filter the data
  try_data_cleaned <- try_data |>
    dplyr::select(
      LastName, 
      FirstName, 
      DatasetID,
      Dataset,
      Reference, 
      sp = AccSpeciesName, 
      wd = StdValue, 
      unit = UnitName
    ) |>
    filter(!is.na(sp) & sp != "unknown" & !is.na(wd) & !is.na(unit) & unit == "g/cm3")
  
  # Calculate mean wood density for each species
  mean_density <- try_data_cleaned |>
    group_by(sp) |>
    summarize(mean_wd = mean(as.numeric(wd), na.rm = TRUE))
  
  # Join the mean density back to the cleaned data
  try_data_final <- try_data_cleaned |>
    left_join(mean_density, by = "sp") |>
    distinct(Dataset, .keep_all = TRUE) |>
    dplyr::select(LastName, FirstName, DatasetID, Dataset) |>
    dplyr::mutate(No. = row_number()) |>
    dplyr::select(No., everything())
  
  # Convert to tibble and return the final data
  try_data_selected <- as_tibble(try_data_final)
  return(try_data_selected)
}


#' @title Generate Combined Community and Species-Level Plot for Tree Height and Crown Radius allometries
#' @description
#' This function generates a combined plot for tree height (H) and crown radius (CR) at both community and species levels,
#' using species-specific DBH ranges. It computes fitted values based on beta parameters extracted from model summaries
#' and plots these values along with community-level fitted lines.
#' @param stan_data_nlr_h Data for tree height (H).
#' @param stan_data_nlr_cr Data for crown radius (CR).
#' @param fit_nlr_summary_weibull1_h Model summary for tree height (H) (Weibull).
#' @param fit_nlr_summary_weibull1_cr Model summary for crown radius (CR) (Weibull).
#' @param length_out The number of points in the DBH sequence used to compute fitted values. Default is 100.
#' @return A combined ggplot2 object with plots for tree height and crown radius.
#' @export
generate_com_sp_plot <- function(stan_data_nlr_h, stan_data_nlr_cr, 
                                 fit_nlr_summary_weibull1_h, fit_nlr_summary_weibull1_cr,
                                 length_out = 100) {

  # Calculate the range of DBH for each species
  dbh_data <- data.frame(sp = as.factor(stan_data_nlr_h$jj), DBH = stan_data_nlr_h$x)
  dbh_range <- dbh_data |>
    group_by(sp) |>
    summarize(min_dbh = min(DBH), max_dbh = max(DBH))

  # Function to extract beta parameters and compute fitted values using species-specific DBH range
  extract_and_compute_fitted <- function(data, summary, dbh_range) {
    log_y <- data$log_y
    sp <- as.factor(data$jj)
    
    # Extract beta parameters
    beta <- summary |> filter(grepl("beta\\[", variable))
    log_a <- beta |> filter(str_detect(variable, "beta\\[\\d+,1\\]")) |> pull(q50)
    b <- beta |> filter(str_detect(variable, "beta\\[\\d+,2\\]")) |> pull(q50)
    k <- beta |> filter(str_detect(variable, "beta\\[\\d+,3\\]")) |> pull(q50)
    
    # Create DBH sequence for each species based on actual range
    create_dbh_sequence <- function(species) {
      min_dbh <- dbh_range$min_dbh[dbh_range$sp == species]
      max_dbh <- dbh_range$max_dbh[dbh_range$sp == species]
      seq(min_dbh, max_dbh, length.out = length_out)
    }
    
    # Calculate fitted values for each species using species-specific DBH range
    fitted_list <- lapply(unique(sp), function(s) {
      dbh_seq <- create_dbh_sequence(s)
      i <- which(unique(sp) == s)  # Get the index of the species
      fitted_log_y <- log_a[i] + log(1 - exp(-b[i] * (dbh_seq ^ k[i])))
      fitted_values <- exp(fitted_log_y)
      
      data.frame(
        DBH = dbh_seq,
        value = fitted_values,
        variable = paste0("fitted_sp_", i)
      )
    })
    
    # Combine all species data frames into one
    data_long <- do.call(rbind, fitted_list)
    
    # Extract community-level parameters
    gamma <- summary |> filter(grepl("gamma", variable))
    community_log_a_hat <- gamma |> filter(variable == "gamma[1,1]") |> pull(q50)
    community_b_hat <- gamma |> filter(variable == "gamma[1,2]") |> pull(q50)
    community_k_hat <- gamma |> filter(variable == "gamma[1,3]") |> pull(q50)
    
    community_log_a <- exp(community_log_a_hat)
    community_b <- exp(community_b_hat)
    community_k <- exp(community_k_hat)
    
    community_fitted_log_y <- community_log_a + log(1 - exp(-community_b * (data$x ^ community_k)))
    community_fitted <- exp(community_fitted_log_y)
    
    list(fitted_list = fitted_list, community_fitted = community_fitted, log_y = log_y, x = data$x, sp = sp, data_long = data_long)
  }

  # Apply the extraction function to both H and CR data
  h_params <- extract_and_compute_fitted(stan_data_nlr_h, fit_nlr_summary_weibull1_h, dbh_range)
  cr_params <- extract_and_compute_fitted(stan_data_nlr_cr, fit_nlr_summary_weibull1_cr, dbh_range)

  # Create the plots
  p1 <- ggplot(data.frame(H = exp(h_params$log_y), DBH = h_params$x, sp = h_params$sp), aes(x = DBH, y = H)) +
    geom_line(data = h_params$data_long, aes(x = DBH, y = value, group = variable), color = "#72b6e3", alpha = 0.3, linewidth = 0.15) +
    geom_line(aes(y = h_params$community_fitted), color = "#0f92e9", linewidth = 0.6) +
    labs(x = "DBH (cm)", y = "Tree Height (m)") +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_log10(labels = scales::label_number()) +
    theme_minimal() +
    theme(
    axis.title.x = element_blank(),   # Remove the x-axis title
    axis.text.x = element_blank()     # Remove the x-axis text labels
  )
  
  p2 <- ggplot(data.frame(CR = exp(cr_params$log_y), DBH = cr_params$x, sp = cr_params$sp), aes(x = DBH, y = CR)) +
    geom_line(data = cr_params$data_long, aes(x = DBH, y = value, group = variable), color = "#72b6e3", alpha = 0.3, linewidth = 0.15) +
    geom_line(aes(y = cr_params$community_fitted), color = "#0f92e9", linewidth = 0.6) +
    labs(x = "DBH (cm)", y = "Crown Radius (m)") +
    scale_x_log10(labels = scales::label_number()) +
    scale_y_log10(labels = scales::label_number()) +
    theme_minimal()
  
  # Combine and return the plots
  p <- p1 / p2 + 
    plot_layout(ncol = 1) &  
    theme(
      plot.margin = unit(c(10, 10, 10, 10), "pt"),  
      text = element_text(size = 9),               
      legend.title = element_text(size = 10)       
    )
  
  return(p)
}

#=======================
# SPECIES-LEVEL PLOTS 
#=======================
generate_h_sp_plot <- function(
  tallo_reduced_lr_df_ang_h,
  sp_posterior_h_df,
  num_species = 32) {
  
  h_data <- tallo_reduced_lr_df_ang_h |> 
    filter(!is.na(h))
  
  # Sample a specific number of species
  species_names <- sample(unique(h_data$sp), min(num_species, n_distinct(h_data$sp)))
  
  # Split the dataset by species
  sub_datasets <- h_data |>
    filter(sp %in% species_names) |>
    group_split(sp)
  
  # Assign species names as dataset names
  names(sub_datasets) <- species_names
  
  # Calculate the range of DBH for each species
  dbh_range <- h_data |>
    group_by(sp) |>
    summarize(
      min_dbh = min(dbh, na.rm = TRUE),
      max_dbh = max(dbh, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Function to create DBH sequence for a given species
  create_dbh_sequence <- function(species, length_out = 200) {
    min_dbh <- dbh_range$min_dbh[dbh_range$sp == species]
    max_dbh <- dbh_range$max_dbh[dbh_range$sp == species]
    # seq(min_dbh * 0.5, max_dbh * 1.5, length.out = length_out)
    seq(min_dbh, max_dbh, length.out = length_out)
  }
  
  combined_data <- data.frame()
  
  for (i in seq_along(sub_datasets)) {
    species_data <- sub_datasets[[i]]
    species_name <- names(sub_datasets)[i]
    
    # Extract posterior data for the current species
    species_post <- sp_posterior_h_df |>
      filter(Species == species_name)
    
    if (nrow(species_post) == 0) {
      warning(paste("No posterior data for species:", species_name))
      next
    }
    
    # Extract Weibull parameters
    a <- as.numeric(sub(" \\(.*\\)", "", species_post$a))
    log_a <- log(a)
    b <- as.numeric(sub(" \\(.*\\)", "", species_post$b))
    k <- as.numeric(sub(" \\(.*\\)", "", species_post$k))
    
    # Generate DBH sequence for the current species
    extended_dbh <- create_dbh_sequence(species_name)
    
    # Compute the fitted height values for the extended range
    fitted_log_y <- log_a + log(1 - exp(-b * (extended_dbh ^ k)))
    fitted_height <- exp(fitted_log_y)
    
    # Create a dataframe for the fitted curve
    fitted_curve <- data.frame(
      DBH = extended_dbh,
      H = fitted_height,
      sp = species_name,
      Source = "Fitted"
    )
    
    # Combine observed and fitted data
    combined_species_data <- species_data |>
      mutate(Species = species_name, Source = "Observed") |>
      dplyr::select(DBH = dbh, H = h, sp, Source) |>
      bind_rows(fitted_curve)
    
    # Append to the combined dataset
    combined_data <- bind_rows(combined_data, combined_species_data)
  }
  
  # Plot the combined dataset
  p <- ggplot(combined_data, aes(x = DBH, y = H, color = Source)) +
    geom_point(data = combined_data |> filter(Source == "Observed"), alpha = 0.5, color = "gray") +
    geom_line(data = combined_data |> filter(Source == "Fitted"), aes(group = sp), color = "#0f92e9", linewidth = 0.6) +
    facet_wrap(~ sp, scales = "free", ncol = 4) +
    labs(
      x = "DBH (cm)",
      y = "Tree Height (m)"
    ) +
    # scale_x_log10() +
    # scale_y_log10() + 
    my_theme()
  
  return(p)
}

#=======================================
# SPECIES-LEVEL PLOTS: H-DBH LARGE TREES
#=======================================
generate_h_sp_large_plot <- function(
  tallo_reduced_lr_df_ang_h,
  sp_posterior_h_df,
  num_species = 32) {
  
  # Filter data to ensure H values are not missing and only include species with DBH > 200
  h_data <- tallo_reduced_lr_df_ang_h |> 
    filter(!is.na(h)) |>
    group_by(sp) |>
    filter(any(dbh > 200)) |>
    ungroup()
  
  # Sample a specific number of species
  species_names <- sample(unique(h_data$sp), min(num_species, n_distinct(h_data$sp)))
  
  # Split the dataset by species
  sub_datasets <- h_data |>
    filter(sp %in% species_names) |>
    group_split(sp)
  
  # Assign species names as dataset names
  names(sub_datasets) <- species_names
  
  # Calculate the range of DBH for each species
  dbh_range <- h_data |>
    group_by(sp) |>
    summarize(
      min_dbh = min(dbh, na.rm = TRUE),
      max_dbh = max(dbh, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Function to create DBH sequence for a given species
  create_dbh_sequence <- function(species, length_out = 200) {
    min_dbh <- dbh_range$min_dbh[dbh_range$sp == species]
    max_dbh <- dbh_range$max_dbh[dbh_range$sp == species]
    # seq(min_dbh * 0.5, max_dbh * 1.5, length.out = length_out)
    seq(min_dbh, max_dbh, length.out = length_out)
  }
  
  combined_data <- data.frame()
  
  for (i in seq_along(sub_datasets)) {
    species_data <- sub_datasets[[i]]
    species_name <- names(sub_datasets)[i]
    
    # Extract posterior data for the current species
    species_post <- sp_posterior_h_df |>
      filter(Species == species_name)
    
    if (nrow(species_post) == 0) {
      warning(paste("No posterior data for species:", species_name))
      next
    }
    
    # Extract Weibull parameters
    a <- as.numeric(sub(" \\(.*\\)", "", species_post$a))
    log_a <- log(a)
    b <- as.numeric(sub(" \\(.*\\)", "", species_post$b))
    k <- as.numeric(sub(" \\(.*\\)", "", species_post$k))
    
    # Generate DBH sequence for the current species
    extended_dbh <- create_dbh_sequence(species_name)
    
    # Compute the fitted height values for the extended range
    fitted_log_y <- log_a + log(1 - exp(-b * (extended_dbh ^ k)))
    fitted_height <- exp(fitted_log_y)
    
    # Create a dataframe for the fitted curve
    fitted_curve <- data.frame(
      DBH = extended_dbh,
      H = fitted_height,
      sp = species_name,
      Source = "Fitted"
    )
    
    # Combine observed and fitted data
    combined_species_data <- species_data |>
      mutate(Species = species_name, Source = "Observed") |>
      dplyr::select(DBH = dbh, H = h, sp, Source) |>
      bind_rows(fitted_curve)
    
    # Append to the combined dataset
    combined_data <- bind_rows(combined_data, combined_species_data)
  }
  
  # Plot the combined dataset
  p <- ggplot(combined_data, aes(x = DBH, y = H, color = Source)) +
    geom_point(data = combined_data |> filter(Source == "Observed"), alpha = 0.5, color = "gray") +
    geom_line(data = combined_data |> filter(Source == "Fitted"), aes(group = sp), color = "#0f92e9", linewidth = 0.6) +
    facet_wrap(~ sp, scales = "free", ncol = 4) +
    labs(
      x = "DBH (cm)",
      y = "Tree Height (m)"
    ) +
    # scale_x_log10() +
    # scale_y_log10() + 
    my_theme()
  
  return(p)
}
