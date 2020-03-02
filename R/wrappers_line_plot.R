






# colors_line = NULL
# variable_names = NULL
# xlab = NULL
# ylab = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# legend_title_colors_line = NULL
# facet_label_both = TRUE
# line_size = 1
# line_type = 1
# point_size = 1.5
# point_shape = 1
# point_alpha = 1
# title_size = 12
# strip_text_size = NULL
# strip_position = "top"
# facet_scales = "fixed"
# xlim = NULL
# ylim = NULL
# background_grid_major = "none"








#' Line plot
#' 
#' Generate a signle line plot.
#' 
#' @param data Data frame.
wrapper_core_line_plot <- function(data, x_var, y_var, group_var, color_line_var = NULL, facet_var = NULL, colors_line = NULL, variable_names = NULL, xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, tag = NULL, legend_title_colors_line = NULL, facet_label_both = TRUE, line_size = 1, line_type = 1, point_size = 1.5, point_shape = 1, point_alpha = 1, title_size = 12, strip_text_size = NULL, strip_position = "top", facet_scales = "fixed", xlim = NULL, ylim = NULL, background_grid_major = "none"){
  
  
  # -------------------------------------------------------------------------
  # Checks about data
  # -------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  stopifnot(length(group_var) == 1)
  # stopifnot(is.factor(data[, group_var]) || is.character(data[, group_var]))
  
  if(!is.null(color_line_var)){
    stopifnot(length(color_line_var) == 1)
    stopifnot(is.factor(data[, color_line_var]))
  }
  
  if(!is.null(facet_var)){
    stopifnot(length(facet_var) == 1)
    stopifnot(is.factor(data[, facet_var]))
  }
  
  
  ### Keep non-missing data
  
  ## We do not exclude the missing values in color_line_var
  data <- data[complete.cases(data[, c(x_var, y_var, facet_var)]), , drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Colors
  # -------------------------------------------------------------------------
  
  if(is.null(color_line_var)){
    
    ### Create a dummy variable 
    stopifnot(!"color_line_dummy" %in% colnames(data))
    data[, "color_line_dummy"] <- factor("color_line_dummy")
    color_line_var <- "color_line_dummy"
    
    if(is.null(colors_line)){
      colors_line <- "darkgrey"
    }else{
      colors_line <- colors_line[1]
    }
  }else{
    colors_line <- format_colors(levels = levels(data[, color_line_var]), colors = colors_line)
  }
  
  
  # -------------------------------------------------------------------------
  # Labels
  # -------------------------------------------------------------------------
  
  
  if(is.null(xlab)){
    xlab <- variable_names[x_var]
  }
  if(is.null(ylab)){
    ylab <- variable_names[y_var]
  }
  if(is.null(legend_title_colors_line)){
    legend_title_colors_line <- variable_names[color_line_var]
  }
  
  
  
  # -------------------------------------------------------------------------
  # Make the plot
  # -------------------------------------------------------------------------
  
  legend_show_colors_line <- color_line_var != "color_line_dummy"
  
  
  
  ggpl <- ggplot(data, aes_string(x = x_var, y = y_var, group = group_var)) +
    geom_line(aes_string(color = color_line_var), linetype = line_type, size = line_size, show.legend = legend_show_colors_line) +
    scale_color_manual(name = legend_title_colors_line, values = colors_line, drop = FALSE, na.value = "grey")
  
  
  if(point_shape %in% 21:25){
    
    ggpl <- ggpl +
      geom_point(aes_string(fill = color_line_var), size = point_size, shape = point_shape, alpha = point_alpha, show.legend = legend_show_colors_line) +
      scale_fill_manual(name = legend_title_colors_line, values = colors_line, drop = FALSE, na.value = "grey")
    
    
  }else{
    
    ggpl <- ggpl +
      geom_point(aes_string(color = color_line_var), size = point_size, shape = point_shape, alpha = point_alpha, show.legend = legend_show_colors_line) 
    
  }
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle, tag = tag) + 
    ylab(ylab) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = title_size),
      plot.tag.position = "top",
      plot.tag = element_text(size = title_size, face = "plain")) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
    coord_cartesian(xlim = xlim, ylim = ylim)
  
  
  if(!is.null(facet_var)){
    
    if(facet_label_both){
      labeller <- function(labels, multi_line = FALSE, sep = ": "){
        colnames(labels) <- variable_names[colnames(labels)]
        out <- ggplot2::label_both(labels, multi_line = multi_line, sep = sep)
        out
      }
    }else{
      labeller <- "label_value"
    }
    
    ggpl <- ggpl +
      facet_wrap(as.formula(paste("~", facet_var)), strip.position = strip_position, labeller = labeller, scales = facet_scales) +
      theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = strip_text_size),
        strip.placement = "outside", 
        axis.line = element_blank()) +
      panel_border(colour = "black", linetype = 1, size = 1, remove = FALSE)
    
  }
  
  
  ggpl
  
  
  
}








#' Line plot
#' 
#' Generate line plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_line_plot_strat <- function(data, x_var, y_var, group_var, color_line_var = NULL, facet_var = NULL, strat1_var = NULL, strat2_var = NULL, colors_line = NULL, variable_names = NULL, xlab = NULL, ylab = NULL, title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, legend_title_colors_line = NULL, facet_label_both = TRUE, line_size = 1, line_type = 1, point_size = 1.5, point_shape = 1, point_alpha = 1, title_size = 12, strip_text_size = NULL, strip_position = "top", facet_scales = "fixed", xlim = NULL, ylim = NULL, background_grid_major = "none", strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  if(!is.null(strat1_var)){
    stopifnot(length(strat1_var) == 1)
    stopifnot(is.factor(data[, strat1_var]))
  }else{
    ### Add dummy variable to data
    stopifnot(!"strat1_dummy" %in% colnames(data))
    data[, "strat1_dummy"] <- factor("strat1_dummy")
    strat1_var <- "strat1_dummy"
  }
  
  if(!is.null(strat2_var)){
    stopifnot(length(strat2_var) == 1)
    stopifnot(is.factor(data[, strat2_var]))
  }else{
    ### Add dummy variable to data
    stopifnot(!"strat2_dummy" %in% colnames(data))
    data[, "strat2_dummy"] <- factor("strat2_dummy")
    strat2_var <- "strat2_dummy"
  }
  
  ### Keep non-missing data
  
  data <- data[complete.cases(data[, c(x_var, y_var, strat1_var, strat2_var)]), , drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  # -------------------------------------------------------------------------
  # Scales, xlim, ylim
  # -------------------------------------------------------------------------
  
  if(strat_scales == "fixed" || strat_scales == "free_y"){
    if(is.null(xlim)){
      if(!is.factor(data[, x_var])){
        xlim <- range(data[, x_var])
      }
    }
  }
  if(strat_scales == "fixed" || strat_scales == "free_x"){
    if(is.null(ylim)){
      ylim <- range(data[, y_var])
    }
  }
  
  # -------------------------------------------------------------------------
  # Lapply to make the stratified plots
  # -------------------------------------------------------------------------
  
  
  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  ggpl <- lapply(1:length(strata2_levels), function(j){
    # j = 1
    
    data_strata2 <- data[data[, strat2_var] == strata2_levels[j] & !is.na(data[, strat2_var]), ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    if(strat2_var == "strat2_dummy"){
      tag <- NULL
    }else{
      if(tag_label_both){
        tag <- paste0(variable_names[strat2_var], ": ", strata2_levels[j])
      }else{
        tag <- strata2_levels[j]
      }
    }
    
    
    ggpl <- lapply(1:length(strata1_levels), function(i){
      # i = 1
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      if(strat1_var == "strat1_dummy"){
        subtitle <- NULL
      }else{
        if(subtitle_label_both){
          subtitle <- paste0(variable_names[strat1_var], ": ", strata1_levels[i])
        }else{
          subtitle <- strata1_levels[i]
        }
      }
      
      
      ggpl <- wrapper_core_line_plot(data = data_strata1, x_var = x_var, y_var = y_var, group_var = group_var, color_line_var = color_line_var, facet_var = facet_var, colors_line = colors_line, variable_names = variable_names, xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, tag = tag, legend_title_colors_line = legend_title_colors_line, facet_label_both = facet_label_both, line_size = line_size, line_type = line_type, point_size = point_size, point_shape = point_shape, point_alpha = point_alpha, title_size = title_size, strip_text_size = strip_text_size, strip_position = strip_position, facet_scales = facet_scales, xlim = xlim, ylim = ylim, background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  ggpl
  
  
}

















