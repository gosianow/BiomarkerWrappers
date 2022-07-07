


# colors_line = NULL; 
# variable_names = NULL; 
# xlab = NULL; ylab = NULL; title = NULL; subtitle = NULL;
# legend_colors_line_title = NULL; legend_position = "right"; facet_label_both = TRUE; 
# line_size = 1; line_type = 1; line_alpha = 1;
# smooth_method = "lm"; smooth_formula = y ~ x; smooth_se = FALSE;
# smooth_size = 2; smooth_linetype = 1; 
# point_size = 1.5; point_alpha = 1; 
# title_size = NULL; strip_text_size = NULL; facet_scales = "fixed"; xlim = NULL; ylim = NULL; 
# background_grid_major = "none"





#' Line plot
#' 
#' @param data Data frame.
#' @export
wrapper_line_plot_core <- function(data, x_var, y_var, group_var, color_line_var = NULL, shape_point_var = NULL, 
  facet_var = NULL, 
  colors_line = NULL, shapes_point = NULL,
  variable_names = NULL, 
  title = TRUE, subtitle = TRUE, xlab = TRUE, ylab = TRUE,
  legend_colors_line_title = TRUE, legend_shapes_point_title = TRUE, legend_position = "right", facet_label_both = TRUE, 
  line_size = 1, line_type = 1, line_alpha = 1, 
  smooth = "none", smooth_method = "lm", smooth_formula = y ~ x, smooth_se = FALSE,
  smooth_size = 2, smooth_linetype = 1, 
  point_size = 1.5, point_alpha = 1, 
  title_size = NULL, strip_text_size = NULL, facet_scales = "fixed", xlim = NULL, ylim = NULL,
  scale_y_continuous_custome = scale_y_continuous(), 
  axis_text_x_angle = 0, axis_text_x_vjust = 1, axis_text_x_hjust = 0.5,
  background_grid_major = "none"){
  
  
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
  
  if(!is.null(shape_point_var)){
    stopifnot(length(shape_point_var) == 1)
    stopifnot(is.factor(data[, shape_point_var]))
  }
  
  if(!is.null(facet_var)){
    stopifnot(length(facet_var) == 1)
    stopifnot(is.factor(data[, facet_var]))
  }
  
  stopifnot(length(smooth) == 1)
  stopifnot(smooth %in% c("none", "pooled", "strat"))
  
  ### Keep non-missing data
  data <- data[stats::complete.cases(data[, c(x_var, y_var, facet_var)]), , drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Colors for lines
  # -------------------------------------------------------------------------
  
  if(is.null(color_line_var)){
    
    ### Create a dummy variable 
    stopifnot(!"color_line_dummy" %in% colnames(data))
    data[, "color_line_dummy"] <- factor("color_line_dummy")
    color_line_var <- "color_line_dummy"
    
    if(is.null(colors_line)){
      colors_line <- "grey40"
    }else{
      colors_line <- colors_line[1]
    }
  }else{
    colors_line <- format_colors(levels(data[, color_line_var]), colors = colors_line)
  }
  
  
  # -------------------------------------------------------------------------
  # Shapes for points
  # -------------------------------------------------------------------------
  
  if(is.null(shape_point_var)){
    
    ### Create a dummy variable 
    stopifnot(!"shape_point_dummy" %in% colnames(data))
    data[, "shape_point_dummy"] <- factor("shape_point_dummy")
    shape_point_var <- "shape_point_dummy"
    
    if(is.null(shapes_point)){
      shapes_point <- 16
    }else{
      shapes_point <- shapes_point[1]
    }
  }else{
    shapes_point <- format_shapes(levels(data[, shape_point_var]), shapes = shapes_point)
  }
  
  
  # -------------------------------------------------------------------------
  # Axis and legend labels
  # -------------------------------------------------------------------------
  
  
  if(is.logical(title)){
    title <- NULL
  }
  
  if(is.logical(subtitle)){
    subtitle <- NULL
  }
  
  if(is.logical(xlab)){
    if(xlab){
      xlab <- variable_names[x_var]
    }else{
      xlab <- NULL
    }
  }
  
  if(is.logical(ylab)){
    if(ylab){
      ylab <- variable_names[y_var]
    }else{
      ylab <- NULL
    }
  }
  
  if(is.logical(legend_colors_line_title)){
    if(legend_colors_line_title){
      legend_colors_line_title <- variable_names[color_line_var]
    }else{
      legend_colors_line_title <- NULL
    }
  }
  
  if(is.logical(legend_shapes_point_title)){
    if(legend_shapes_point_title){
      legend_shapes_point_title <- variable_names[shape_point_var]
    }else{
      legend_shapes_point_title <- NULL
    }
  }
  
  # -------------------------------------------------------------------------
  # Make the plot
  # -------------------------------------------------------------------------
  
  legend_show_colors_line <- color_line_var != "color_line_dummy"
  legend_show_shapes_point <- shape_point_var != "shape_point_dummy"
  
  
  ggpl <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]], group = .data[[group_var]])) +
    geom_line(aes(color = .data[[color_line_var]]), linetype = line_type, size = line_size, alpha = line_alpha) +
    scale_color_manual(name = legend_colors_line_title, values = colors_line, drop = FALSE, na.value = "grey") +
    guides(color = ifelse(legend_show_colors_line, guide_legend(), "none"))
  
  
  if(all(shapes_point %in% 21:25)){
    
    ggpl <- ggpl +
      geom_point(aes(fill = .data[[color_line_var]], shape = .data[[shape_point_var]]), size = point_size, alpha = point_alpha) +
      scale_fill_manual(name = legend_colors_line_title, values = colors_line, drop = FALSE, na.value = "grey") +
      scale_shape_manual(name = legend_shapes_point_title, values = shapes_point) +
      guides(fill = ifelse(legend_show_colors_line, guide_legend(), "none"), shape = ifelse(legend_show_shapes_point, guide_legend(), "none"))
    
    
  }else{
    
    ggpl <- ggpl +
      geom_point(aes(color = .data[[color_line_var]], shape = .data[[shape_point_var]]), size = point_size, alpha = point_alpha) +
      scale_shape_manual(name = legend_shapes_point_title, values = shapes_point) +
      guides(shape = ifelse(legend_show_shapes_point, guide_legend(), "none"))
    
  }
  
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle) + 
    ylab(ylab) +
    xlab(xlab) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = title_size),
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),
      panel.border = element_rect(colour = "black", size = 0.8),
      legend.position = legend_position) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.15) +
    coord_cartesian(xlim = xlim, ylim = ylim) +
    scale_y_continuous_custome
  
  
  if(is.factor(data[, x_var])){
    
    ggpl <- ggpl + 
      scale_x_discrete(drop = FALSE)
    
  }
  
  
  
  if(smooth == "pooled"){
    
    if(is.factor(data[, x_var])){
      
      ggpl <- ggpl + 
        stat_summary(aes(group = 1), 
          geom = "line", fun = "mean", color = "dodgerblue", linetype = smooth_linetype, size = smooth_size)
      
    }else{
      
      ggpl <- ggpl + 
        geom_smooth(aes(group = 1),
          method = smooth_method, formula = smooth_formula, se = smooth_se, 
          linetype = smooth_linetype, size = smooth_size)
      
    }
    
  }else if(smooth == "strat"){
    
    if(is.factor(data[, x_var])){
      
      ggpl <- ggpl + 
        stat_summary(aes(group = .data[[color_line_var]], color = .data[[color_line_var]]), 
          geom = "line", fun = "mean", color = "dodgerblue", linetype = smooth_linetype, size = smooth_size)
      
    }else{
      
      ggpl <- ggpl + 
        geom_smooth(aes(group = .data[[color_line_var]], color = .data[[color_line_var]]), 
          method = smooth_method, formula = smooth_formula, se = smooth_se, 
          linetype = smooth_linetype, size = smooth_size)
      
    }
    
  }
  
  
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
      facet_wrap(stats::as.formula(paste("~", facet_var)), labeller = labeller, scales = facet_scales) +
      theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = strip_text_size))
    
  }
  
  
  ggpl
  
  
  
}









#' @rdname wrapper_line_plot_core
#' @param strat1_var Name of the first stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' @export
wrapper_line_plot_core_strat <- function(data, x_var, y_var, group_var, color_line_var = NULL, shape_point_var = NULL, facet_var = NULL, 
  strat1_var = NULL, strat2_var = NULL, 
  colors_line = NULL, shapes_point = NULL,
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, ylab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_line_title = TRUE, legend_shapes_point_title = TRUE, legend_position = "right", facet_label_both = TRUE, 
  line_size = 1, line_type = 1, line_alpha = 1,
  smooth = "none", smooth_method = "lm", smooth_formula = y ~ x, smooth_se = FALSE,
  smooth_size = 2, smooth_linetype = 1, 
  point_size = 1.5, point_alpha = 1, 
  title_size = NULL, strip_text_size = NULL, facet_scales = "fixed", xlim = NULL, ylim = NULL, 
  scale_y_continuous_custome = scale_y_continuous(), 
  axis_text_x_angle = 0, axis_text_x_vjust = 1, axis_text_x_hjust = 0.5,
  background_grid_major = "none", 
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1, less_legends = FALSE){
  
  
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
  
  data <- data[stats::complete.cases(data[, c(x_var, y_var, strat1_var, strat2_var)]), , drop = FALSE]
  
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
    
    data_strata2 <- data[data[, strat2_var] %in% strata2_levels[j], ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    ### subtitle_strat2
    if(strat2_var == "strat2_dummy"){
      subtitle_strat2 <- NULL
    }else{
      if(strat2_label_both){
        subtitle_strat2 <- paste0(variable_names[strat2_var], ": ", strata2_levels[j])
      }else{
        subtitle_strat2 <- strata2_levels[j]
      }
    }
    
    
    ggpl <- lapply(1:length(strata1_levels), function(i){
      # i = 1
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] %in% strata1_levels[i], ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      ### subtitle_strat1
      if(strat1_var == "strat1_dummy"){
        subtitle_strat1 <- NULL
      }else{
        if(strat1_label_both){
          subtitle_strat1 <- paste0(variable_names[strat1_var], ": ", strata1_levels[i])
        }else{
          subtitle_strat1 <- strata1_levels[i]
        }
      }
      
      subtitle <- paste0(c(subtitle_strat2, subtitle_strat1), collapse = "\n")
      if(subtitle == ""){
        subtitle <- NULL
      }
      
      
      ggpl <- wrapper_line_plot_core(data = data_strata1, x_var = x_var, y_var = y_var, group_var = group_var, color_line_var = color_line_var, shape_point_var = shape_point_var, facet_var = facet_var, 
        colors_line = colors_line, shapes_point = shapes_point, 
        variable_names = variable_names, 
        xlab = xlab, ylab = ylab, title = title, subtitle = subtitle,
        legend_colors_line_title = legend_colors_line_title, legend_shapes_point_title = legend_shapes_point_title, legend_position = legend_position, facet_label_both = facet_label_both, 
        line_size = line_size, line_type = line_type, line_alpha = line_alpha,
        smooth = smooth, smooth_method = smooth_method, smooth_formula = smooth_formula, smooth_se = smooth_se,
        smooth_size = smooth_size, smooth_linetype = smooth_linetype, 
        point_size = point_size, point_alpha = point_alpha, 
        title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, xlim = xlim, ylim = ylim,
        scale_y_continuous_custome = scale_y_continuous_custome,
        background_grid_major = background_grid_major,
        axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust)
      
      
      return(ggpl)
      
      
    })
    
    if(less_legends && legend_position == "right"){
      
      ggpl_non_empty <- ggpl[!sapply(ggpl, is.null)]
      
      # Extract the legend from one of the plots
      legend <- get_legend(
        # Create some space to the left of the legend
        ggpl_non_empty[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
      )
      
      # Remove legends in the plots
      ggpl <- lapply(ggpl, function(x) x + theme(legend.position = "none"))
      
      # Append the legend panel 
      ggpl[["legend"]] <- legend
      
    }
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  ggpl
  
  
}

















