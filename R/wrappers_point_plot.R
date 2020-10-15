







#' Scatter plot
#' 
#' @param data Data frame.
#' @examples 
#' 
#' data(bdata)
#' data <- bdata
#' 
#' x_var <- "GeneA"
#' y_var <- "GeneB"
#' 
#' 
#' wrapper_point_plot_core(data = data, x_var = x_var, y_var = y_var)
#' 
#' wrapper_point_plot_core(data = bdata, x_var = x_var, y_var = y_var, color_point_var = "Sex", smooth = "strat")
#' 
#' wrapper_point_plot_core(data = bdata, x_var = x_var, y_var = y_var, color_point_var = "Age")
#' 
#' @export
wrapper_point_plot_core <- function(data, x_var, y_var, color_point_var = NULL, facet_var = NULL, 
  colors_point = NULL, scale_gradient = "gradientn", color_low_point = '#42399B', color_mid_point = "white", color_high_point = '#D70131', midpoint = 0,
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, 
  legend_colors_point_title = NULL, legend_position = "right", facet_label_both = TRUE, 
  point_size = 1.5, point_shape = 20, point_alpha = 1, point_stroke = 0.8,
  smooth = "none", smooth_method = "auto", smooth_formula = y ~ x, smooth_se = FALSE,
  smooth_size = 2, smooth_linetype = 1, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", xlim = NULL, ylim = NULL, 
  background_grid_major = "none"){
  
  
  # -------------------------------------------------------------------------
  # Checks about data
  # -------------------------------------------------------------------------
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.numeric(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  if(!is.null(color_point_var)){
    stopifnot(length(color_point_var) == 1)
    stopifnot(is.factor(data[, color_point_var]) || is.numeric(data[, color_point_var]))
  }
  
  if(!is.null(facet_var)){
    stopifnot(length(facet_var) == 1)
    stopifnot(is.factor(data[, facet_var]))
  }
  
  stopifnot(length(smooth) == 1)
  stopifnot(smooth %in% c("none", "pooled", "strat"))
  
  stopifnot(length(scale_gradient) == 1)
  stopifnot(scale_gradient %in% c("gradientn", "gradient2", "gradient"))
  
  
  ### Keep non-missing data
  
  ## We do not exclude the missing values in color_point_var
  data <- data[stats::complete.cases(data[, c(x_var, y_var, facet_var)]), , drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Colors
  # -------------------------------------------------------------------------
  
  if(is.null(color_point_var)){
    
    ### Create a dummy variable 
    stopifnot(!"color_point_dummy" %in% colnames(data))
    data[, "color_point_dummy"] <- factor("color_point_dummy")
    color_point_var <- "color_point_dummy"
    
    if(is.null(colors_point)){
      colors_point <- "grey40"
    }else{
      colors_point <- colors_point[1]
    }
    
  }else{
    
    if(is.factor(data[, color_point_var])){
      colors_point <- format_colors(levels = levels(data[, color_point_var]), colors = colors_point)
    }else{
      colors_point <- rev(RColorBrewer::brewer.pal(11, "Spectral"))
    }
    
  }
  
  if(is.factor(data[, color_point_var])){
    scale <- "manual"
  }else{
    scale <- scale_gradient
    if(smooth == "strat"){
      smooth <- "pooled"
    }
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
  if(is.null(legend_colors_point_title)){
    legend_colors_point_title <- variable_names[color_point_var]
  }
  
  
  
  # -------------------------------------------------------------------------
  # Make the plot
  # -------------------------------------------------------------------------
  
  legend_show_colors_point <- color_point_var != "color_point_dummy"
  
  
  ggpl <- ggplot(data, aes(x = .data[[x_var]], y = .data[[y_var]])) 
  
  
  if(point_shape %in% 21:25){
    
    ggpl <- ggpl +
      geom_point(aes(fill = .data[[color_point_var]]), size = point_size, shape = point_shape, alpha = point_alpha, stroke = point_stroke, 
        show.legend = legend_show_colors_point) 
    
    if(scale == "manual"){
      ggpl <- ggpl + scale_fill_manual(name = legend_colors_point_title, values = colors_point, drop = FALSE)
    }else if(scale == "gradientn") {
      ggpl <- ggpl + scale_fill_gradientn(name = legend_colors_point_title, colors = colors_point)
    }else if(scale == "gradient2"){
      ggpl <- ggpl + scale_fill_gradient2(name = legend_colors_point_title, low = color_low_point, mid = color_mid_point, high = color_high_point, midpoint = midpoint)
    }else if(scale == "gradient"){
      ggpl <- ggpl + scale_fill_gradient(name = legend_colors_point_title, low = color_low_point, high = color_high_point)
    }
    
  }else{
    
    ggpl <- ggpl +
      geom_point(aes(color = .data[[color_point_var]]), size = point_size, shape = point_shape, alpha = point_alpha, stroke = point_stroke,
        show.legend = legend_show_colors_point)
    
    if(scale == "manual"){
      ggpl <- ggpl + scale_color_manual(name = legend_colors_point_title, values = colors_point, drop = FALSE)
    }else if(scale == "gradientn") {
      ggpl <- ggpl + scale_color_gradientn(name = legend_colors_point_title, colors = colors_point)
    }else if(scale == "gradient2"){
      ggpl <- ggpl + scale_color_gradient2(name = legend_colors_point_title, low = color_low_point, mid = color_mid_point, high = color_high_point, midpoint = midpoint)
    }else if(scale == "gradient"){
      ggpl <- ggpl + scale_color_gradient(name = legend_colors_point_title, low = color_low_point, high = color_high_point)
    }
    
  }
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle) + 
    ylab(ylab) +
    xlab(xlab) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = title_size),
      axis.line = element_blank(),
      axis.ticks = element_line(color = "black", size = 0.5),
      panel.border = element_rect(colour = "black", size = 0.8),
      legend.position = legend_position) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.15) +
    coord_cartesian(xlim = xlim, ylim = ylim)
  
  
  if(smooth == "pooled"){
    
    ggpl <- ggpl + 
      geom_smooth(aes(group = 1), 
        method = smooth_method, formula = smooth_formula, se = smooth_se, 
        linetype = smooth_linetype, size = smooth_size)
    
  }else if(smooth == "strat"){
    
    ggpl <- ggpl + 
      geom_smooth(aes(group = .data[[color_point_var]], color = .data[[color_point_var]]), 
        method = smooth_method, formula = smooth_formula, se = smooth_se, 
        linetype = smooth_linetype, size = smooth_size)
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




#' @rdname wrapper_point_plot_core
#' @param strat1_var Name of the first stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' @export
wrapper_point_plot_core_strat <- function(data, x_var, y_var, color_point_var = NULL, facet_var = NULL, 
  strat1_var = NULL, strat2_var = NULL, 
  colors_point = NULL, scale_gradient = "gradientn", color_low_point = '#42399B', color_mid_point = "white", color_high_point = '#D70131', midpoint = 0,
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_point_title = NULL, legend_position = "right", facet_label_both = TRUE, 
  point_size = 1.5, point_shape = 20, point_alpha = 1, point_stroke = 0.8,
  smooth = "none", smooth_method = "auto", smooth_formula = y ~ x, smooth_se = FALSE,
  smooth_size = 2, smooth_linetype = 1, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", xlim = NULL, ylim = NULL, 
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
  
  data <- data[stats::complete.cases(data[, c(x_var, y_var, facet_var, strat1_var, strat2_var)]), , drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  # -------------------------------------------------------------------------
  # Scales, xlim, ylim
  # -------------------------------------------------------------------------
  
  if(strat_scales == "fixed" || strat_scales == "free_y"){
    if(is.null(xlim)){
      xlim <- range(data[, x_var])
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
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
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
      
      
      ggpl <- wrapper_point_plot_core(data = data_strata1, x_var = x_var, y_var = y_var, color_point_var = color_point_var, facet_var = facet_var, 
        colors_point = colors_point, scale_gradient = scale_gradient, color_low_point = color_low_point, color_mid_point = color_mid_point, color_high_point = color_high_point, midpoint = midpoint,
        variable_names = variable_names, 
        xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, 
        legend_colors_point_title = legend_colors_point_title, legend_position = legend_position, facet_label_both = facet_label_both, 
        point_size = point_size, point_shape = point_shape, point_alpha = point_alpha, point_stroke = point_stroke,
        smooth = smooth, smooth_method = smooth_method, smooth_formula = smooth_formula, smooth_se = smooth_se,
        smooth_size = smooth_size, smooth_linetype = smooth_linetype,
        title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, xlim = xlim, ylim = ylim, 
        background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    if(less_legends && legend_position == "right"){
      
      # Extract the legend from one of the plots
      legend <- get_legend(
        # Create some space to the left of the legend
        ggpl[[1]] + theme(legend.box.margin = margin(0, 0, 0, 12))
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

















