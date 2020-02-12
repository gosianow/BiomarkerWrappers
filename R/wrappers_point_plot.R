




#' Scatterplot
#' 
#' Generate a signle scatter.
#' 
#' @param data Data frame.
wrapper_core_point_plot <- function(data, x_var, y_var, color_var = NULL, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, point_size = 1.5, title.size = 12, xlim = NULL, ylim = NULL, background_grid_major = "none"){
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.numeric(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  stopifnot(length(color_var) == 1)
  stopifnot(is.factor(data[, color_var]))
  
  
  if(!is.null(color_var)){
    colors <- format_colors(levels = levels(data[, color_var]), colors = colors)
  }else{
    if(is.null(colors)){
      colors <- "lightgrey"
    }
  }
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  ### Keep non-missing data
  
  data <- data[complete.cases(data[, c(x_var, y_var)]), , drop = FALSE]
  
  xlab <- variable_names[x_var]
  ylab <- variable_names[y_var]
  scale_fill_name <- variable_names[color_var]
  
  
  ### Make the plot
  
  
  if(!is.null(color_var)){
    
    
    ggpl <- ggplot(data, aes_string(x = x_var, y = y_var, fill = color_var)) +
      geom_point(size = point_size, shape = 21) +
      labs(title = title, subtitle = subtitle, tag = tag) + 
      ylab(ylab) +
      xlab(xlab) +
      theme_cowplot(12) +
      theme(plot.title = element_text(size = title.size, face = "bold"),
        plot.subtitle = element_text(size = title.size),
        plot.tag.position = "top",
        plot.tag = element_text(size = title.size, face = "plain")) +
      background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
      scale_fill_manual(scale_fill_name, values = colors, drop = FALSE) +
      coord_cartesian(xlim = xlim, ylim = ylim)
    
    
  }else{
    
    ggpl <- ggplot(data, aes_string(x = x_var, y = y_var)) +
      geom_point(size = point_size, shape = 21, color = colors) +
      labs(title = title, subtitle = subtitle, tag = tag) + 
      ylab(ylab) +
      xlab(xlab) +
      theme_cowplot(12) +
      theme(plot.title = element_text(size = title.size, face = "bold"),
        plot.subtitle = element_text(size = title.size),
        plot.tag.position = "top",
        plot.tag = element_text(size = title.size, face = "plain")) +
      background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
      coord_cartesian(xlim = xlim, ylim = ylim)
    
    
  }
  
  
  
  
  
  return(ggpl)
  
  
}




#' Scatterplot
#' 
#' Generate scatter plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_point_plot_strat <- function(data, x_var, y_var, color_var = NULL, strat1_var = NULL, strat2_var = NULL, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, point_size = 1.5, title.size = 12, xlim = NULL, ylim = NULL, background_grid_major = "none", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  if(!is.null(strat1_var)){
    stopifnot(length(strat1_var) == 1)
    stopifnot(is.factor(data[, strat1_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat1_dummy"] <- factor("strat1_dummy")
    strat1_var <- "strat1_dummy"
  }
  
  if(!is.null(strat2_var)){
    stopifnot(length(strat2_var) == 1)
    stopifnot(is.factor(data[, strat2_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat2_dummy"] <- factor("strat2_dummy")
    strat2_var <- "strat2_dummy"
  }
  
  ### Keep non-missing data
  
  data <- data[complete.cases(data[, c(x_var, y_var, strat1_var, strat2_var)]), , drop = FALSE]
  
  
  if(is.null(xlim)){
    xlim <- range(data[, x_var])
  }
  if(is.null(ylim)){
    ylim <- range(data[, y_var])
  }
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
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
      tag <- strata2_levels[j]
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
        subtitle <- strata1_levels[i]
      }
      
      
      ggpl <- wrapper_core_point_plot(data = data_strata1, x_var = x_var, y_var = y_var, color_var = color_var, colors = colors, variable_names = variable_names, title = title, subtitle = subtitle, tag = tag, point_size = point_size, title.size = title.size, xlim = xlim, ylim = ylim, background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  return(ggpl)
  
  
}

















