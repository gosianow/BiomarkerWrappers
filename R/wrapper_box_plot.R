



# data <- data_goya
# data$Cell_Of_Origin[data$Cell_Of_Origin == "ABC"] <- NA
# 
# x_var <- "Cell_Of_Origin"
# y_var <- "FCGR3A"
# colors = NULL
# variable_names = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# show_total_counts = TRUE
# point_size = 1.5
# title.size = 12
# ylim = NULL
# axis.text.x.angle = 0
# axis.text.x.vjust = 0
# axis.text.x.hjust = 0.5
# background_grid_major = "none"



#' Boxplot
#' 
#' Generate a signle boxplot.
#' 
#' @param data Data frame.
wrapper_core_box_plot <- function(data, x_var, y_var, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, show_total_counts = TRUE, point_size = 1.5, title.size = 12, ylim = NULL, axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, background_grid_major = "none"){
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  
  colors <- format_colors(levels = levels(data[, x_var]), colors = colors)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  ### Keep non-missing data
  
  data <- data[!is.na(data[, x_var]) & !is.na(data[, y_var]), ]
  
  
  ### Calculate counts per subgroup
  if(show_total_counts){
    
    tbl <- table(data[, x_var])
    
    data[, x_var] <- factor(data[, x_var], levels = levels(data[, x_var]), labels = paste0(levels(data[, x_var]), "\n(", tbl, ")"))
    
    ## Update names for colors
    
    names(colors) <- levels(data[, x_var])
    
  }
  
  
  xlab <- variable_names[x_var]
  ylab <- variable_names[y_var]
  
  
  ### Make the plot
  
  
  ggpl <- ggplot(data, aes_string(x = x_var, y = y_var)) +
    geom_boxplot(aes_string(fill = x_var), outlier.color = NA) +
    geom_quasirandom(width = 0.25, size = point_size, shape = 1) +
    labs(title = title, subtitle = subtitle, tag = tag) + 
    ylab(ylab) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = title.size, face = "bold"),
      plot.subtitle = element_text(size = title.size),
      axis.text.x = element_text(angle = axis.text.x.angle, vjust = axis.text.x.vjust, hjust = axis.text.x.hjust),
      legend.position = "none",
      plot.tag.position = "top",
      plot.tag = element_text(size = title.size, face = "plain")) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
    scale_fill_manual(values = colors, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = ylim)
  
  
  
  return(ggpl)
  
  
}



# data <- data_goya
# 
# x_var <- "BCL2_cat2"
# y_var <- "FCGR3A"
# 
# 
# strat1_var = "Treatment_Arm"
# strat2_var = "Cell_Of_Origin2"
# strat2_var = NULL
# 
# colors = NULL
# variable_names = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# show_total_counts = TRUE
# point_size = 1.5
# title.size = 12
# ylim = NULL
# axis.text.x.angle = 0
# axis.text.x.vjust = 0
# axis.text.x.hjust = 0.5
# background_grid_major = "none"
# 
# 
# strat1_nrow = 1
# strat1_ncol = NULL
# strat2_nrow = NULL
# strat2_ncol = 1


#' Boxplot
#' 
#' Generate box plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_box_plot_strat <- function(data, x_var, y_var, strat1_var = NULL, strat2_var = NULL, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, show_total_counts = TRUE, point_size = 1.5, title.size = 12, ylim = NULL, axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, background_grid_major = "none", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
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
  
  data <- data[!is.na(data[, x_var]) & !is.na(data[, y_var]) & !is.na(data[, strat1_var]) & !is.na(data[, strat2_var]), ]
  
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
      
      
      ggpl <- wrapper_core_box_plot(data = data_strata1, x_var = x_var, y_var = y_var, colors = colors, variable_names = variable_names, title = title, subtitle = subtitle, tag = tag, show_total_counts = show_total_counts, point_size = point_size, title.size = title.size, ylim = ylim, axis.text.x.angle = axis.text.x.angle, axis.text.x.vjust = axis.text.x.vjust, axis.text.x.hjust = axis.text.x.hjust, background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  return(ggpl)
  
  
}

















