




#' Error line plot
#' 
#' @param data Data frame.
#' @export
wrapper_error_line_plot_core <- function(data = NULL, x_var, y_var, color_line_var, 
  colors_line = NULL, 
  variable_names = NULL, 
  title = TRUE, subtitle = TRUE, xlab = TRUE, ylab = TRUE,
  ggdata = NULL){
  
  
  # -------------------------------------------------------------------------
  # Axis and legend labels
  # -------------------------------------------------------------------------
  
  
  if(is.logical(title)){
    title <- NULL
  }
  
  if(is.logical(subtitle)){
    subtitle <- NULL
  }
  
  
  
  if(is.null(ggdata)){
    ggdata <- data %>% 
      group_by_at(c(x_var, color_line_var)) %>%
      summarise_at(y_var, list(value = ~ median(., na.rm = TRUE), first_quartile = ~ quantile(., probs = 0.25, na.rm = TRUE), third_quartile = ~ quantile(., probs = 0.75, na.rm = TRUE))) %>%
      data.frame()
  }
  
  # list(value = ~ mean(., na.rm = TRUE), error = ~ sd(., na.rm = TRUE))
  
  ggplot(ggdata, aes(x = .data[[x_var]], y = .data[["value"]], color = .data[[color_line_var]])) +
    geom_point(position = position_dodge(0.5)) +
    geom_errorbar(aes(ymin = first_quartile, ymax = third_quartile), position = position_dodge(0.5), width = 0.3) +
    geom_line() +
    ylab(y_var) +
    background_grid(major = "xy", minor = "none", size.major = 0.15) +
    labs(title = title, subtitle = subtitle)
  
  
}








#' @rdname wrapper_error_line_plot_core
#' @param strat1_var Name of the first stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' @export
wrapper_error_line_plot_core_strat <- function(data, x_var, y_var, color_line_var, 
  strat1_var = NULL, strat2_var = NULL, 
  colors_line = NULL, 
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, ylab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_line_title = TRUE, legend_position = "right",
  xlim = NULL, ylim = NULL, 
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
  
  data <- data[stats::complete.cases(data[, c(x_var, y_var, color_line_var, strat1_var, strat2_var)]), , drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Compute summary data for ggplot
  # -------------------------------------------------------------------------
  
  ggdata <- data %>% 
    group_by_at(c(x_var, color_line_var, strat1_var, strat2_var)) %>%
    summarise_at(y_var, list(value = ~ median(., na.rm = TRUE), first_quartile = ~ quantile(., probs = 0.25, na.rm = TRUE), third_quartile = ~ quantile(., probs = 0.75, na.rm = TRUE))) %>%
    data.frame()
  
  
  data <- ggdata
  
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
      ylim <- range(data[, c("value", "first_quartile", "third_quartile")])
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
      
      ggpl <- wrapper_error_line_plot_core(data = NULL, x_var = x_var, y_var = y_var, color_line_var = color_line_var, ggdata = data_strata1,
        title = title, subtitle = subtitle)
      
      
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















