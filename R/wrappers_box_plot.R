



#' Boxplot
#' 
#' Generate a signle boxplot.
#' 
#' @param data Data frame.
wrapper_core_box_plot <- function(data, x_var, y_var, color_point_var = NULL, dodge_var = NULL, facet_var = NULL, 
  colors_box = NULL, colors_point = NULL, 
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, tag = NULL, 
  legend_title_colors_box = NULL, legend_title_colors_point = NULL, legend_position = "right", facet_label_both = TRUE, 
  show_total_counts = TRUE, show_median = TRUE, 
  point_size = 1, point_shape = 1, point_alpha = 1, 
  label_size = 3.5, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
  background_grid_major = "none"){
  
  
  # -------------------------------------------------------------------------
  # Checks about data
  # -------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  if(!is.null(color_point_var)){
    stopifnot(length(color_point_var) == 1)
    stopifnot(is.factor(data[, color_point_var]))
  }
  
  if(!is.null(dodge_var)){
    stopifnot(length(dodge_var) == 1)
    stopifnot(is.factor(data[, dodge_var]))
  }
  
  if(!is.null(facet_var)){
    stopifnot(length(facet_var) == 1)
    stopifnot(is.factor(data[, facet_var]))
  }
  
  
  ### Keep non-missing data
  
  ## We do not exclude the missing values in color_line_var
  data <- data[complete.cases(data[, c(x_var, y_var, dodge_var, facet_var)]), , drop = FALSE]
  
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
      colors_point <- "black"
    }else{
      colors_point <- colors_point[1]
    }
  }else{
    colors_point <- format_colors(levels = levels(data[, color_point_var]), colors = colors_point)
  }
  
  
  
  if(is.null(dodge_var)){
    if(length(colors_box) == 1){
      colors_box <- rep(colors_box, nlevels(data[, x_var]))
      names(colors_box) <- levels(data[, x_var])
    }else{
      colors_box <- format_colors(levels = levels(data[, x_var]), colors = colors_box)
    }
  }else{
    colors_box <- format_colors(levels = levels(data[, dodge_var]), colors = colors_box)
  }
  
  
  
  
  # -------------------------------------------------------------------------
  # Axis and legend labels
  # -------------------------------------------------------------------------
  
  
  if(is.null(xlab)){
    xlab <- variable_names[x_var]
  }
  if(is.null(ylab)){
    ylab <- variable_names[y_var]
  }
  if(is.null(legend_title_colors_point)){
    legend_title_colors_point <- variable_names[color_point_var]
  }
  if(is.null(legend_title_colors_box)){
    legend_title_colors_box <- variable_names[dodge_var]
  }
  
  
  # --------------------------------------------------------------------------
  ### Calculate counts per subgroup
  ### Calculate median per subgroup
  # --------------------------------------------------------------------------
  
  subgroup_vars <- c(x_var, facet_var, dodge_var)
  stopifnot(sum(duplicated(subgroup_vars)) == 0)
  
  
  N <- aggregate(data[, y_var], lapply(subgroup_vars, function(x) data[, x]), FUN = length, drop = FALSE)
  colnames(N) <- c(subgroup_vars, "N")
  
  
  Median <- aggregate(data[, y_var], lapply(subgroup_vars, function(x) data[, x]), FUN = median, na.rm = TRUE, drop = FALSE)
  colnames(Median) <- c(subgroup_vars, "Median")
  
  
  ggdata_summ <- N %>% 
    left_join(Median, by = subgroup_vars)
  
  
  ggdata_summ[, x_var] <- factor(ggdata_summ[, x_var], levels = levels(data[, x_var]))
  
  
  if(!is.null(dodge_var)){
    ggdata_summ[, dodge_var] <- factor(ggdata_summ[, dodge_var], levels = levels(data[, dodge_var]))
  }
  if(!is.null(facet_var)){
    ggdata_summ[, facet_var] <- factor(ggdata_summ[, facet_var], levels = levels(data[, facet_var]))
  }
  
  
  ### Prepare labels
  if(show_total_counts){
    if(show_median){
      ggdata_summ$Label <- paste0("(", paste(ggdata_summ$N, round(ggdata_summ$Median, 2), sep = ", "), ")")
    }else{
      ggdata_summ$Label <- paste0("(", paste(ggdata_summ$N, sep = ", "), ")")
    }
    
  }else{
    if(show_median){
      ggdata_summ$Label <- paste0("(", paste(round(ggdata_summ$Median, 2), sep = ", "), ")")
    }else{
      ggdata_summ$Label <- NULL
    }
  }
  
  
  ### Remove NA cases so they are not displayed as label
  ggdata_summ <- ggdata_summ[!is.na(ggdata_summ$N), ]
  
  
  
  # --------------------------------------------------------------------------
  # Make the plot
  # --------------------------------------------------------------------------
  
  ### We use NA instead of TRUE, because if TRUE legend for x-axis is displayed.
  legend_show_colors_point <- ifelse(color_point_var == "color_point_dummy", FALSE, NA)
  
  
  if(is.null(dodge_var)){
    
    ggpl <- ggplot(data, aes_string(x = x_var, y = y_var)) +
      geom_boxplot(aes_string(fill = x_var), outlier.color = NA, show.legend = FALSE) +
      scale_fill_manual(name = legend_title_colors_box, values = colors_box, drop = FALSE)
    
    
    if(point_shape %in% 21:25){
      
      ggpl <- ggpl +
        new_scale_fill() +
        geom_jitter(aes_string(fill = color_point_var), width = 0.3, size = point_size, shape = point_shape, alpha = point_alpha, show.legend = legend_show_colors_point) +
        scale_fill_manual(name = legend_title_colors_point, values = colors_point, drop = FALSE, na.value = "grey") 
      
      
    }else{
      
      ggpl <- ggpl +
        geom_jitter(aes_string(color = color_point_var), width = 0.3, size = point_size, shape = point_shape, alpha = point_alpha, show.legend = legend_show_colors_point) +
        scale_color_manual(name = legend_title_colors_point, values = colors_point, drop = FALSE, na.value = "grey") 
      
    }
    
    
  }else{
    
    
    ggpl <- ggplot(data, aes_string(x = x_var, y = y_var, fill = dodge_var)) +
      geom_boxplot(outlier.color = NA, position = position_dodge2(preserve = "single", width = 0.75)) +
      scale_fill_manual(name = legend_title_colors_box, values = colors_box, drop = FALSE)
    
    
    if(point_shape %in% 21:25){
      
      ggpl <- ggpl +
        new_scale_fill() +
        geom_jitter(aes_string(fill = color_point_var, group = dodge_var), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = point_size, shape = point_shape, alpha = point_alpha, show.legend = legend_show_colors_point) +
        scale_fill_manual(name = legend_title_colors_point, values = colors_point, drop = FALSE, na.value = "grey") 
      
      
    }else{
      
      ## The group determines dodging 
      ggpl <- ggpl +
        geom_jitter(aes_string(color = color_point_var, group = dodge_var), position = position_jitterdodge(jitter.width = 0.3, dodge.width = 0.75), size = point_size, shape = point_shape, alpha = point_alpha, show.legend = legend_show_colors_point) +
        scale_color_manual(name = legend_title_colors_point, values = colors_point, drop = FALSE, na.value = "grey") 
      
    }
    
    
  }
  
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle, tag = tag) + 
    ylab(ylab) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = title_size),
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      plot.tag.position = "top",
      plot.tag = element_text(size = title_size, face = "plain"),
      legend.position = legend_position) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = ylim)
  
  
  
  
  ### Facet
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
      facet_wrap(as.formula(paste("~", facet_var)), labeller = labeller, scales = facet_scales) +
      theme(strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = strip_text_size),
        axis.line = element_blank()) +
      panel_border(colour = "black", linetype = 1, size = 1, remove = FALSE)
    
  }
  
  
  
  
  ### Labels
  if(show_total_counts || show_median){
    ## ggplot2 doesn't know you want to give the labels the same virtual width as the bars. So tell it. 
    ## You can't nudge and dodge text, so instead adjust the y position.
    
    
    ### Calculate the place at the bottom where labels should be placed
    if(!is.null(ylim)){
      ymin <- ylim[1]
    }else{
      ymin <- min(data[, y_var], na.rm = TRUE)  
    }
    
    ### The default expand is 5% of the range, so let's nudge by 2.5%
    ### Calculate the nudge
    if(!is.null(ylim)){
      yrange <- range(ylim)
    }else{
      yrange <- range(data[, y_var], na.rm = TRUE)  
    }
    yrange <- yrange[2] - yrange[1]
    ynudge <- yrange * 0.025
    
    
    if(is.null(dodge_var)){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_summ, aes_string(x = x_var, y = ymin - ynudge, label = "Label"), size = label_size, vjust = 0.5)
      
    }else{
      
      ggpl <- ggpl +
        geom_text(data = ggdata_summ, aes_string(x = x_var, y = ymin - ynudge, group = dodge_var, label = "Label"), size = label_size, vjust = 0.5, position = position_dodge(preserve = "total", width = 0.9))
      
      
    }
    
    
  }
  
  
  ggpl
  
  
  
}






#' Boxplot
#' 
#' Generate box plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_box_plot_strat <- function(data, x_var, y_var, color_point_var = NULL, dodge_var = NULL, facet_var = NULL,
  strat1_var = NULL, strat2_var = NULL, 
  colors_box = NULL, colors_point = NULL,
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors_box = NULL, legend_title_colors_point = NULL, legend_position = "right", facet_label_both = TRUE, 
  show_total_counts = TRUE, show_median = TRUE, 
  point_size = 1, point_shape = 1, point_alpha = 1, 
  label_size = 3.5,
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
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
  
  data <- data[complete.cases(data[, c(x_var, y_var, dodge_var, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Scales, xlim, ylim
  # -------------------------------------------------------------------------
  
  if(strat_scales == "fixed"){
    if(is.null(ylim)){
      ylim <- range(data[, y_var])
    }
  }
  
  # -------------------------------------------------------------------------
  # Lapply to make the stratified plots
  # -------------------------------------------------------------------------
  
  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  legend_positions <- rep(legend_position, length(strata1_levels))
  if(less_legends){
    if(!is.null(strat1_nrow)){
      strat1_ncol <- ceiling(length(strata1_levels) / strat1_nrow)
    }
    indx_keep <- seq(strat1_ncol, length(strata1_levels), by = strat1_ncol)
    legend_positions[-indx_keep] <- "none"
  }

  
  
  ggpl <- lapply(1:length(strata2_levels), function(j){
    # j = 1
    
    data_strata2 <- data[data[, strat2_var] == strata2_levels[j] & !is.na(data[, strat2_var]), ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    ### Tag
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
      
      ### Subtitle
      if(strat1_var == "strat1_dummy"){
        subtitle <- NULL
      }else{
        if(subtitle_label_both){
          subtitle <- paste0(variable_names[strat1_var], ": ", strata1_levels[i])
        }else{
          subtitle <- strata1_levels[i]
        }
      }
      
      legend_position <- legend_positions[i]
      
      ggpl <- wrapper_core_box_plot(data = data_strata1, x_var = x_var, y_var = y_var, color_point_var = color_point_var, dodge_var = dodge_var, facet_var = facet_var, 
        colors_box = colors_box, colors_point = colors_point, 
        variable_names = variable_names, 
        xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, tag = tag, 
        legend_title_colors_box = legend_title_colors_box, legend_title_colors_point = legend_title_colors_point, legend_position = legend_position, facet_label_both = facet_label_both, 
        show_total_counts = show_total_counts, show_median = show_median, 
        point_size = point_size, point_shape = point_shape, point_alpha = point_alpha, 
        title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, ylim = ylim, 
        axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust, 
        label_size = label_size, 
        background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  ggpl
  
  
}








#' Boxplots
#' 
#' Generate boxplots for multiple variables in one faceted or dodged panel.
#' 
#' @param data Data frame.
wrapper_core_box_plot_yvars_strat <- function(data, y_vars, x_var = NULL, color_point_var = NULL, dodge_var = NULL, facet_var = NULL, 
  strat1_var = NULL, strat2_var = NULL, 
  colors_box = NULL, colors_point = NULL, 
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors_box = NULL, legend_title_colors_point = NULL, legend_position = "right", facet_label_both = TRUE, 
  show_total_counts = TRUE, show_median = TRUE, 
  point_size = 1, point_shape = 1, point_alpha = 1, 
  label_size = 3.5, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
  background_grid_major = "none", 
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1, less_legends = FALSE,
  names_to = "name", values_to = "value"){
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "numeric"))
  
  ### We use either x_var or facet_var to present the multiple y_vars. That is why x_var and facet_var cannot be both defined. At least one of them must be NULL.
  stopifnot(is.null(x_var) || is.null(facet_var))
  
  
  data <- data[, c(y_vars, x_var, color_point_var, dodge_var, facet_var), drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # pivot_longer the data from y_vars
  # --------------------------------------------------------------------------
  
  data_longer <- tidyr::pivot_longer(data, y_vars, names_to = names_to, values_to = values_to) %>% 
    as.data.frame()
  
  data_longer[, names_to] <- factor(data_longer[, names_to], levels = y_vars, labels = variable_names[y_vars])
  
  y_var <- values_to
  
  if(is.null(x_var)){
    x_var <- names_to
  }else if(is.null(facet_var)){
    facet_var <- names_to
    facet_label_both <- FALSE
  }
  
  
  ggpl <- wrapper_core_box_plot_strat(data = data_longer, x_var = x_var, y_var = y_var, color_point_var = color_point_var, dodge_var = dodge_var, facet_var = facet_var, 
    colors_box = colors_box, colors_point = colors_point, 
    strat1_var = strat1_var, strat2_var = strat2_var, 
    variable_names = variable_names, 
    xlab = xlab, ylab = ylab, title = title, subtitle_label_both = subtitle_label_both, tag_label_both = tag_label_both,
    legend_title_colors_box = legend_title_colors_box, legend_title_colors_point = legend_title_colors_point, legend_position = legend_position, facet_label_both = facet_label_both, 
    show_total_counts = show_total_counts, show_median = show_median, 
    point_size = point_size, point_shape = point_shape, point_alpha = point_alpha, 
    title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, ylim = ylim, 
    axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust, 
    label_size = label_size, 
    background_grid_major = background_grid_major, 
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol, less_legends = less_legends)
  
  
  
  ggpl
  
  
}






































