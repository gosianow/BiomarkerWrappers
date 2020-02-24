







# ggplot(mtcars, aes(factor(cyl), fill = factor(vs))) +
#   geom_bar(position = position_dodge2(preserve = "single"))
# 
# 
# ggplot(mtcars, aes(factor(cyl), fill = factor(vs))) +
#   geom_bar(position = position_dodge(preserve = "single", width = 1), width = 0.9)



# ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(vs))) +
#   geom_boxplot(outlier.color = NA, position = position_dodge2(preserve = "single"))
# 
# 
# ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(vs))) +
#   geom_boxplot(outlier.color = NA, position = position_dodge(preserve = "single", width = 0.9))
# 
# 
# 
# 
# ### Points and boxplots are not aligned!
# ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(vs))) +
#   geom_boxplot(outlier.color = NA, position = position_dodge(preserve = "single", width = 0.9)) +
#   geom_jitter(shape = 21, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.9))
# 
# 
# ### This works
# ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(vs))) +
#   geom_boxplot(outlier.color = NA, position = position_dodge2(preserve = "single")) +
#   geom_jitter(shape = 21, position = position_jitterdodge(jitter.width = 0.25))
# 
# 
# ggplot(mtcars, aes(x = factor(cyl), y = mpg, fill = factor(vs))) +
#   geom_boxplot(outlier.color = NA, position = position_dodge(preserve = "total"))
# 
# 
# 
# # https://github.com/tidyverse/ggplot2/issues/2712
# # I didn't expect that one needs to pass different positions to the two geoms for them to line up, but I guess it makes sense to always use dodge2 for geoms that have width and dodge for geoms that don't.
# 
# ggplot(mtcars, aes(factor(gear), mpg, fill = factor(am))) +
#   geom_boxplot(position = position_dodge2(0.75, preserve = 'single')) +
#   geom_point(position = position_dodge(0.75, preserve = 'total'))
# 
# 
# ### The points are not centered when the width is different than 0.75 
# ggplot(mtcars, aes(factor(gear), mpg, fill = factor(am))) +
#   geom_boxplot(position = position_dodge2(0.9, preserve = 'single')) +
#   geom_point(position = position_dodge(0.9, preserve = 'total'))
# 
# 
# ggplot(mtcars, aes(factor(gear), mpg, fill = factor(am))) +
#   geom_boxplot(position = position_dodge2(0.5, preserve = 'single')) +
#   geom_point(position = position_dodge(0.5, preserve = 'total'))


# Additionally, ggplot2 doesn't know you want to give the labels the same virtual width as the bars. So tell it. You can't nudge and dodge text, so instead adjust the y position.

### Based on that we use:
# geom_boxplot with position_dodge2(preserve = "single", width = 0.75) 
# geom_jitter with position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75)
# geom_text with position_dodge(preserve = "total", width = 0.75)




# colors_box = NULL
# colors_point = NULL
# variable_names = NULL
# xlab = NULL
# ylab = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# show_total_counts = TRUE
# show_median = TRUE
# point_size = 1
# point_shape = 1
# title.size = 12
# ylim = NULL
# axis.text.x.angle = 0
# axis.text.x.vjust = 0
# axis.text.x.hjust = 0.5
# geom_text_size = 3.5
# background_grid_major = "none"
# strip.text.size = NULL


#' Boxplot
#' 
#' Generate a signle boxplot.
#' 
#' @param data Data frame.
wrapper_core_box_plot <- function(data, x_var, y_var, facet_var = NULL, fill_var = NULL, color_var = NULL, colors_box = NULL, colors_point = NULL, variable_names = NULL, xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, tag = NULL, show_total_counts = TRUE, show_median = TRUE, point_size = 1, point_shape = 1, title.size = 12, ylim = NULL, axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, geom_text_size = 3.5, background_grid_major = "none", strip.text.size = NULL){
  
  
  stopifnot(is.data.frame(data))
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  if(!is.null(facet_var)){
    stopifnot(length(facet_var) == 1)
    stopifnot(is.factor(data[, facet_var]))
  }
  
  if(!is.null(color_var)){
    stopifnot(length(color_var) == 1)
    stopifnot(is.factor(data[, color_var]))
  }
  
  if(!is.null(fill_var)){
    stopifnot(length(fill_var) == 1)
    stopifnot(is.factor(data[, fill_var]))
  }
  
  
  if(!is.null(fill_var)){
    colors_box <- format_colors(levels = levels(data[, fill_var]), colors = colors_box)
    legend_name_fill <- variable_names[fill_var]
  }else{
    if(length(colors_box) == 1){
      colors_box <- rep(colors_box, nlevels(data[, x_var]))
      names(colors_box) <- levels(data[, x_var])
    }else{
      colors_box <- format_colors(levels = levels(data[, x_var]), colors = colors_box)
    }
    legend_name_fill <- NULL
  }
  
  
  if(!is.null(color_var)){
    colors_point <- format_colors(levels = levels(data[, color_var]), colors = colors_point)
    legend_name_color <- variable_names[color_var]
  }
  
  
  ### Keep non-missing data
  data <- data[complete.cases(data[, c(x_var, y_var, facet_var, fill_var)]), ]
  
  
  # data_expand <- expand.grid(lapply(data[, c(x_var, fill_var)], levels))
  # data_expand[, y_var] <- NA
  # data <- rbind.fill(data, data_expand)
  
  
  # --------------------------------------------------------------------------
  ### Calculate counts per subgroup
  ### Calculate median per subgroup
  # --------------------------------------------------------------------------
  
  
  N <- aggregate(data[, y_var], lapply(c(x_var, facet_var, fill_var), function(x) data[, x]), FUN = length, drop = FALSE)
  colnames(N) <- c(x_var, facet_var, fill_var, "N")
  
  
  Median <- aggregate(data[, y_var], lapply(c(x_var, facet_var, fill_var), function(x) data[, x]), FUN = median, na.rm = TRUE, drop = FALSE)
  colnames(Median) <- c(x_var, facet_var, fill_var, "Median")
  
  
  ggdata_summ <- N %>% 
    left_join(Median, by = c(x_var, facet_var, fill_var))
  
  
  ggdata_summ[, x_var] <- factor(ggdata_summ[, x_var], levels = levels(data[, x_var]))
  
  
  if(!is.null(fill_var)){
    ggdata_summ[, fill_var] <- factor(ggdata_summ[, fill_var], levels = levels(data[, fill_var]))
  }
  if(!is.null(facet_var)){
    ggdata_summ[, facet_var] <- factor(ggdata_summ[, facet_var], levels = levels(data[, facet_var]))
  }
  
  
  ### Prepare label
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
  ### Make the plot
  # --------------------------------------------------------------------------
  
  if(is.null(xlab)){
    xlab <- variable_names[x_var]
  }
  if(is.null(ylab)){
    ylab <- variable_names[y_var]
  }
  
  
  if(!is.null(fill_var)){
    
    ggpl <- ggplot(data, aes_string(x = x_var, y = y_var, fill = fill_var)) +
      geom_boxplot(outlier.color = NA, position = position_dodge2(preserve = "single", width = 0.75))
    
    
    if(!is.null(color_var)){
      ## The group determines dodging 
      ggpl <- ggpl +
        geom_jitter(aes_string(color = color_var, group = fill_var), size = point_size, shape = point_shape, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75)) +
        scale_color_manual(name = legend_name_color, values = colors_point, drop = FALSE) 
      
    }else{
      ggpl <- ggpl +
        geom_jitter(size = point_size, shape = point_shape, position = position_jitterdodge(jitter.width = 0.25, dodge.width = 0.75))
      
    }
    
    
  }else{
    
    ggpl <- ggplot(data, aes_string(x = x_var, y = y_var)) +
      geom_boxplot(aes_string(fill = x_var), outlier.color = NA, show.legend = FALSE)
    
    if(!is.null(color_var)){
      ggpl <- ggpl +
        geom_jitter(aes_string(color = color_var), size = point_size, shape = point_shape, width = 0.25) +
        scale_color_manual(name = legend_name_color, values = colors_point, drop = FALSE) 
    }else{
      ggpl <- ggpl +
        geom_jitter(size = point_size, shape = point_shape, width = 0.25)
    }
    
    
  }
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle, tag = tag) + 
    ylab(ylab) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = title.size, face = "bold"),
      plot.subtitle = element_text(size = title.size),
      axis.text.x = element_text(angle = axis.text.x.angle, vjust = axis.text.x.vjust, hjust = axis.text.x.hjust),
      plot.tag.position = "top",
      plot.tag = element_text(size = title.size, face = "plain"),
      strip.background = element_rect(colour = "white", fill = "white"),
      strip.text = element_text(size = strip.text.size)) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
    scale_fill_manual(name = legend_name_fill, values = colors_box, drop = FALSE) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = ylim)
  
  
  if(!is.null(facet_var)){
    
    ggpl <- ggpl +
      facet_wrap(as.formula(paste("~", facet_var))) +
      theme(axis.line = element_blank()) +
      panel_border(colour = "black", linetype = 1, size = 1, remove = FALSE)
    
  }
  
  
  if(show_total_counts || show_median){
    
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
    
    
    if(!is.null(fill_var)){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_summ, aes_string(x = x_var, y = ymin - ynudge, group = fill_var, label = "Label"), size = geom_text_size, vjust = 1, position = position_dodge(preserve = "total", width = 0.9))
      
      
    }else{
      
      ggpl <- ggpl +
        geom_text(data = ggdata_summ, aes_string(x = x_var, y = ymin - ynudge, label = "Label"), size = geom_text_size, vjust = 1)
      
    }
    
    
  }
  
  
  return(ggpl)
  
}





#' Boxplot
#' 
#' Generate box plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_box_plot_strat <- function(data, x_var, y_var, facet_var = NULL, fill_var = NULL, color_var = NULL, colors_box = NULL, colors_point = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, tag = NULL, show_total_counts = TRUE, show_median = TRUE, point_size = 1, point_shape = 1, title.size = 12, ylim = NULL, axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, background_grid_major = "none", strip.text.size = NULL, strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
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
  
  data <- data[complete.cases(data[, c(x_var, y_var, fill_var, strat1_var, strat2_var)]), ]
  
  
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
      
      
      ggpl <- wrapper_core_box_plot(data = data_strata1, x_var = x_var, y_var = y_var, facet_var = facet_var, color_var = color_var, fill_var = fill_var, colors_box = colors_box, colors_point = colors_point, variable_names = variable_names, xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, tag = tag, show_total_counts = show_total_counts, show_median = show_median, point_size = point_size, point_shape = point_shape, title.size = title.size, ylim = ylim, axis.text.x.angle = axis.text.x.angle, axis.text.x.vjust = axis.text.x.vjust, axis.text.x.hjust = axis.text.x.hjust, background_grid_major = background_grid_major, strip.text.size = strip.text.size)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  return(ggpl)
  
  
}






# x_var = NULL
# facet_var = NULL
# fill_var = NULL
# color_var = NULL
# colors_box = NULL
# colors_point = NULL
# strat1_var = NULL
# strat2_var = NULL
# variable_names = NULL
# xlab = NULL
# ylab = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# show_total_counts = TRUE
# show_median = TRUE
# point_size = 1
# point_shape = 1
# title.size = 12
# ylim = NULL
# axis.text.x.angle = 0
# axis.text.x.vjust = 0
# axis.text.x.hjust = 0.5
# background_grid_major = "none"
# strip.text.size = NULL
# strat1_nrow = 1
# strat1_ncol = NULL
# strat2_nrow = NULL
# strat2_ncol = 1





#' Boxplot
#' 
#' Generate a signle boxplot.
#' 
#' @param data Data frame.
wrapper_core_box_plot_yvars_strat <- function(data, y_vars, x_var = NULL, facet_var = NULL, fill_var = NULL, color_var = NULL, colors_box = NULL, colors_point = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, tag = NULL, show_total_counts = TRUE, show_median = TRUE, point_size = 1, point_shape = 1, title.size = 12, ylim = NULL, axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, background_grid_major = "none", strip.text.size = NULL, strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "numeric"))
  
  stopifnot(is.null(x_var) || is.null(facet_var))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  data <- data[, c(y_vars, x_var, facet_var, fill_var, color_var), drop = FALSE]
  
  
  # --------------------------------------------------------------------------
  # pivot_longer the data from y_vars
  # --------------------------------------------------------------------------
  
  data_longer <- tidyr::pivot_longer(data, y_vars) %>% 
    as.data.frame()
  
  data_longer[, "name"] <- factor(data_longer[, "name"], levels = y_vars, labels = variable_names[y_vars])
  
  y_var <- "value"
  
  if(is.null(x_var)){
    x_var <- "name"
  }else if(is.null(facet_var)){
    facet_var <- "name"
  }
  
  
  ggpl <- wrapper_core_box_plot_strat(data_longer, x_var = x_var, y_var = y_var, facet_var = facet_var, fill_var = fill_var, color_var = color_var, colors_box = colors_box, colors_point = colors_point, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, tag = tag, show_total_counts = show_total_counts, show_median = show_median, point_size = point_size, point_shape = point_shape, title.size = title.size, ylim = ylim, axis.text.x.angle = axis.text.x.angle, axis.text.x.vjust = axis.text.x.vjust, axis.text.x.hjust = axis.text.x.hjust, background_grid_major = background_grid_major, strip.text.size = strip.text.size, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  
  return(ggpl)
  
  
}






































