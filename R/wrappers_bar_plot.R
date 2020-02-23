

# data <- data_goya
# x_var <- "Cell_Of_Origin"
# y_var <- "Ann_Arbor_Stage"
# colors = NULL
# variable_names = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# skip_levels = NULL
# show_total_counts = FALSE
# show_proportions = TRUE
# show_total_proportions = TRUE
# title.size = 12
# ylim = c(0, 100)
# axis.text.x.angle = 0
# axis.text.x.vjust = 0
# axis.text.x.hjust = 0.5
# geom_text_size = 3
# geom_text_vjust = 0.5
# background_grid_major = "none"



#' Barplot
#' 
#' Generate a signle barplot.
#' 
#' @param data Data frame.
wrapper_core_bar_plot <- function(data, x_var, y_var, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, skip_levels = NULL, show_total_counts = TRUE, show_proportions = TRUE, show_total_proportions = FALSE, title.size = 12, ylim = c(0, 100), axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, geom_text_size = 3.5, geom_text_vjust = 0.5, background_grid_major = "none"){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.factor(data[, y_var]))
  
  keep_levels <- setdiff(levels(data[, y_var]), skip_levels)
  stopifnot(length(keep_levels) > 0)
  
  ### min must be zero
  if(!is.null(ylim)){
    ylim[1] <- 0
  }
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  ### Calculate counts and proportions
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, y_var], data[, x_var])
  
  prop <- prop.table(tbl, margin = 2) * 100
  
  
  countdf <- as.data.frame.matrix(tbl)
  countdf <- data.frame(Observation = rownames(countdf), countdf, check.names = FALSE, row.names = NULL)
  
  propdf <- as.data.frame.matrix(prop)
  propdf <- data.frame(Observation = rownames(propdf), propdf, check.names = FALSE, row.names = NULL)
  
  
  ### Prepare data for ggplot
  
  ggdata <- cbind(gather(propdf, key = "Subgroup", value = "Proportion", -Observation), 
    gather(countdf, key = "Subgroup", value = "Count", -Observation)[, "Count", drop = FALSE])
  
  
  ### Remove data with zero counts so the ggplot warning with missing values is not displayed
  ggdata <- ggdata[ggdata$Count > 0, ]
  
  
  ### Define lables 
  # ggdata$Label <- paste0(ggdata$Count, " (", round(ggdata$Proportion, 1), "%)")
  ggdata$Label <- paste0(round(ggdata$Proportion, 1), "% (", ggdata$Count, ")")
  
  ggdata$Subgroup <- factor(ggdata$Subgroup, levels = levels(data[, x_var]))
  
  ggdata$Observation <- factor(ggdata$Observation, levels = levels(data[, y_var]))
  
  
  # --------------------------------------------------------------------------
  # Calculate total proportions and counts
  # --------------------------------------------------------------------------
  
  ## Total counts
  
  sum_count_total <- aggregate(ggdata[, "Count"], list(Subgroup = ggdata[, "Subgroup"]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count_total)[2] <- "Count_Total"
  
  ## For zero counts, there is NA. Let's fix it to zero
  sum_count_total[is.na(sum_count_total[, 2]), 2] <- 0
  
  
  ### Skip unwanted levels
  
  if(!is.null(skip_levels)){
    
    ggdata <- ggdata[ggdata$Observation %in% keep_levels, , drop = FALSE]
    
    ggdata$Observation <- factor(ggdata$Observation, keep_levels)
    
  }
  
  
  ## Recalculate total counts for the remaining levels
  
  sum_prop <- aggregate(ggdata[, "Proportion"], list(Subgroup = ggdata[, "Subgroup"]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_prop)[2] <- "Proportion"
  
  sum_count <- aggregate(ggdata[, "Count"], list(Subgroup = ggdata[, "Subgroup"]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count)[2] <- "Count"
  
  
  
  ggdata_total <- sum_count_total %>% 
    left_join(sum_count, by = "Subgroup") %>% 
    left_join(sum_prop, by = "Subgroup")
  
  
  ggdata_total$Subgroup <- factor(ggdata_total$Subgroup, levels = levels(ggdata$Subgroup))
  
  
  ggdata_total$Label <- paste0(round(ggdata_total$Proportion, 1), "% (", ggdata_total$Count, ")")
  
  ggdata_total$Label_Total <- paste0("(", ggdata_total$Count_Total, ")")
  
  
  # --------------------------------------------------------------------------
  ### Make the plot
  # --------------------------------------------------------------------------
  
  
  colors <- format_colors(levels = levels(data[, y_var]), colors = colors)
  
  xlab <- variable_names[x_var]
  ylab <- paste0(variable_names[y_var], "\nProportion (%)")
  
  
  ggpl <- ggplot() +
    geom_col(data = ggdata, aes(x = Subgroup, y = Proportion, fill = Observation), color = "black") +
    labs(title = title, subtitle = subtitle, tag = tag) + 
    ylab(ylab) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = title.size, face = "bold"),
      plot.subtitle = element_text(size = title.size),
      axis.text.x = element_text(angle = axis.text.x.angle, vjust = axis.text.x.vjust, hjust = axis.text.x.hjust),
      legend.title = element_blank(),
      plot.tag.position = "top",
      plot.tag = element_text(size = title.size, face = "plain")) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
    scale_fill_manual(values = colors, drop = FALSE) +
    coord_cartesian(ylim = ylim) 
  # scale_y_continuous(expand = expand_scale(mult = 0.05, add = 0))
  
  
  ### The default expand is 5% of the range, so let's nudge by 2.5%
  ### Calculate the place at the bottom where labels should be placed
  if(!is.null(ylim)){
    ymax <- max(ylim)
  }else{
    ymax <- max(ggdata_total$Proportion, na.rm = TRUE)  
  }
  ynudge <- ymax * 0.025
  
  
  if(show_total_counts){
    
    ggpl <- ggpl +
      geom_text(data = ggdata_total, aes(x = Subgroup, y = 0, label = Label_Total), size = geom_text_size, nudge_y = -ynudge, vjust = 0.5) #  nudge_y = 0, vjust = 1.5
    
  }
  
  
  if(show_proportions){
    
    ggpl <- ggpl +
      geom_text(data = ggdata, aes(x = Subgroup, y = Proportion, group = Observation, label = Label), 
        position = position_stack(vjust = geom_text_vjust), size = geom_text_size, na.rm = TRUE)
    
  }
  
  
  if(show_total_proportions){
    
    ggpl <- ggpl +
      geom_text(data = ggdata_total, aes(x = Subgroup, y = Proportion , label = Label), size = geom_text_size, nudge_y = ynudge, vjust = 0.5) #  nudge_y = 0, vjust = -0.5
    
    
  }
  
  
  
  return(ggpl)
  
  
}









#' Barplot
#' 
#' Generate bar plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_bar_plot_strat <- function(data, x_var, y_var, strat1_var = NULL, strat2_var = NULL, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, skip_levels = NULL, show_total_counts = TRUE, show_proportions = TRUE, show_total_proportions = FALSE, title.size = 12, ylim = c(0, 100), axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, geom_text_size = 3, geom_text_vjust = 0.5, background_grid_major = "none", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
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
      
      ggpl <- wrapper_core_bar_plot(data = data_strata1, x_var = x_var, y_var = y_var, colors = colors, variable_names = variable_names, title = title, subtitle = subtitle, tag = tag, skip_levels = skip_levels, show_total_counts = show_total_counts, show_proportions = show_proportions, show_total_proportions = show_total_proportions, title.size = title.size, ylim = ylim, axis.text.x.angle = axis.text.x.angle, axis.text.x.vjust = axis.text.x.vjust, axis.text.x.hjust = axis.text.x.hjust, geom_text_size = geom_text_size, geom_text_vjust = geom_text_vjust, background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  return(ggpl)
  
  
}






# data <- data_goya
# x_var <- "Cell_Of_Origin"
# y_vars <- c("FCGR1A_cat2", "FCGR2A_cat2", "FCGR3A_cat2")
# y_value <- "Expression"
# colors = NULL
# variable_names = NULL
# title = NULL
# subtitle = NULL
# tag = NULL
# skip_levels = "<=MEDIAN"
# # skip_levels = NULL
# show_proportions = TRUE
# show_total_proportions = TRUE
# title.size = 12
# ylim = c(0, 100)
# axis.text.x.angle = 0
# axis.text.x.vjust = 0
# axis.text.x.hjust = 0.5
# geom_text_size = 3
# geom_text_vjust = 0.5
# background_grid_major = "none"
# 
# method = "facet"
# strip.text.size = NULL




#' Barplot
#' 
#' Generate a signle barplot.
#' 
#' @param data Data frame.
wrapper_core_bar_plot_yvars <- function(data, x_var, y_vars, y_value = NULL, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, skip_levels = NULL, show_total_counts = TRUE, show_proportions = TRUE, show_total_proportions = FALSE, title.size = 12, ylim = c(0, 100), axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, geom_text_size = 3, geom_text_vjust = 0.5, background_grid_major = "none", method = "facet", strip.text.size = NULL){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(length(method) == 1)
  stopifnot(method %in% c("facet", "dodge"))
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(x_var, y_vars), drop = FALSE]
  
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_vars) >= 1)
  
  yvars_class <- sapply(data[, y_vars], class)
  stopifnot(all(yvars_class %in% c("factor")))
  
  ## All y_vars must have the same levels 
  all_yvars_levels <- unique(unlist(lapply(1:length(y_vars), function(i){
    levels(data[, y_vars[i]])
  })))
  
  for(i in 1:length(y_vars)){
    stopifnot(levels(data[, y_vars[i]]) == all_yvars_levels)
  }
  
  
  keep_levels <- setdiff(all_yvars_levels, skip_levels)
  stopifnot(length(keep_levels) > 0)
  
  ## Dodge method cannot be used when there are more than one levels 
  if(length(keep_levels) > 1){
    method <- "facet"
  }
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions
  # --------------------------------------------------------------------------
  
  
  ggdata <- lapply(1:length(y_vars), function(i){
    # i = 1
    
    y_var <- y_vars[i]
    
    
    tbl <- table(data[, y_var], data[, x_var])
    
    prop <- prop.table(tbl, margin = 2) * 100
    
    
    countdf <- as.data.frame.matrix(tbl)
    countdf <- data.frame(Observation = rownames(countdf), countdf, check.names = FALSE, row.names = NULL)
    
    propdf <- as.data.frame.matrix(prop)
    propdf <- data.frame(Observation = rownames(propdf), propdf, check.names = FALSE, row.names = NULL)
    
    
    ### Prepare data for ggplot
    
    ggdata <- cbind(gather(propdf, key = "Subgroup", value = "Proportion", -Observation), 
      gather(countdf, key = "Subgroup", value = "Count", -Observation)[, "Count", drop = FALSE])
    
    
    ### Remove data with zero counts so the ggplot warning with missing values is not displayed
    ggdata <- ggdata[ggdata$Count > 0, ]
    
    
    ### Define lables 
    # ggdata$Label <- paste0(ggdata$Count, " (", round(ggdata$Proportion, 1), "%)")
    ggdata$Label <- paste0(round(ggdata$Proportion, 1), "% (", ggdata$Count, ")")
    
    ggdata$Subgroup <- factor(ggdata$Subgroup, levels = levels(data[, x_var]))
    
    ggdata$Observation <- factor(ggdata$Observation, levels = levels(data[, y_var]))
    
    
    ggdata$y_var <- variable_names[y_var]
    
    
    return(ggdata)
    
    
  })
  
  
  ggdata <- plyr::rbind.fill(ggdata)
  
  ggdata$y_var <- factor(ggdata$y_var, levels = variable_names[y_vars])
  
  
  # --------------------------------------------------------------------------
  # Calculate total proportions and counts
  # --------------------------------------------------------------------------
  
  sum_count_total <- aggregate(ggdata[, "Count"], list(Subgroup = ggdata[, "Subgroup"], y_var = ggdata[, "y_var"]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count_total)[3] <- "Count_Total"
  
  
  ### Skip unwanted levels
  
  if(!is.null(skip_levels)){
    
    ggdata <- ggdata[ggdata$Observation %in% keep_levels, , drop = FALSE]
    
    ggdata$Observation <- factor(ggdata$Observation, keep_levels)
    
  }
  
  
  sum_prop <- aggregate(ggdata[, "Proportion"], list(Subgroup = ggdata[, "Subgroup"], y_var = ggdata[, "y_var"]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_prop)[3] <- "Proportion"
  
  
  sum_count <- aggregate(ggdata[, "Count"], list(Subgroup = ggdata[, "Subgroup"], y_var = ggdata[, "y_var"]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count)[3] <- "Count"
  
  
  ggdata_total <- sum_count_total %>% 
    left_join(sum_count, by = c("Subgroup", "y_var")) %>% 
    left_join(sum_prop, by = c("Subgroup", "y_var"))
  
  
  
  ggdata_total$Subgroup <- factor(ggdata_total$Subgroup, levels = levels(ggdata$Subgroup))
  
  ggdata_total$y_var <- factor(ggdata_total$y_var, levels = levels(ggdata$y_var))
  
  ggdata_total$Label <- paste0(round(ggdata_total$Proportion, 1), "% (", ggdata_total$Count, ")")
  
  ggdata_total$Label_Total <- paste0("(", ggdata_total$Count_Total, ")")
  
  
  
  # --------------------------------------------------------------------------
  ### Make the plot
  # --------------------------------------------------------------------------
  
  if(method == "facet"){
    
    colors <- format_colors(levels = levels(data[, y_vars[1]]), colors = colors)
    
    xlab <- variable_names[x_var]
    ylab <- paste0(ifelse(is.null(y_value), "", paste0(y_value, "\n")), "Proportion (%)")
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = Subgroup, y = Proportion, fill = Observation), color = "black") +
      labs(title = title, subtitle = subtitle, tag = tag) + 
      ylab(ylab) +
      xlab(xlab) +
      theme_cowplot(12) +
      theme(plot.title = element_text(size = title.size, face = "bold"),
        plot.subtitle = element_text(size = title.size),
        axis.text.x = element_text(angle = axis.text.x.angle, vjust = axis.text.x.vjust, hjust = axis.text.x.hjust),
        legend.title = element_blank(),
        plot.tag.position = "top",
        plot.tag = element_text(size = title.size, face = "plain"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = strip.text.size)) +
      background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
      scale_fill_manual(values = colors, drop = FALSE) +
      coord_cartesian(ylim = ylim) +
      facet_wrap(~ y_var)
    
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = Subgroup, y = 0, label = Label_Total), size = geom_text_size, nudge_y = -1, vjust = 1)
      
    }
    
    if(show_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata, aes(x = Subgroup, y = Proportion, group = Observation, label = Label), 
          position = position_stack(vjust = geom_text_vjust), size = geom_text_size, na.rm = TRUE)
      
    }
    
    
    
    if(show_total_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = Subgroup, y = Proportion + 1, label = Label), size = geom_text_size, nudge_y = 1)
      
      
    }
    
    
  }
  
  
  if(method == "dodge"){
    
    colors <- format_colors(levels = levels(data[, x_var]), colors = colors)
    
    xlab <- NULL
    ylab <- paste0(ifelse(is.null(y_value), "", paste0(y_value, " : ")), keep_levels, "\nProportion (%)")
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = y_var, y = Proportion, fill = Subgroup), color = "black", position = "dodge2") +
      labs(title = title, subtitle = subtitle, tag = tag) + 
      ylab(ylab) +
      xlab(xlab) +
      theme_cowplot(12) +
      theme(plot.title = element_text(size = title.size, face = "bold"),
        plot.subtitle = element_text(size = title.size),
        axis.text.x = element_text(angle = axis.text.x.angle, vjust = axis.text.x.vjust, hjust = axis.text.x.hjust),
        legend.title = element_blank(),
        plot.tag.position = "top",
        plot.tag = element_text(size = title.size, face = "plain"),
        strip.background = element_rect(colour = "white", fill = "white"),
        strip.text = element_text(size = strip.text.size)) +
      background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
      scale_fill_manual(values = colors, drop = FALSE) +
      coord_cartesian(ylim = ylim)
    
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = y_var, y = -1, group = Subgroup, label = Label_Total), size = geom_text_size, vjust = 1, position = position_dodge(0.9))
      
    }
    
    ### For one level left in the dodge setup show_proportions is equivalent to show_total_proportions
    
    if(show_proportions || show_total_proportions){
      
      # ggplot2 doesn't know you want to give the labels the same virtual width as the bars.
      # So tell it.
      # You can't nudge and dodge text, so instead adjust the y position.
      
      ggpl <- ggpl +
        geom_text(data = ggdata, aes(x = y_var, y = Proportion + 2, group = Subgroup, label = Label), 
          position = position_dodge(0.9), size = geom_text_size, na.rm = TRUE)
      
    }
    
    
  }
  
  
  return(ggpl)
  
  
}














