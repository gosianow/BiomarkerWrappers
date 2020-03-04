




# facet_var = NULL; 
# colors_bar = NULL; 
# variable_names = NULL; 
# xlab = NULL; ylab = NULL; title = NULL; subtitle = NULL; tag = NULL; 
# legend_title_colors_bar = NULL; facet_label_both = TRUE; 
# skip_levels = NULL; method = "facet"; 
# show_proportions = TRUE; show_counts = TRUE; show_total_proportions = FALSE; show_total_counts = TRUE; 
# label_size = 3.5; label_vjust = 0.5; 
# title_size = 12; strip_text_size = NULL; facet_scales = "fixed"; ylim = NULL; 
# axis_text_x_angle = 0; axis_text_x_vjust = 0; axis_text_x_hjust = 0.5; 
# background_grid_major = "none"




#' Barplot
#' 
#' Generate a signle barplot.
#' 
#' @param data Data frame.
wrapper_core_bar_plot <- function(data, x_var, y_var, facet_var = NULL, 
  colors_bar = NULL, 
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle = NULL, tag = NULL, 
  legend_title_colors_bar = NULL, facet_label_both = TRUE, 
  skip_levels = NULL, method = "facet", 
  show_proportions = TRUE, show_counts = TRUE, show_total_proportions = FALSE, show_total_counts = TRUE, 
  label_size = 3.5, label_vjust = 0.5, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
  background_grid_major = "none"){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.factor(data[, y_var]))
  
  
  if(!is.null(facet_var)){
    stopifnot(length(facet_var) == 1)
    stopifnot(is.factor(data[, facet_var]))
  }else{
    ### Add facet dummy variable to data. 
    ### We do it here because we calculate proportions per facet using lapply. 
    stopifnot(!"facet_dummy" %in% colnames(data))
    data[, "facet_dummy"] <- factor("facet_dummy")
    facet_var <- "facet_dummy"
  }
  
  ### min must be zero
  if(!is.null(ylim)){
    ylim[1] <- 0
  }
  
  keep_levels <- setdiff(levels(data[, y_var]), skip_levels)
  stopifnot(length(keep_levels) > 0)
  
  stopifnot(length(method) == 1)
  stopifnot(method %in% c("facet", "dodge"))
  
  ## Dodge method cannot be used when there are more than one levels 
  if(length(keep_levels) > 1){
    method <- "facet"
  }
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(x_var, y_var, facet_var), drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Colors
  # -------------------------------------------------------------------------
  
  if(method == "facet"){
    colors_bar <- format_colors(levels = levels(data[, y_var]), colors = colors_bar)
  }
  
  if(method == "dodge"){
    colors_bar <- format_colors(levels = levels(data[, x_var]), colors = colors_bar)
  }
  
  # -------------------------------------------------------------------------
  # Axis and legend labels
  # -------------------------------------------------------------------------
  
  if(method == "facet"){
    
    ### By default we do not display legend title as the description of the variable is on the y-axis 
    legend_title_colors_bar <- legend_title_colors_bar
    
    if(is.null(xlab)){
      xlab <- variable_names[x_var]
    }
    if(is.null(ylab)){
      ylab <- paste0(variable_names[y_var], "\nProportion (%)")
    }
  }
  
  
  if(method == "dodge"){
    
    if(is.null(legend_title_colors_bar)){
      legend_title_colors_bar <- variable_names[x_var]
    }
    
    if(is.null(xlab)){
      xlab <- variable_names[facet_var]
    }
    if(is.null(ylab)){
      ylab <- paste0(variable_names[y_var], ": ", keep_levels, "\nProportion (%)")
    }
  }
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions
  # --------------------------------------------------------------------------
  
  facet_levels <- levels(data[, facet_var])
  
  
  ggdata <- lapply(1:length(facet_levels), function(i){
    # i = 1
    
    data_sub <- data[data[, facet_var] == facet_levels[i] & !is.na(data[, facet_var]), ]
    
    if(nrow(data_sub) == 0){
      return(NULL)
    }
    
    
    tbl <- table(data_sub[, y_var], data_sub[, x_var])
    
    prop <- prop.table(tbl, margin = 2) * 100
    
    
    countdf <- as.data.frame.matrix(tbl)
    countdf <- data.frame(Observation = rownames(countdf), countdf, check.names = FALSE, row.names = NULL)
    
    propdf <- as.data.frame.matrix(prop)
    propdf <- data.frame(Observation = rownames(propdf), propdf, check.names = FALSE, row.names = NULL)
    
    
    ### Prepare data for ggplot
    
    ggdata <- cbind(gather(propdf, key = "Subgroup", value = "Proportion", -Observation), 
      gather(countdf, key = "Subgroup", value = "Count", -Observation)[, "Count", drop = FALSE])
    
    
    
    ### Prepare labels
    if(show_proportions){
      if(show_counts){
        ggdata$Label <- paste0(round(ggdata$Proportion, 1), "% (", ggdata$Count, ")")
      }else{
        ggdata$Label <- paste0(round(ggdata$Proportion, 1), "%")
      }
      
    }else{
      if(show_counts){
        ggdata$Label <- paste0("(", ggdata$Count, ")")
      }else{
        ggdata$Label <- NULL
      }
    }
    
    
    ### Remove data with zero counts so the ggplot warning with missing values is not displayed
    # ggdata <- ggdata[ggdata$Count > 0, ]
    
    ### Set label to NA for zero counts so it is not displayed
    ggdata$Label[ggdata$Count == 0] <- NA
    
    
    ggdata$Subgroup <- factor(ggdata$Subgroup, levels = levels(data_sub[, x_var]))
    
    ggdata$Observation <- factor(ggdata$Observation, levels = levels(data_sub[, y_var]))
    
    
    ggdata[, facet_var] <- facet_levels[i]
    
    
    return(ggdata)
    
    
  })
  
  
  ggdata <- plyr::rbind.fill(ggdata)
  
  
  ggdata[, facet_var] <- factor(ggdata[, facet_var], levels = levels(data[, facet_var]))
  
  
  # --------------------------------------------------------------------------
  # Calculate total proportions and counts
  # --------------------------------------------------------------------------
  
  
  sum_count_total <- aggregate(ggdata[, "Count"], lapply(c("Subgroup", facet_var), function(x) ggdata[, x]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count_total) <- c("Subgroup", facet_var, "Count_Total")
  
  ## For zero counts, there is NA. Let's fix it to zero
  sum_count_total[is.na(sum_count_total[, "Count_Total"]), "Count_Total"] <- 0
  
  
  ### Skip unwanted levels
  
  if(!is.null(skip_levels)){
    
    ggdata <- ggdata[ggdata$Observation %in% keep_levels, , drop = FALSE]
    
    ggdata$Observation <- factor(ggdata$Observation, keep_levels)
    
  }
  
  
  sum_prop <- aggregate(ggdata[, "Proportion"], lapply(c("Subgroup", facet_var), function(x) ggdata[, x]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_prop) <- c("Subgroup", facet_var, "Proportion")
  
  
  sum_count <- aggregate(ggdata[, "Count"], lapply(c("Subgroup", facet_var), function(x) ggdata[, x]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count) <- c("Subgroup", facet_var, "Count")
  
  
  ggdata_total <- sum_count_total %>% 
    left_join(sum_count, by = c("Subgroup", facet_var)) %>% 
    left_join(sum_prop, by = c("Subgroup", facet_var))
  
  
  
  ggdata_total$Subgroup <- factor(ggdata_total$Subgroup, levels = levels(ggdata$Subgroup))
  
  ggdata_total[, facet_var] <- factor(ggdata_total[, facet_var], levels = levels(ggdata[, facet_var]))
  
  ggdata_total$Label <- paste0(round(ggdata_total$Proportion, 1), "%")
  
  ggdata_total$Label_Total <- paste0("(", ggdata_total$Count_Total, ")")
  
  
  
  # --------------------------------------------------------------------------
  # Make the plot
  # --------------------------------------------------------------------------
  
  
  ### The default expand is 5% of the range, so let's nudge by 2.5%
  ### Calculate the place at the top where labels should be placed. At the bottom we place them at zero.
  if(!is.null(ylim)){
    ymax <- max(ylim)
  }else{
    ymax <- max(ggdata_total$Proportion, na.rm = TRUE)  
  }
  ynudge <- ymax * 0.025
  
  
  
  
  if(method == "facet"){
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = Subgroup, y = Proportion, fill = Observation), color = "black") +
      scale_fill_manual(name = legend_title_colors_bar, values = colors_bar, drop = FALSE) +
      theme_cowplot(12)
    
    
    ### Facet
    if(facet_var != "facet_dummy"){
      
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
          strip.text = element_text(size = strip_text_size))
      
    }
    
    
    
    ### Labels
    if(show_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata, aes(x = Subgroup, y = Proportion, group = Observation, label = Label), 
          position = position_stack(vjust = label_vjust), size = label_size)
      
    }
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = Subgroup, y = 0 - ynudge, label = Label_Total), size = label_size, vjust = 0.5)
      
    }
    
    if(show_total_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = Subgroup, y = Proportion + ynudge, label = Label), size = label_size, vjust = 0.5)
      
    }
    
  }
  
  
  
  
  if(method == "dodge"){
    
    ## ggplot2 doesn't know you want to give the labels the same virtual width as the bars. So tell it. 
    ## You can't nudge and dodge text, so instead adjust the y position.
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes_string(x = facet_var, y = "Proportion", fill = "Subgroup"), color = "black", position = position_dodge2(preserve = "single", width = 0.75)) +
      scale_fill_manual(name = legend_title_colors_bar, values = colors_bar, drop = FALSE) +
      theme_cowplot(12)
    
    
    ### Labels
    
    ## For dodge, zero total counts must be removed
    # ggdata_total <- ggdata_total[ggdata_total$Count_Total != 0, , drop = FALSE]
    
    ggdata_total$bottom_position <- 0 - ynudge
    ggdata_total$top_position <- ggdata_total$Proportion + ynudge
    
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes_string(x = facet_var, y = "bottom_position", group = "Subgroup", label = "Label_Total"), size = label_size, vjust = 0.5, position = position_dodge(preserve = "total", width = 0.9))
      
    }
    
    ### For one level left in the dodge setup show_proportions is equivalent to show_total_proportions
    if(show_proportions || show_total_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes_string(x = facet_var, y = "top_position", group = "Subgroup", label = "Label"), 
          size = label_size, vjust = 0.5, position = position_dodge(preserve = "total", width = 0.9))
      
    }
    
  }
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle, tag = tag) + 
    ylab(ylab) +
    xlab(xlab) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = title_size),
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      plot.tag.position = "top",
      plot.tag = element_text(size = title_size, face = "plain")) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.2) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = ylim)
  
  
  
  ggpl
  
  
}




# facet_var = NULL;
# strat1_var = NULL; strat2_var = NULL;
# colors_bar = NULL;
# variable_names = NULL;
# xlab = NULL; ylab = NULL; title = NULL; subtitle_label_both = TRUE; tag_label_both = TRUE;
# legend_title_colors_bar = NULL; facet_label_both = TRUE;
# skip_levels = NULL; method = "facet";
# show_proportions = TRUE; show_counts = TRUE; show_total_proportions = FALSE; show_total_counts = TRUE;
# label_size = 3.5; label_vjust = 0.5;
# title_size = 12; strip_text_size = NULL; facet_scales = "fixed"; ylim = NULL;
# axis_text_x_angle = 0; axis_text_x_vjust = 0; axis_text_x_hjust = 0.5;
# background_grid_major = "none";
# strat_scales = "fixed"; strat1_nrow = 1; strat1_ncol = NULL; strat2_nrow = NULL; strat2_ncol = 1




#' Barplot
#' 
#' Generate bar plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
wrapper_core_bar_plot_strat <- function(data, x_var, y_var, facet_var = NULL, 
  strat1_var = NULL, strat2_var = NULL,
  colors_bar = NULL, 
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors_bar = NULL, facet_label_both = TRUE, 
  skip_levels = NULL, method = "facet", 
  show_proportions = TRUE, show_counts = TRUE, show_total_proportions = FALSE, show_total_counts = TRUE, 
  label_size = 3.5, label_vjust = 0.5, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
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
  
  data <- data[complete.cases(data[, c(x_var, y_var, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Scales, xlim, ylim
  # -------------------------------------------------------------------------
  
  ### It is tricky to pre-calculate ylim for proportions before stratification. I would need to create a function for that part. So for now we set the fixed ylim to 0-100. 
  if(strat_scales == "fixed"){
    if(is.null(ylim)){
      ylim <- c(0, 100)
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
      
      
      ggpl <- wrapper_core_bar_plot(data = data_strata1, x_var = x_var, y_var = y_var, facet_var = facet_var, 
        colors_bar = colors_bar, 
        variable_names = variable_names, 
        xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, tag = tag, 
        legend_title_colors_bar = legend_title_colors_bar, facet_label_both = facet_label_both, 
        skip_levels = skip_levels, method = method, 
        show_proportions = show_proportions, show_counts = show_counts, show_total_proportions = show_total_proportions, show_total_counts = show_total_counts, 
        label_size = label_size, label_vjust = label_vjust, 
        title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, ylim = ylim, 
        axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust, 
        background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  ggpl
  
  
}






#' Barplots
#' 
#' Generate barplots for multiple variables in one faceted or dodged panel.
#' 
#' @param data Data frame.
wrapper_core_bar_plot_yvars_strat <- function(data, x_var, y_vars, 
  strat1_var = NULL, strat2_var = NULL,
  colors_bar = NULL, 
  variable_names = NULL, 
  xlab = NULL, ylab = NULL, title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors_bar = NULL, facet_label_both = TRUE, 
  skip_levels = NULL, method = "facet", 
  show_proportions = TRUE, show_counts = TRUE, show_total_proportions = FALSE, show_total_counts = TRUE, 
  label_size = 3.5, label_vjust = 0.5, 
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1,
  names_to = "name", values_to = "value"){
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  ### y_vars must be all factors
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "factor"))
  
  ### y_vars must have the same levels
  yvars_levels <- unique(unlist(lapply(data[, y_vars], levels)))
  stopifnot(all(unlist(lapply(data[, y_vars], function(x){levels(x) == yvars_levels}))))
  
  
  data <- data[, c(x_var, y_vars, strat1_var, strat2_var), drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Pivot to longer the data from y_vars
  # --------------------------------------------------------------------------
  
  data_longer <- tidyr::pivot_longer(data, y_vars, names_to = names_to, values_to = values_to) %>% 
    as.data.frame()
  
  data_longer[, names_to] <- factor(data_longer[, names_to], levels = y_vars, labels = variable_names[y_vars])
  data_longer[, values_to] <- factor(data_longer[, values_to], levels = yvars_levels)
  
  x_var <- x_var
  y_var <- values_to
  facet_var <- names_to
  
  
  facet_label_both <- FALSE
  
  
  ggpl <- wrapper_core_bar_plot_strat(data = data_longer, x_var = x_var, y_var = y_var, facet_var = facet_var, 
    strat1_var = strat1_var, strat2_var = strat2_var,
    colors_bar = colors_bar, 
    variable_names = variable_names, 
    xlab = xlab, ylab = ylab, title = title, subtitle_label_both = subtitle_label_both, tag_label_both = tag_label_both, 
    legend_title_colors_bar = legend_title_colors_bar, facet_label_both = facet_label_both, 
    skip_levels = skip_levels, method = method, 
    show_proportions = show_proportions, show_counts = show_counts, show_total_proportions = show_total_proportions, show_total_counts = show_total_counts, 
    label_size = label_size, label_vjust = label_vjust, 
    title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, ylim = ylim, 
    axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust, 
    background_grid_major = background_grid_major,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  ggpl
  
  
}














