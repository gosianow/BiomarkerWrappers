




# facet_var = NULL;
# colors_bar = NULL;
# variable_names = NULL;
# xlab = NULL; ylab = NULL; title = NULL; subtitle = NULL;
# legend_colors_title = NULL; facet_label_both = TRUE;
# skip_levels = NULL; method = "facet";
# show_proportions = TRUE; show_counts = TRUE; show_total_proportions = FALSE; show_total_counts = TRUE;
# label_size = 3.5; label_angle = 0;
# title_size = 12; strip_text_size = NULL; facet_scales = "fixed"; ylim = NULL;
# axis_text_x_angle = 0; axis_text_x_vjust = 0; axis_text_x_hjust = 0.5;
# background_grid_major = "none"




#' Barplot
#' 
#' Generate a signle barplot.
#' 
#' @param data Data frame.
#' 
#' @examples 
#'  
#' data(bdata)
#' 
#' data <- bdata
#' data$GeneA_cat2 <- wrapper_cut_2groups(data$GeneA)
#' 
#' x_var = "GeneA_cat2"
#' y_var = "Response"
#' 
#' wrapper_bar_plot_core(data = data, x_var = x_var, y_var = y_var)
#' 
#' @export
wrapper_bar_plot_core <- function(data, x_var, y_var, facet_var = NULL, 
  colors_bar = NULL, 
  variable_names = NULL, 
  title = TRUE, subtitle = TRUE, xlab = TRUE, ylab = TRUE,
  legend_colors_title = TRUE, legend_position = "right", facet_label_both = TRUE, 
  skip_levels = NULL, method = "facet", 
  show_proportions = TRUE, show_counts = TRUE, show_total_proportions = FALSE, show_total_counts = FALSE, 
  label_size = 3.5, label_angle = 0, label_nudge = 0.025,
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
  stopifnot(method %in% c("facet", "dodge_facet", "dodge_facet2", "facet2", "dodge"))
  
  
  ## Dodge method cannot be used when there is more than one level kept
  if(length(keep_levels) > 1 && method %in% c("dodge", "facet2")){
    stop("Method 'dodge' or 'facet2' cannot be used when there is more than one level kept")
  }
  
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(x_var, y_var, facet_var), drop = FALSE]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # Colors
  # -------------------------------------------------------------------------
  
  if(method %in% c("facet", "dodge_facet")){
    colors_bar <- format_colors(levels = levels(data[, y_var]), colors = colors_bar)
  }
  
  if(method %in% c("dodge", "facet2", "dodge_facet2")){
    colors_bar <- format_colors(levels = levels(data[, x_var]), colors = colors_bar)
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
  
  
  if(method %in% c("facet", "dodge_facet")){
    
    if(is.logical(legend_colors_title)){
      if(legend_colors_title){
        legend_colors_title <- variable_names[y_var]
      }else{
        legend_colors_title <- NULL
      }
    }
    
    if(is.logical(ylab)){
      if(ylab){
        ylab <- paste0(variable_names[y_var], "\nProportion (%)")
      }else{
        ylab <- NULL
      }
    }
    
    if(is.logical(xlab)){
      if(xlab){
        xlab <- variable_names[x_var]
      }else{
        xlab <- NULL
      }
    }
    
  }
  
  
  if(method %in% c("dodge_facet2")){
    
    if(is.logical(legend_colors_title)){
      if(legend_colors_title){
        legend_colors_title <- variable_names[x_var]
      }else{
        legend_colors_title <- NULL
      }
    }
    
    if(is.logical(ylab)){
      if(ylab){
        ylab <- paste0(variable_names[y_var], "\nProportion (%)")
      }else{
        ylab <- NULL
      }
    }
    
    if(is.logical(xlab)){
      if(xlab){
        xlab <- variable_names[y_var]
      }else{
        xlab <- NULL
      }
    }
    
  }
  
  
  if(method %in% c("facet2")){
    
    if(is.logical(legend_colors_title)){
      if(legend_colors_title){
        legend_colors_title <- variable_names[x_var]
      }else{
        legend_colors_title <- NULL
      }
    }
    
    if(is.logical(ylab)){
      if(ylab){
        ylab <- paste0(variable_names[y_var], ": ", keep_levels, "\nProportion (%)")
      }else{
        ylab <- NULL
      }
    }
    
    if(is.logical(xlab)){
      if(xlab){
        xlab <- variable_names[x_var]
      }else{
        xlab <- NULL
      }
    }
    
  }
  
  
  
  if(method == "dodge"){
    
    if(is.logical(legend_colors_title)){
      if(legend_colors_title){
        legend_colors_title <- variable_names[x_var]
      }else{
        legend_colors_title <- NULL
      }
    }
    
    if(is.logical(xlab)){
      if(xlab){
        xlab <- variable_names[facet_var]
      }else{
        xlab <- NULL
      }
    }
    
    if(is.logical(ylab)){
      if(ylab){
        ylab <- paste0(variable_names[y_var], ": ", keep_levels, "\nProportion (%)")
      }else{
        ylab <- NULL
      }
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
    
    ggdata <- dplyr::left_join(tidyr::pivot_longer(propdf, cols = colnames(tbl), names_to = "Subgroup", values_to = "Proportion"),
      tidyr::pivot_longer(countdf, cols = colnames(tbl), names_to = "Subgroup", values_to = "Count"), by = c("Observation", "Subgroup")) %>% 
      data.frame()
    
    
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
  
  
  sum_count_total <- stats::aggregate(ggdata[, "Count"], lapply(c("Subgroup", facet_var), function(x) ggdata[, x]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count_total) <- c("Subgroup", facet_var, "Count_Total")
  
  ## For zero counts, there is NA. Let's fix it to zero
  sum_count_total[is.na(sum_count_total[, "Count_Total"]), "Count_Total"] <- 0
  
  
  ### Skip unwanted levels
  
  if(!is.null(skip_levels)){
    
    ggdata <- ggdata[ggdata$Observation %in% keep_levels, , drop = FALSE]
    
    ggdata$Observation <- factor(ggdata$Observation, keep_levels)
    
  }
  
  
  sum_prop <- stats::aggregate(ggdata[, "Proportion"], lapply(c("Subgroup", facet_var), function(x) ggdata[, x]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_prop) <- c("Subgroup", facet_var, "Proportion")
  
  
  sum_count <- stats::aggregate(ggdata[, "Count"], lapply(c("Subgroup", facet_var), function(x) ggdata[, x]), FUN = sum, na.rm = TRUE, drop = FALSE)
  colnames(sum_count) <- c("Subgroup", facet_var, "Count")
  
  
  ggdata_total <- sum_count_total %>% 
    dplyr::left_join(sum_count, by = c("Subgroup", facet_var)) %>% 
    dplyr::left_join(sum_prop, by = c("Subgroup", facet_var))
  
  
  
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
    if(method %in% c("dodge_facet", "dodge_facet2")){
      ymax <- max(ggdata$Proportion, na.rm = TRUE)
    }else{
      ymax <- max(ggdata_total$Proportion, na.rm = TRUE)  
    }
  }
  ynudge <- ymax * label_nudge
  
  
  
  
  if(method == "facet"){
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = .data$Subgroup, y = .data$Proportion, fill = .data$Observation), color = "black") +
      scale_fill_manual(name = legend_colors_title, values = colors_bar, drop = FALSE)
    
    ### Labels
    if(show_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata, aes(x = .data$Subgroup, y = .data$Proportion, group = .data$Observation, label = .data$Label), 
          position = position_stack(vjust = 0.5), size = label_size, angle = label_angle)
      
    }
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        # geom_point(data = ggdata_total, aes(x = .data$Subgroup, y = 0 - ynudge)) + 
        geom_text(data = ggdata_total, aes(x = .data$Subgroup, y = 0 - ynudge, label = .data$Label_Total), size = label_size, angle = label_angle, vjust = 0.5)
      
    }
    
    if(show_total_proportions){
      
      ggpl <- ggpl +
        # geom_point(data = ggdata_total, aes(x = .data$Subgroup, y = .data$Proportion + ynudge)) +
        geom_text(data = ggdata_total, aes(x = .data$Subgroup, y = .data$Proportion + ynudge, label = .data$Label), size = label_size, angle = label_angle, vjust = 0.5)
      
    }
    
  }
  
  
  
  
  if(method == "dodge_facet"){
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = .data$Subgroup, y = .data$Proportion, fill = .data$Observation), color = "black", position = position_dodge2(preserve = "single", width = 0.75)) +
      scale_fill_manual(name = legend_colors_title, values = colors_bar, drop = FALSE)
    
    
    ### Labels
    if(show_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata, aes(x = .data$Subgroup, y = .data$Proportion + ynudge, group = .data$Observation, label = .data$Label), 
          position = position_dodge(preserve = "total", width = 0.9), size = label_size, angle = label_angle, vjust = 0.5)
      
    }
    
    
  }
  
  
  
  if(method == "dodge_facet2"){
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = .data$Observation, y = .data$Proportion, fill = .data$Subgroup), color = "black", position = position_dodge2(preserve = "single", width = 0.75)) +
      scale_fill_manual(name = legend_colors_title, values = colors_bar, drop = FALSE)
    
    
    ### Labels
    if(show_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata, aes(x = .data$Observation, y = .data$Proportion + ynudge, group = .data$Subgroup, label = .data$Label), 
          position = position_dodge(preserve = "total", width = 0.9), size = label_size, angle = label_angle, vjust = 0.5)
      
    }
    
    
  }
  
  
  
  if(method == "facet2"){
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = .data$Subgroup, y = .data$Proportion, fill = .data$Subgroup), color = "black") +
      scale_fill_manual(name = legend_colors_title, values = colors_bar, drop = FALSE)
    
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        # geom_point(data = ggdata_total, aes(x = .data$Subgroup, y = 0 - ynudge)) + 
        geom_text(data = ggdata_total, aes(x = .data$Subgroup, y = 0 - ynudge, label = .data$Label_Total), size = label_size, angle = label_angle, vjust = 0.5)
      
    }
    
    if(show_total_proportions){
      
      ggpl <- ggpl +
        # geom_point(data = ggdata_total, aes(x = .data$Subgroup, y = .data$Proportion + ynudge)) +
        geom_text(data = ggdata_total, aes(x = .data$Subgroup, y = .data$Proportion + ynudge, label = .data$Label), size = label_size, angle = label_angle, vjust = 0.5)
      
    }
    
  }
  
  
  
  
  if(method %in% c("facet", "dodge_facet", "dodge_facet2", "facet2")){
    
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
        facet_wrap(stats::as.formula(paste("~", facet_var)), labeller = labeller, scales = facet_scales, nrow = 1) +
        theme(strip.background = element_rect(colour = "white", fill = "white"),
          strip.text = element_text(size = strip_text_size))
      
    }
    
  }
  
  
  
  
  if(method == "dodge"){
    
    
    ## ggplot2 doesn't know you want to give the labels the same virtual width as the bars. So tell it. 
    ## You can't nudge and dodge text, so instead adjust the y position.
    
    
    ggpl <- ggplot() +
      geom_col(data = ggdata, aes(x = .data[[facet_var]], y = .data$Proportion, fill = .data$Subgroup), color = "black", position = position_dodge2(preserve = "single", width = 0.75)) +
      scale_fill_manual(name = legend_colors_title, values = colors_bar, drop = FALSE)
    
    
    ### Labels
    
    ## For dodge, zero total counts must be removed
    # ggdata_total <- ggdata_total[ggdata_total$Count_Total != 0, , drop = FALSE]
    
    ggdata_total$bottom_position <- 0 - ynudge
    ggdata_total$top_position <- ggdata_total$Proportion + ynudge
    
    
    if(show_total_counts){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = .data[[facet_var]], y = .data$bottom_position, group = .data$Subgroup, label = .data$Label_Total), size = label_size, angle = label_angle, vjust = 0.5, position = position_dodge(preserve = "total", width = 0.9))
      
    }
    
    ### For one level left in the dodge setup show_proportions is equivalent to show_total_proportions
    if(show_proportions || show_total_proportions){
      
      ggpl <- ggpl +
        geom_text(data = ggdata_total, aes(x = .data[[facet_var]], y = .data$top_position, group = .data$Subgroup, label = .data$Label), 
          size = label_size, angle = label_angle, vjust = 0.5, position = position_dodge(preserve = "total", width = 0.9))
      
    }
    
    
    
  }
  
  
  ggpl <- ggpl +
    labs(title = title, subtitle = subtitle) + 
    ylab(ylab) +
    xlab(xlab) +
    theme(plot.title = element_text(size = title_size, face = "bold"),
      plot.subtitle = element_text(size = title_size),
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      # axis.line = element_blank(),
      # axis.ticks = element_line(color = "black", size = 0.5),
      # panel.border = element_rect(colour = "black", size = 1),
      legend.position = legend_position) +
    background_grid(major = background_grid_major, minor = "none", size.major = 0.15) +
    scale_x_discrete(drop = FALSE) +
    coord_cartesian(ylim = ylim)
  
  
  
  ggpl
  
  
}




# facet_var = NULL;
# strat1_var = NULL; strat2_var = NULL;
# colors_bar = NULL;
# variable_names = NULL;
# xlab = NULL; ylab = NULL; title = NULL; strat1_label_both = TRUE; strat2_label_both = TRUE;
# legend_colors_title = NULL; facet_label_both = TRUE;
# skip_levels = NULL; method = "facet";
# show_proportions = TRUE; show_counts = TRUE; show_total_proportions = FALSE; show_total_counts = TRUE;
# label_size = 3.5; label_angle = 0;
# title_size = 12; strip_text_size = NULL; facet_scales = "fixed"; ylim = NULL;
# axis_text_x_angle = 0; axis_text_x_vjust = 0; axis_text_x_hjust = 0.5;
# background_grid_major = "none";
# strat_scales = "fixed"; strat1_nrow = 1; strat1_ncol = NULL; strat2_nrow = NULL; strat2_ncol = 1





#' @rdname wrapper_bar_plot_core
#' @param strat1_var Name of the first stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' @export
wrapper_bar_plot_core_strat <- function(data, x_var, y_var, facet_var = NULL, 
  strat1_var = NULL, strat2_var = NULL,
  colors_bar = NULL, 
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, ylab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_title = TRUE, legend_position = "right", facet_label_both = TRUE, 
  skip_levels = NULL, method = "facet", 
  show_proportions = TRUE, show_counts = TRUE, show_total_proportions = FALSE, show_total_counts = FALSE, 
  label_size = 3.5, label_angle = 0, label_nudge = 0.025,
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
  
  data <- data[stats::complete.cases(data[, c(x_var, y_var, strat1_var, strat2_var)]), ]
  
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
      
      ggpl <- wrapper_bar_plot_core(data = data_strata1, x_var = x_var, y_var = y_var, facet_var = facet_var, 
        colors_bar = colors_bar, 
        variable_names = variable_names, 
        xlab = xlab, ylab = ylab, title = title, subtitle = subtitle, 
        legend_colors_title = legend_colors_title, legend_position = legend_position, facet_label_both = facet_label_both, 
        skip_levels = skip_levels, method = method, 
        show_proportions = show_proportions, show_counts = show_counts, show_total_proportions = show_total_proportions, show_total_counts = show_total_counts, 
        label_size = label_size, label_angle = label_angle, label_nudge = label_nudge,
        title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, ylim = ylim, 
        axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust, 
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






#' Barplots
#' 
#' Generate barplots for multiple variables in one faceted or dodged panel.
#' 
#' @inheritParams wrapper_bar_plot_core_strat
#' @export
wrapper_bar_plot_yvars_core_strat <- function(data, x_var, y_vars, 
  strat1_var = NULL, strat2_var = NULL,
  colors_bar = NULL, 
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, ylab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_title = TRUE, legend_position = "right", facet_label_both = TRUE, 
  skip_levels = NULL, method = "facet", 
  show_proportions = TRUE, show_counts = TRUE, show_total_proportions = FALSE, show_total_counts = FALSE, 
  label_size = 3.5, label_angle = 0, label_nudge = 0.025,
  title_size = 12, strip_text_size = NULL, facet_scales = "fixed", ylim = NULL, 
  axis_text_x_angle = 0, axis_text_x_vjust = 0, axis_text_x_hjust = 0.5, 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1, less_legends = FALSE,
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
  
  
  ggpl <- wrapper_bar_plot_core_strat(data = data_longer, x_var = x_var, y_var = y_var, facet_var = facet_var, 
    strat1_var = strat1_var, strat2_var = strat2_var,
    colors_bar = colors_bar, 
    variable_names = variable_names, 
    xlab = xlab, ylab = ylab, title = title, strat1_label_both = strat1_label_both, strat2_label_both = strat2_label_both, 
    legend_colors_title = legend_colors_title, legend_position = legend_position, facet_label_both = facet_label_both, 
    skip_levels = skip_levels, method = method, 
    show_proportions = show_proportions, show_counts = show_counts, show_total_proportions = show_total_proportions, show_total_counts = show_total_counts, 
    label_size = label_size, label_angle = label_angle, label_nudge = label_nudge, 
    title_size = title_size, strip_text_size = strip_text_size, facet_scales = facet_scales, ylim = ylim, 
    axis_text_x_angle = axis_text_x_angle, axis_text_x_vjust = axis_text_x_vjust, axis_text_x_hjust = axis_text_x_hjust, 
    background_grid_major = background_grid_major,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol, less_legends = less_legends)
  
  
  ggpl
  
  
}







