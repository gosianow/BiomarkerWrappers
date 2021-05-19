




#' KM plot
#' 
#' @param data Data frame.
#' @param tte_var Name of the variable containing time-to-event data.
#' @param censor_var Name of the variable containing censoring information. Censor variable must be numeric and encode 1 for event and 0 for censor.
#' @param covariate_var Name of variable that defines the subgroups where the survival is calculated. This variable must be a factor.
#' @param level_mapping Named vector with level mapping. The names correspond to the original levels of covariate_var. 
#' 
#' @examples 
#' 
#' data(bdata)
#' 
#' data <- bdata
#' 
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' covariate_var <- "Treatment_Arm"
#' 
#' 
#' wrapper_KM_plot_core(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var)
#' 
#' @export
wrapper_KM_plot_core <- function(data, tte_var, censor_var, covariate_var, 
  colors = NULL, linetypes = 1, 
  variable_names = NULL, 
  title = TRUE, subtitle = TRUE, xlab = TRUE,
  legend_colors_title = TRUE, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, surv_median_line = "none",
  ggtheme = ggplot2::theme_classic(12), 
  line_size = 1, title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none", level_mapping = NULL){
  
  
  # -------------------------------------------------------------------------
  # Checks about data
  # -------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_var)]), ]
  
  stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)

  
  ## Time to event variable must be numeric
  stopifnot(length(tte_var) == 1)
  stopifnot(is.numeric(data[, tte_var]))
  
  ## Censor variable must be numeric and encode 1 for event and 0 for censor
  stopifnot(length(censor_var) == 1)
  stopifnot(is.numeric(data[, censor_var]) && all(data[, censor_var] %in% c(0, 1)))
  
  stopifnot(length(covariate_var) == 1)
  stopifnot(is.factor(data[, covariate_var]))
  
  
  # -------------------------------------------------------------------------
  # Colors
  # -------------------------------------------------------------------------
  
  colors <- format_colors(levels(data[, covariate_var]), colors = colors, allow_duplicated = TRUE)
  
  ### Because colors are taken in a row from the beginning of the vector to have consistent coloring we have to remove colors for the levels with zero counts. For the ggsurvplot function and in ggplot adjustment colors cannot have names. Otherwise, it does not work. 
  
  tbl <- table(data[, covariate_var])
  
  colors <- colors[tbl != 0]
  
  colors <- as.character(colors)
  
  # -------------------------------------------------------------------------
  # Line types
  # -------------------------------------------------------------------------
  
  
  stopifnot(length(linetypes) == 1 || length(linetypes) == nlevels(data[, covariate_var]))
  
  if(length(linetypes) == 1){
    linetypes <- rep(linetypes, nlevels(data[, covariate_var]))
  }
  
  linetypes <- linetypes[tbl != 0]
  
  
  ## Remove unused levels from the data too
  
  data[, covariate_var] <- factor(data[, covariate_var])
  
  
  ## Map levels 
  
  if(!is.null(level_mapping)){
    
    data[, covariate_var] <- factor(data[, covariate_var], labels = level_mapping[levels(data[, covariate_var])])
    
  }
  
  
  # -------------------------------------------------------------------------
  # Labels
  # -------------------------------------------------------------------------
  
  if(is.logical(title)){
    if(title){
      title <- variable_names[covariate_var]
    }else{
      title <- NULL
    }
  }
  
  if(is.logical(subtitle)){
    subtitle <- NULL
  }
  
  if(is.logical(xlab)){
    if(xlab){
      xlab <- variable_names[tte_var]
    }else{
      xlab <- NULL
    }
  }
  
  if(is.logical(legend_colors_title)){
    if(legend_colors_title){
      legend_colors_title <- variable_names[covariate_var]
    }else{
      legend_colors_title <- NULL
    }
  }
  
  
  # -------------------------------------------------------------------------
  # break_time_by and max_tte
  # -------------------------------------------------------------------------
  
  
  ### The default is to have about 10 breaks
  
  if(is.null(break_time_by)){
    break_time_by <- calculate_break_time(data[, tte_var], n_breaks = 10)
  }
  # print(break_time_by)
  
  if(is.null(max_tte)){
    max_tte <- max(data[, tte_var], na.rm = TRUE)
  }
  # print(max_tte)
  
  ### To make sure that no data is cut off the range of the plot, extend the x-axis
  ## Extend to the next break time point - It is usually too much and gives zeros in risk table. 
  # max_tte <- ceiling(max_tte / break_time_by) * break_time_by
  ## Extend by 5% of the break time
  max_tte <- max_tte + break_time_by * 0.05
  
  
  # -------------------------------------------------------------------------
  # Fit the survival model
  # -------------------------------------------------------------------------
  
  
  ### Define the model formula
  f <- stats::as.formula(paste0("Surv(", tte_var, ",", censor_var,") ~ ", covariate_var))
  
  ### Fit the model
  fit <- survival::survfit(f, data, conf.type = "plain")
  
  ## Overwrite the formula. Otherwise, it does not work!!! 
  fit$call$formula <- f
  
  
  # -------------------------------------------------------------------------
  # Generate Kaplan-Meier plot
  # -------------------------------------------------------------------------
  
  
  ## palette must be a non-named vector. Otherwise, it does not work. For each subplot has to have unique values. If a level has zero counts, it is not plotted. Because colors are taken in a row from the beginning of the vector to have consistent coloring we have to remove colors for the levels with zero counts.
  
  ### Use that trick until the issue https://github.com/kassambara/survminer/issues/519 is resolved 
  
  if(sum(duplicated(colors)) > 0){
    
    ggpl_plot <- survminer::ggsurvplot(fit, data = data, palette = colors, linetype = linetypes, conf.int = conf_int, surv.median.line = surv_median_line, risk.table = FALSE, ggtheme = ggtheme, xlab = xlab, break.time.by = break_time_by, xlim = c(0, max_tte), fontsize = label_size, size = line_size) 
    
    ggpl_table <- survminer::ggsurvplot(fit, data = data, palette = rep("black", length(colors)), linetype = linetypes, conf.int = conf_int, surv.median.line = surv_median_line, risk.table = TRUE, ggtheme = ggtheme, xlab = xlab, break.time.by = break_time_by, xlim = c(0, max_tte), fontsize = label_size, size = line_size)
    
    ggpl <- list()
    
    ggpl$plot <- ggpl_plot$plot
    ggpl$table <- ggpl_table$table
    
  }else{
    
    ggpl <- survminer::ggsurvplot(fit, data = data, palette = colors, linetype = linetypes, conf.int = conf_int, surv.median.line = surv_median_line, risk.table = risk_table, ggtheme = ggtheme, xlab = xlab, break.time.by = break_time_by, xlim = c(0, max_tte), fontsize = label_size, size = line_size) 
    
  }
  
  
  ### Customize the plot
  suppressMessages(ggpl_plot <- ggpl$plot +
      labs(title = title, subtitle = subtitle) +
      theme(plot.title = element_text(size = title_size, face = "bold"),
        plot.subtitle = element_text(size = title_size),
        legend.position = legend_position,
        legend.justification = legend_justification,
        legend.background = element_rect(fill = NA)) +
      scale_color_manual(name = legend_colors_title, labels = levels(data[, covariate_var]), values = colors) +
      scale_fill_manual(name = legend_colors_title, labels = levels(data[, covariate_var]), values = colors) +
      scale_linetype_manual(name = legend_colors_title, labels = levels(data[, covariate_var]), values = linetypes) +
      coord_cartesian(xlim = c(0, max_tte)) +
      background_grid(major = background_grid_major, minor = "none", size.major = 0.15))
  
  
  if(risk_table){
    
    suppressMessages(ggpl_table <- ggpl$table +
        theme(plot.title = element_text(size = title_size),
          axis.title = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks = element_blank(), 
          axis.line = element_blank()) +
        scale_y_discrete(labels = rev(levels(data[, covariate_var]))) +
        coord_cartesian(xlim = c(0, max_tte)))
    
    
    ggpl_new <- cowplot::plot_grid(ggpl_plot, ggpl_table, ncol = 1, nrow = 2, align = 'v', axis = 'lr', rel_heights = rel_heights)
    
    
  }else{
    
    ggpl_new <- ggpl_plot
    
  }
  
  
  ggpl_new
  
  
}






#' @rdname wrapper_KM_plot_core
#' @param strat1_var Name of the firts stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' @examples 
#' 
#' data(bdata)
#' 
#' data <- bdata
#' data$GeneA_cat2 <- wrapper_cut_2groups(data$GeneA)
#' 
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' covariate_var <- "GeneA_cat2"
#' 
#' strat1_var = "Treatment_Arm"
#' strat2_var = "Cell_Of_Origin"
#' 
#' 
#' wrapper_KM_plot_core_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, strat1_var = strat1_var, strat2_var = strat2_var)
#' 
#' @export
wrapper_KM_plot_core_strat <- function(data, tte_var, censor_var, covariate_var, 
  strat1_var = NULL, strat2_var = NULL,
  colors = NULL, linetypes = 1, 
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_title = TRUE, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, surv_median_line = "none",
  ggtheme = ggplot2::theme_classic(12),
  line_size = 1, title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none", level_mapping = NULL,
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
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_var, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # -------------------------------------------------------------------------
  # break_time_by and max_tte
  # -------------------------------------------------------------------------
  
  
  ### The default is to have about 10 breaks
  
  if(is.null(break_time_by)){
    break_time_by <- calculate_break_time(data[, tte_var], n_breaks = 10)
  }
  
  
  if(strat_scales == "fixed"){
    if(is.null(max_tte)){
      max_tte <- max(data[, tte_var], na.rm = TRUE)
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
      
      
      ggpl <- wrapper_KM_plot_core(data = data_strata1, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, 
        colors = colors, linetypes = linetypes, 
        variable_names = variable_names, 
        title = title, subtitle = subtitle, xlab = xlab,
        legend_colors_title = legend_colors_title, legend_position = legend_position, legend_justification = legend_justification,
        break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, surv_median_line = surv_median_line,
        ggtheme = ggtheme,
        line_size = line_size, title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
        background_grid_major = background_grid_major, level_mapping = level_mapping)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  
  ggpl
  
  
}








#' KM plot with curves per biomarker and treatment in a single panel
#' 
#' @inheritParams wrapper_KM_plot_core_strat
#' @param colors Vector with colors for treatment X biomarker levels.
#' @export
wrapper_KM_plot_interaction <- function(data, tte_var, censor_var, biomarker_var, treatment_var, 
  strat1_var = NULL, strat2_var = NULL,
  colors = NULL, linetypes = 1, 
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_title = TRUE, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, surv_median_line = "none",
  ggtheme = ggplot2::theme_classic(12),
  line_size = 1, title_size = 12, label_size = 3, rel_heights = c(4, 1), 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  
  ## biomarker_var and treatment_var must be factors for the color definition
  stopifnot(length(biomarker_var) == 1)
  stopifnot(is.factor(data[, biomarker_var]))
  stopifnot(length(treatment_var) == 1)
  stopifnot(is.factor(data[, treatment_var]))
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(is.logical(title)){
    if(title){
      title <- variable_names[biomarker_var]
    }else{
      title <- NULL
    }
  }
  
  
  # -------------------------------------------------------------------------
  # Generate the treatment-biomarker interaction covariate
  # -------------------------------------------------------------------------
  
  stopifnot(!"treatment_biomarker_interaction" %in% colnames(data))
  
  data$treatment_biomarker_interaction <- interaction(data[, treatment_var], data[, biomarker_var], lex.order = TRUE, sep = ", ")
  
  covariate_var <- "treatment_biomarker_interaction"
  
  variable_names[[covariate_var]] <- structure(paste0(variable_names[treatment_var], ", ", variable_names[biomarker_var]), names = covariate_var)
  
  
  # -------------------------------------------------------------------------
  # Colors
  # -------------------------------------------------------------------------
  
  if(is.null(colors)){
    
    colors <- format_colors_cat_strata(levels(data[, biomarker_var]), strata = levels(data[, treatment_var]))
    
    # barplot(rep(1, length(colors)), col = colors)
    
  }else{
    
    colors <- format_colors(levels(data[, covariate_var]), colors = colors, allow_duplicated = TRUE)
    
  }
  
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  
  ggpl <- wrapper_KM_plot_core_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var,
    strat1_var = strat1_var, strat2_var = strat2_var, 
    colors = colors, linetypes = linetypes, 
    variable_names = variable_names, 
    title = title, xlab = xlab, strat1_label_both = strat1_label_both, strat2_label_both = strat2_label_both, 
    legend_colors_title = legend_colors_title, legend_position = legend_position, legend_justification = legend_justification,
    break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, surv_median_line = surv_median_line,
    ggtheme = ggtheme,
    line_size = line_size, title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
    background_grid_major = background_grid_major,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  return(ggpl)
  
  
}





#' KM plots with biomarker effect per treatment arm
#' 
#' @inheritParams wrapper_KM_plot_interaction
#' @param colors Vector with colors for treatment X biomarker levels.
#' @export
wrapper_KM_plot_biomarker <- function(data, tte_var, censor_var, biomarker_var, treatment_var = NULL, 
  strat2_var = NULL,
  colors = NULL, linetypes = 1, 
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_title = TRUE, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, surv_median_line = "none",
  ggtheme = ggplot2::theme_classic(12),
  line_size = 1, title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){


  ## biomarker_var and treatment_var must be factors for the color definition
  stopifnot(length(biomarker_var) == 1)
  stopifnot(is.factor(data[, biomarker_var]))
  
  if(!is.null(treatment_var)){
    stopifnot(length(treatment_var) == 1)
    stopifnot(is.factor(data[, treatment_var]))
  }
  
  

  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(is.logical(title)){
    if(title){
      title <- variable_names[biomarker_var]
    }else{
      title <- NULL
    }
  }
  
  
  if(!is.null(treatment_var)){
    
    
    # -------------------------------------------------------------------------
    # Generate the treatment-biomarker interaction covariate
    # -------------------------------------------------------------------------
    
    stopifnot(!"treatment_biomarker_interaction" %in% colnames(data))
    
    data$treatment_biomarker_interaction <- interaction(data[, treatment_var], data[, biomarker_var], lex.order = TRUE, sep = ", ")
    
    covariate_var <- "treatment_biomarker_interaction"
    
    variable_names[[covariate_var]] <- structure(variable_names[biomarker_var], names = covariate_var)
    
    
    level_mapping <- rep(levels(data[, biomarker_var]), times = nlevels(data[, treatment_var]))
    names(level_mapping) <- levels(data$treatment_biomarker_interaction)
    
    
    # -------------------------------------------------------------------------
    # Colors
    # -------------------------------------------------------------------------
    
    
    if(is.null(colors)){
      
      colors <- format_colors_cat_strata(levels(data[, biomarker_var]), strata = levels(data[, treatment_var]))
      
      # barplot(rep(1, length(colors)), col = colors)
      
    }else{
      
      colors <- format_colors(levels(data[, covariate_var]), colors = colors, allow_duplicated = TRUE)
      
    }
    
    strat1_var <- treatment_var
    
    
  }else{
    
    covariate_var <- biomarker_var
    
    colors <- format_colors(levels(data[, covariate_var]), colors = colors, allow_duplicated = TRUE)
    
    strat1_var <- NULL
    
    
  }
  
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  
  ggpl <- wrapper_KM_plot_core_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var,
    strat1_var = strat1_var, strat2_var = strat2_var, 
    colors = colors, linetypes = linetypes, 
    variable_names = variable_names, 
    title = title, xlab = xlab, strat1_label_both = strat1_label_both, strat2_label_both = strat2_label_both, 
    legend_colors_title = legend_colors_title, legend_position = legend_position, legend_justification = legend_justification,
    break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, surv_median_line = surv_median_line,
    ggtheme = ggtheme,
    line_size = line_size, title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
    background_grid_major = background_grid_major, level_mapping = level_mapping,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  
  ggpl
  
  
}












#' KM plot with treatment effect per biomarker subgroup
#' 
#' @inheritParams wrapper_KM_plot_interaction
#' @param colors Vector with colors for treatment X biomarker levels.
#' @export
wrapper_KM_plot_treatment <- function(data, tte_var, censor_var, treatment_var, biomarker_var = NULL,
  strat2_var = NULL,
  colors = NULL, linetypes = 1,
  variable_names = NULL, 
  title = TRUE, xlab = TRUE, strat1_label_both = TRUE, strat2_label_both = TRUE, 
  legend_colors_title = TRUE, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, surv_median_line = "none",
  ggtheme = ggplot2::theme_classic(12),
  line_size = 1, title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  ## biomarker_var and treatment_var must be factors for the color definition
  if(!is.null(biomarker_var)){
    stopifnot(length(biomarker_var) == 1)
    stopifnot(is.factor(data[, biomarker_var]))
  }
  
  stopifnot(length(treatment_var) == 1)
  stopifnot(is.factor(data[, treatment_var]))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(is.logical(title)){
    if(title){
      title <- variable_names[biomarker_var]
    }else{
      title <- NULL
    }
  }
  
  if(!is.null(biomarker_var)){
    
    # -------------------------------------------------------------------------
    # Generate the treatment-biomarker interaction covariate
    # -------------------------------------------------------------------------
    
    stopifnot(!"treatment_biomarker_interaction" %in% colnames(data))
    
    data$treatment_biomarker_interaction <- interaction(data[, treatment_var], data[, biomarker_var], lex.order = TRUE, sep = ", ")
    
    covariate_var <- "treatment_biomarker_interaction"
    
    variable_names[[covariate_var]] <- structure(variable_names[treatment_var], names = covariate_var)
    
    
    level_mapping <- rep(levels(data[, treatment_var]), each = nlevels(data[, biomarker_var]))
    names(level_mapping) <- levels(data$treatment_biomarker_interaction)
    
    
    # -------------------------------------------------------------------------
    # Colors
    # -------------------------------------------------------------------------
    
    
    if(is.null(colors)){
      
      # colors <- format_colors_cat_strata(levels(data[, biomarker_var]), strata = levels(data[, treatment_var]))
      
      colors <- format_colors_cat_strata("one_level", strata = levels(data[, treatment_var]))
      
      colors <- rep(colors, each = nlevels(data[, biomarker_var]))
      
      names(colors) <- levels(data$treatment_biomarker_interaction)
      
      # barplot(rep(1, length(colors)), col = colors)
      
    }else{
      
      colors <- format_colors(levels(data[, covariate_var]), colors = colors, allow_duplicated = TRUE)
      
    }
    
    
    strat1_var <- biomarker_var
    
  }else{
    
    covariate_var <- treatment_var
    
    colors <- format_colors(levels(data[, covariate_var]), colors = colors, allow_duplicated = TRUE)
    
    strat1_var <- NULL
    
  }
  
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  
  ggpl <- wrapper_KM_plot_core_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var,
    strat1_var = strat1_var, strat2_var = strat2_var, 
    colors = colors, linetypes = linetypes, 
    variable_names = variable_names, 
    title = title, xlab = xlab, strat1_label_both = strat1_label_both, strat2_label_both = strat2_label_both, 
    legend_colors_title = legend_colors_title, legend_position = legend_position, legend_justification = legend_justification,
    break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, surv_median_line = surv_median_line,
    ggtheme = ggtheme,
    line_size = line_size, title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
    background_grid_major = background_grid_major, level_mapping = level_mapping, 
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  
  ggpl
  
  
}




























