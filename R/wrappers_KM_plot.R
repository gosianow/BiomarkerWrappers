




#' KM plot
#' 
#' Generate a signle KM plot.
#' 
#' @param data Data frame.
#' @export
wrapper_core_KM_plot <- function(data, tte_var, censor_var, covariate_var, 
  colors = NULL, 
  variable_names = NULL, 
  title = NULL, subtitle = NULL, tag = NULL, 
  legend_title_colors = NULL, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, 
  title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none"){
  
  
  # -------------------------------------------------------------------------
  # Checks about data
  # -------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_var)]), ]
  
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
  
  colors <- format_colors(levels = levels(data[, covariate_var]), colors = colors, allow_duplicated = FALSE)
  
  ### Because colors are taken in a row from the beginning of the vector to have consistent coloring we have to remove colors for the levels with zero counts. For the ggsurvplot function and in ggplot adjustment colors cannot have names. Otherwise, it does not work. 
  
  tbl <- table(data[, covariate_var])
  
  colors <- colors[tbl != 0]
  
  colors <- as.character(colors)
  
  ## Removed unused levels from the data too
  
  data[, covariate_var] <- factor(data[, covariate_var])
  
  
  # -------------------------------------------------------------------------
  # Labels
  # -------------------------------------------------------------------------
  
  
  xlab <- variable_names[tte_var]
  
  
  if(is.null(legend_title_colors)){
    legend_title_colors <- variable_names[covariate_var]
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
  fit <- survival::survfit(f, data)
  
  ## Fix a bug. Otherwise, it does not work!!! 
  fit$call$formula <- f
  
  
  # -------------------------------------------------------------------------
  # Generate Kaplan-Meier plot
  # -------------------------------------------------------------------------
  
  
  ## palette must be a non-named vector. Otherwise, it does not work. For each subplot has to have unique values. If a level has zero counts, it is not plotted. Because colors are taken in a row from the beginning of the vector to have consistent coloring we have to remove colors for the levels with zero counts.
  
  ggpl <- survminer::ggsurvplot(fit, data = data, palette = colors, linetype = 1, conf.int = conf_int, risk.table = risk_table, ggtheme = theme_classic(), xlab = xlab, break.time.by = break_time_by, xlim = c(0, max_tte), fontsize = label_size) 
  
  
  ### Customize the plot
  suppressMessages(ggpl_plot <- ggpl$plot +
      labs(title = title, subtitle = subtitle, tag = tag) +
      theme(plot.title = element_text(size = title_size, face = "bold"),
        plot.subtitle = element_text(size = title_size),
        legend.position = legend_position,
        legend.justification = legend_justification,
        legend.background = element_rect(fill = NA),
        plot.tag.position = "top",
        plot.tag = element_text(size = title_size, face = "plain")) +
      scale_color_manual(name = legend_title_colors, labels = levels(data[, covariate_var]), values = colors) +
      scale_fill_manual(name = legend_title_colors, labels = levels(data[, covariate_var]), values = colors) +
      coord_cartesian(xlim = c(0, max_tte)) +
      background_grid(major = background_grid_major, minor = "none", size.major = 0.15))
  
  
  if(risk_table){
    
    suppressMessages(ggpl_table <- ggpl$table +
        theme(plot.title = element_text(size = title_size),
          axis.title = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks = element_blank(), 
          axis.line=element_blank()) +
        scale_y_discrete(labels = rev(levels(data[, covariate_var]))) +
        coord_cartesian(xlim = c(0, max_tte)))
    
    
    ggpl_new <- cowplot::plot_grid(ggpl_plot, ggpl_table, ncol = 1, nrow = 2, align = 'v', axis = 'lr', rel_heights = rel_heights)
    
    
  }else{
    
    ggpl_new <- ggpl_plot
    
  }
  
  
  ggpl_new
  
  
}











#' KM plot
#' 
#' Generate KM plots for each subgroup defined by two stratification variables.
#' 
#' @param data Data frame.
#' @export
wrapper_core_KM_plot_strat <- function(data, tte_var, censor_var, covariate_var, 
  strat1_var = NULL, strat2_var = NULL,
  colors = NULL, 
  variable_names = NULL, 
  title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors = NULL, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, 
  title_size = 12, label_size = 3, rel_heights = c(5, 1), 
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
      
      
      ggpl <- wrapper_core_KM_plot(data = data_strata1, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, 
        colors = colors, 
        variable_names = variable_names, 
        title = title, subtitle = subtitle, tag = tag, 
        legend_title_colors = legend_title_colors, legend_position = legend_position, legend_justification = legend_justification,
        break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, 
        title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
        background_grid_major = background_grid_major)
      
      
      return(ggpl)
      
      
    })
    
    
    ggpl <- plot_grid(plotlist = ggpl, nrow = strat1_nrow, ncol = strat1_ncol)
    
    return(ggpl)
    
  })
  
  
  ggpl <- plot_grid(plotlist = ggpl, nrow = strat2_nrow, ncol = strat2_ncol)
  
  
  
  ggpl
  
  
}








#' KM plot - All curves per biomarker and treatment in one panel
#' 
#' @param data Data frame.
#' @param colors A list of length equal to treatment levels.
#' @export
wrapper_KM_plot_interaction <- function(data, tte_var, censor_var, biomarker_var, treatment_var, 
  strat1_var = NULL, strat2_var = NULL,
  colors = NULL, 
  variable_names = NULL, 
  title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors = NULL, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, 
  title_size = 12, label_size = 3, rel_heights = c(4, 1), 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(is.null(title)){
    title <- variable_names[biomarker_var]
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
  
  
  nlevels_biomarker <- nlevels(data[, biomarker_var])
  levels_biomarker <- levels(data[, biomarker_var])
  nlevels_treatment <- nlevels(data[, treatment_var])
  levels_treatment <- levels(data[, treatment_var])
  
  
  if(is.null(colors)){
    ## Some default colors that work for max 4 treatment levels
    
    if(nlevels_treatment > 4){
      stop("There are no default colors available when number of treatment levels is higher than 4. Please, provide the colors.")
    }
    
    default_colors_per_treatment <- list(c("#E69F00", "#B22222"), c("#56B4E9", "#104E8B"), c("#84d44b", "#007800"), c("#836FFF", "#50449d"))
    
    # palette <- unlist(default_colors_per_treatment)
    # barplot(rep(1, length(palette)), col = palette)
    
    
    colors <- unlist(lapply(1:nlevels_treatment, function(i){
      format_colors(levels = paste0(levels_treatment[i], ", ", levels_biomarker), palette = default_colors_per_treatment[[i]], allow_duplicated = FALSE)
    }))
    
    
  }else{
    
    colors <- unlist(lapply(1:nlevels_treatment, function(i){
      # i = 1
      out <- format_colors(levels = levels_biomarker, colors = colors[[i]], allow_duplicated = FALSE)
      names(out) <- paste0(levels_treatment[i], ", ", levels_biomarker)
      return(out)
    }))
    
  }
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  
  ggpl <- wrapper_core_KM_plot_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var,
    strat1_var = strat1_var, strat2_var = strat2_var, 
    colors = colors, 
    variable_names = variable_names, 
    title = title, subtitle_label_both = subtitle_label_both, tag_label_both = tag_label_both, 
    legend_title_colors = legend_title_colors, legend_position = legend_position, legend_justification = legend_justification,
    break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, 
    title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
    background_grid_major = background_grid_major,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  return(ggpl)
  
  
}





#' KM plot - Biomarker effect per treatment arm
#' 
#' @param data Data frame.
#' @param colors If treatment = NULL, then a vector. Otherwise, a list of length equal to tratment levels.
#' @export
wrapper_KM_plot_biomarker <- function(data, tte_var, censor_var, biomarker_var, treatment_var = NULL, 
  strat2_var = NULL,
  colors = NULL, 
  variable_names = NULL, 
  title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors = NULL, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, 
  title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  ### TODO Do not display treatment in the legend
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(is.null(title)){
    title <- variable_names[biomarker_var]
  }
  
  
  if(!is.null(treatment_var)){
    
    
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
    
    
    nlevels_biomarker <- nlevels(data[, biomarker_var])
    levels_biomarker <- levels(data[, biomarker_var])
    nlevels_treatment <- nlevels(data[, treatment_var])
    levels_treatment <- levels(data[, treatment_var])
    
    
    if(is.null(colors)){
      ## Some default colors that work for max 4 treatment levels
      
      if(nlevels_treatment > 4){
        stop("There are no default colors available when number of treatment levels is higher than 4. Please, provide the colors.")
      }
      
      default_colors_per_treatment <- list(c("#E69F00", "#B22222"), c("#56B4E9", "#104E8B"), c("#84d44b", "#007800"), c("#836FFF", "#50449d"))
      
      colors <- unlist(lapply(1:nlevels_treatment, function(i){
        format_colors(levels = paste0(levels_treatment[i], ", ", levels_biomarker), palette = default_colors_per_treatment[[i]], allow_duplicated = FALSE)
      }))
      
      
    }else{
      
      colors <- unlist(lapply(1:nlevels_treatment, function(i){
        # i = 1
        out <- format_colors(levels = levels_biomarker, colors = colors[[i]], allow_duplicated = FALSE)
        names(out) <- paste0(levels_treatment[i], ", ", levels_biomarker)
        return(out)
      }))
      
    }
    
    
    strat1_var <- treatment_var
    
  }else{
    
    covariate_var <- biomarker_var
    strat1_var <- NULL
    
    
  }
  
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  
  ggpl <- wrapper_core_KM_plot_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var,
    strat1_var = strat1_var, strat2_var = strat2_var, 
    colors = colors, 
    variable_names = variable_names, 
    title = title, subtitle_label_both = subtitle_label_both, tag_label_both = tag_label_both, 
    legend_title_colors = legend_title_colors, legend_position = legend_position, legend_justification = legend_justification,
    break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, 
    title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
    background_grid_major = background_grid_major,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  
  ggpl
  
  
}












#' KM plot - Treatment effect per biomarker subgroup
#' 
#' @param data Data frame.
#' @param colors If biomarker_var = NULL, then a vector. Otherwise, a list of length equal to tratment levels.
#' @export
wrapper_KM_plot_treatment <- function(data, tte_var, censor_var, biomarker_var = NULL, treatment_var, 
  strat2_var = NULL,
  colors = NULL, 
  variable_names = NULL, 
  title = NULL, subtitle_label_both = TRUE, tag_label_both = TRUE, 
  legend_title_colors = NULL, legend_position = c(0.03, 0.03), legend_justification = c(0, 0),
  break_time_by = NULL, max_tte = NULL, risk_table = TRUE, conf_int = FALSE, 
  title_size = 12, label_size = 3, rel_heights = c(5, 1), 
  background_grid_major = "none",
  strat_scales = "fixed", strat1_nrow = 1, strat1_ncol = NULL, strat2_nrow = NULL, strat2_ncol = 1){
  
  
  ### TODO Do not display biomarker in the legend
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(is.null(title)){
    title <- variable_names[biomarker_var]
  }
  
  
  
  if(!is.null(biomarker_var)){
    
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
    
    
    nlevels_biomarker <- nlevels(data[, biomarker_var])
    levels_biomarker <- levels(data[, biomarker_var])
    nlevels_treatment <- nlevels(data[, treatment_var])
    levels_treatment <- levels(data[, treatment_var])
    
    
    if(is.null(colors)){
      ## Some default colors that work for max 4 treatment levels
      
      if(nlevels_treatment > 4){
        stop("There are no default colors available when number of treatment levels is higher than 4. Please, provide the colors.")
      }
      
      default_colors_per_treatment <- list(c("#E69F00", "#B22222"), c("#56B4E9", "#104E8B"), c("#84d44b", "#007800"), c("#836FFF", "#50449d"))
      
      colors <- unlist(lapply(1:nlevels_treatment, function(i){
        format_colors(levels = paste0(levels_treatment[i], ", ", levels_biomarker), palette = default_colors_per_treatment[[i]], allow_duplicated = FALSE)
      }))
      
      
    }else{
      
      colors <- unlist(lapply(1:nlevels_treatment, function(i){
        # i = 1
        out <- format_colors(levels = levels_biomarker, colors = colors[[i]], allow_duplicated = FALSE)
        names(out) <- paste0(levels_treatment[i], ", ", levels_biomarker)
        return(out)
      }))
      
    }
    
    
    strat1_var <- biomarker_var
    
  }else{
    
    covariate_var <- treatment_var
    strat1_var <- NULL
    
  }
  
  
  # -------------------------------------------------------------------------
  # Plot
  # -------------------------------------------------------------------------
  
  
  ggpl <- wrapper_core_KM_plot_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var,
    strat1_var = strat1_var, strat2_var = strat2_var, 
    colors = colors, 
    variable_names = variable_names, 
    title = title, subtitle_label_both = subtitle_label_both, tag_label_both = tag_label_both, 
    legend_title_colors = legend_title_colors, legend_position = legend_position, legend_justification = legend_justification,
    break_time_by = break_time_by, max_tte = max_tte, risk_table = risk_table, conf_int = conf_int, 
    title_size = title_size, label_size = label_size, rel_heights = rel_heights, 
    background_grid_major = background_grid_major,
    strat_scales = strat_scales, strat1_nrow = strat1_nrow, strat1_ncol = strat1_ncol, strat2_nrow = strat2_nrow, strat2_ncol = strat2_ncol)
  
  
  
  ggpl
  
  
}




























