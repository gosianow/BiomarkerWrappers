



# data <- data
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
wrapper_core_bar_plot <- function(data, x_var, y_var, colors = NULL, variable_names = NULL, title = NULL, subtitle = NULL, tag = NULL, skip_levels = NULL, show_total_counts = TRUE, show_proportions = TRUE, show_total_proportions = FALSE, title.size = 12, ylim = c(0, 100), axis.text.x.angle = 0, axis.text.x.vjust = 0, axis.text.x.hjust = 0.5, geom_text_size = 3, geom_text_vjust = 0.5, background_grid_major = "none"){
  
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.factor(data[, y_var]))
  
  
  colors <- format_colors(levels = levels(data[, y_var]), colors = colors)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  ### Calculate counts and proportions
  
  tbl <- table(data[, x_var], data[, y_var])
  
  prop <- prop.table(tbl, margin = 1) * 100
  
  
  countdf <- as.data.frame.matrix(tbl)
  countdf <- data.frame(Subgroup = rownames(countdf), countdf, check.names = FALSE, row.names = NULL)
  
  propdf <- as.data.frame.matrix(prop)
  propdf <- data.frame(Subgroup = rownames(propdf), propdf, check.names = FALSE, row.names = NULL)
  
  
  ### Prepare data for ggplot
  
  ggdata <- cbind(gather(propdf, key = "Observation", value = "Proportion", -Subgroup), 
    gather(countdf, key = "Observation", value = "Count", -Subgroup)[, "Count", drop = FALSE])
  
  
  ### Remove data with zero counts so the ggplot warning with missing values is not displayed
  ggdata <- ggdata[ggdata$Count > 0, ]
  
  
  ### Define lables 
  ggdata$Label <- paste0(ggdata$Count, " (", round(ggdata$Proportion, 1), "%)")
  
  
  
  ### Calculate marginal counts per Subgroup
  if(show_total_counts){
    ggdata$Subgroup <- factor(ggdata$Subgroup, levels = levels(data[, x_var]), labels = paste0(levels(data[, x_var]), "\n", rowSums(tbl)))
  }else{
    ggdata$Subgroup <- factor(ggdata$Subgroup, levels = levels(data[, x_var]))
  }
  
  
  ggdata$Observation <- factor(ggdata$Observation, levels = levels(data[, y_var]))
  
  
  ### Skip unwanted levels
  
  keep_levels <- levels(ggdata$Observation)
  
  if(!is.null(skip_levels)){
    
    keep_levels <- setdiff(levels(ggdata$Observation), skip_levels)
    
    stopifnot(length(keep_levels) > 0)
    
    ggdata <- ggdata[ggdata$Observation %in% keep_levels, , drop = FALSE]
    
    ggdata$Observation <- factor(ggdata$Observation, keep_levels)
    
  }
  
  
  xlab <- variable_names[x_var]
  ylab <- paste0(variable_names[y_var], "\nProportion (%)")
  
  
  ### Make the plot
  
  
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
  
  
  if(show_proportions){
    
    ggpl <- ggpl +
      geom_text(data = ggdata, aes(x = Subgroup, y = Proportion, group = Observation, label = Label), 
        position = position_stack(vjust = geom_text_vjust), size = geom_text_size, na.rm = TRUE)
    
  }
  
  
  if(show_total_proportions){
    
    ggdata_total_prop <- data.frame(Subgroup = factor(rownames(prop), levels = rownames(prop), labels = levels(ggdata$Subgroup)), 
      Total_Proportion = rowSums(prop[, keep_levels]))
    
    ggdata_total_prop$Label <- paste0(round(ggdata_total_prop$Total_Proportion, 1), "%")
    
    
    ggpl <- ggpl +
      geom_text(data = ggdata_total_prop, aes(x = Subgroup, y = Total_Proportion + 1, label = Label), size = geom_text_size, nudge_y = 1)
    
    
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

















