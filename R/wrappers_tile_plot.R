






#' Tile plot
#' 
#' @param data Data frame.
#' @param x_var Name of a variable that defined objects unique IDs.
#' @param y_vars Names of variables to be plotted as tiles. Those variables must be factors.
#' @param colors Named list with colors for 'y_vars'.
#' @param variable_names Named vector with nicer variable names.
#' @param skip_NAs Logical. Whether to skip NAs.
#' @export
wrapper_tile_plot1_core <- function(data, x_var, y_vars, colors = NULL, variable_names = NULL, skip_NAs = FALSE){
  
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "factor"))
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(skip_NAs){
    data <- data[complete.cases(data), , drop = FALSE]
  }
  
  data <- data[order2(data[, y_vars, drop = FALSE], decreasing = TRUE), , drop = FALSE]
  
  data[, x_var] <- factor(data[, x_var], levels = data[, x_var])
  
  
  plotlist <- lapply(seq_along(y_vars), function(i){
    # i = 1
    
    y_var <- y_vars[i]
    
    colors_tmp <- format_colors(levels = levels(data[, y_var]), colors = colors[[y_var]])
    
    data$dummy <- variable_names[y_var]
    
    
    # theme_cowplot()$plot.margin
    
    ggplot(data) +
      geom_tile(aes(x = .data[[x_var]], y = .data[["dummy"]], fill = .data[[y_var]])) +
      theme(axis.text.x = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "line"),
        plot.margin = margin(t = 1, r = 7, b = 1, l = 7, unit = "pt")) +
      scale_fill_manual(values = colors_tmp, na.value = "gray95") +
      scale_y_discrete(expand = c(0,0)) +
      guides(fill = guide_legend(nrow = 4, reverse = TRUE))
    
    
  })
  
  
  plot_grid(plotlist = plotlist, ncol = 1, align = "v", axis = "lr") 
  
  
}






#' @rdname wrapper_tile_plot1_core
#' @export
wrapper_tile_plot2_core <- function(data, x_var, y_vars, colors = NULL, variable_names = NULL, skip_NAs = FALSE){
  
  
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "factor"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  if(skip_NAs){
    data <- data[complete.cases(data), , drop = FALSE]
  }else{
    
    ### To calculate proportions of interaction factors we have to add NA as a factor level
    
    levels_original <- list()
    
    for(i in seq_along(y_vars)){
      # i = 2
      
      y_var <- y_vars[i]
      levels_original[[y_var]] <- levels(data[, y_var])
        
      if(any(is.na(data[, y_var]))){
        data[, y_var] <- factor(data[, y_var], exclude = NULL, levels = c(NA, levels(data[, y_var])), labels = c("NA", levels(data[, y_var])))
      }
      
    }
    
    
  }
  
  
  
  ### Compute tile size 
  
  data_tile_size_list <- lapply(seq_along(y_vars), function(i){
    # i = 2
    
    value_interaction <- interaction(data[, y_vars[1:i], drop = FALSE], sep = " /dummy-dummy/ ", lex.order = TRUE)
    
    tbl <- table(value_interaction)
    prop <- tbl / nrow(data) * 100
    
    
    levels_y_var <- limma::strsplit2(names(tbl), split = " /dummy-dummy/ ")[, i]
    
    
    out <- data.frame(y_var = y_vars[i], 
      names_y_var = variable_names[y_vars[i]],
      value_interaction = factor(names(tbl), levels = levels(value_interaction)), 
      levels_y_var = factor(levels_y_var, levels = levels(data[, y_vars[i]])), 
      prop = as.numeric(prop), 
      stringsAsFactors = FALSE, row.names = NULL)
    
    
    return(out)
    
  })
  
  
  plotlist <- lapply(seq_along(data_tile_size_list), function(i){
    # i = 1
    
    data_tile_size <- data_tile_size_list[[i]]
    
    y_var <- data_tile_size$y_var[1]
    names_y_var <- data_tile_size$names_y_var[1]
    

    ### Update colors 
    
    colors_tmp <- format_colors(levels = levels_original[[y_var]], colors = colors[[y_var]])
    colors_tmp <- c(colors_tmp, "NA" = "gray95")

    
    
    ggplot(data_tile_size) +
      geom_col(aes(x = names_y_var, y = prop, group = value_interaction, fill = levels_y_var)) +
      coord_flip() +
      theme(axis.text.x = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        legend.title = element_blank(), 
        legend.text = element_text(size = 8),
        legend.key.size = unit(0.5, "line"),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 7, unit = "pt")) +
      scale_fill_manual(values = colors_tmp) +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      guides(fill = guide_legend(nrow = 4, reverse = TRUE))
    
    
  })
  
  
  plot_grid(plotlist = plotlist, ncol = 1, align = "v", axis = "lr") 
  
  
  
}














