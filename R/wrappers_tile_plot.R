






#' Tile plot
#' 
#' @param data Data frame.
#' @param y_vars Names of variables to be plotted as tiles. Those variables must be factors.
#' @param colors Named list with colors for 'y_vars'.
#' @param variable_names Named vector with nicer variable names.
#' @param skip_NAs Logical vector. Whether to skip NAs.
#' @export
wrapper_tile_plot1_core <- function(data, y_vars, colors = NULL, variable_names = NULL, skip_NAs = FALSE, rev = FALSE, order = TRUE, return_plotlist = FALSE, nrow_legend = 3){
  
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "factor"))
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  colors <- lapply(seq_along(y_vars), function(i){
    # i = 1
    y_var <- y_vars[i]
    format_colors(levels(data[, y_var]), colors = colors[[y_var]])
  })
  names(colors) <- y_vars
  
  
  
  
  if(length(rev) == 1){
    rev <- rep(rev, length(y_vars))
  }
  
  stopifnot(length(rev) == length(y_vars))
  
  for(i in seq_along(rev)){
    # i = 1
    
    if(!rev[i]){
      
      data[, y_vars[i]] <- factor(data[, y_vars[i]], levels = rev(levels(data[, y_vars[i]])))
      
    }
    
  }
  
  
  
  if(any(skip_NAs)){
    
    if(length(skip_NAs) == 1){
      skip_NAs <- rep(skip_NAs, length(y_vars))
    }
    
    data <- data[complete.cases(data[, y_vars[skip_NAs], drop = FALSE]), , drop = FALSE]
  }
  
  
  if(order){
    data <- data[order2(data[, y_vars, drop = FALSE], decreasing = FALSE), , drop = FALSE]
  }
  
  
  x_var <- "dummy_x_var"
  
  data[, x_var] <- factor(paste0("Sample", seq_len(nrow(data))), levels = paste0("Sample", seq_len(nrow(data))))
  
  
  plotlist <- lapply(seq_along(y_vars), function(i){
    # i = 1
    
    y_var <- y_vars[i]
    
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
        # legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "line"),
        plot.margin = margin(t = 1, r = 7, b = 1, l = 7, unit = "pt")) +
      scale_fill_manual(values = colors[[y_var]], na.value = "gray95") +
      scale_y_discrete(expand = c(0,0)) +
      guides(fill = guide_legend(nrow = nrow_legend, reverse = FALSE))
    
    
  })
  
  
  if(return_plotlist){
    plotlist
  }else{
    plot_grid(plotlist = plotlist, ncol = 1, align = "v", axis = "lr")  
  }
  
  
  
  
}








#' @rdname wrapper_tile_plot1_core
#' @export
wrapper_tile_plot2_core <- function(data, y_vars, colors = NULL, variable_names = NULL, skip_NAs = FALSE, rev = FALSE, return_plotlist = FALSE, nrow_legend = 3){
  
  
  stopifnot(length(y_vars) >= 1)
  stopifnot(all(sapply(data[, y_vars], class) == "factor"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  if(any(rev == TRUE)){
    
    if(length(rev) == 1){
      rev <- rep(rev, length(y_vars))
    }
    
    stopifnot(length(rev) == length(y_vars))
    
    for(i in seq_along(rev)){
      # i = 1
      
      if(rev[i]){
        
        data[, y_vars[i]] <- factor(data[, y_vars[i]], levels = rev(levels(data[, y_vars[i]])))
        
      }
      
    }
    
  }
  
  
  
  levels_original <- lapply(data[, y_vars, drop = FALSE], levels)
  names(levels_original) <- y_vars
  
  
  if(any(skip_NAs)){
    
    if(length(skip_NAs) == 1){
      skip_NAs <- rep(skip_NAs, length(y_vars))
    }
    
    data <- data[complete.cases(data[, y_vars[skip_NAs], drop = FALSE]), , drop = FALSE]
    
  }
  
  if(any(skip_NAs == FALSE)){
    
    ### To calculate proportions of interaction factors we have to add NA as a factor level
    
    for(i in seq_along(y_vars)){
      # i = 2
      
      y_var <- y_vars[i]
      
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
    
    ### Information about marginal proportions
    
    tbl_y_var <- table(data[, y_vars[i]])
    prop_y_var <- tbl_y_var / nrow(data) * 100
    
    
    
    out <- data.frame(y_var = y_vars[i], 
      names_y_var = variable_names[y_vars[i]],
      value_interaction = factor(names(tbl), levels = levels(value_interaction)), 
      levels_y_var = factor(levels_y_var, levels = levels(data[, y_vars[i]])), 
      levels_with_prop_y_var = factor(levels_y_var, levels = levels(data[, y_vars[i]]), 
        labels = paste0(levels(data[, y_vars[i]]), " ", format_props(prop_y_var, digits = 0))), 
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
    
    colors_tmp <- format_colors(levels_original[[y_var]], colors = colors[[y_var]])
    colors_tmp <- c(colors_tmp, "NA" = "gray95")
    
    ### Use names with proportions
    
    prop_tbl <- unique(data_tile_size[, c("levels_y_var", "levels_with_prop_y_var")])
    
    mm <- match(names(colors_tmp), prop_tbl$levels_y_var)
    names(colors_tmp) <- prop_tbl$levels_with_prop_y_var[mm]
    
    
    
    ggplot(data_tile_size) +
      geom_col(aes(x = names_y_var, y = prop, group = value_interaction, fill = levels_with_prop_y_var)) +
      coord_flip() +
      theme(axis.text.x = element_blank(), 
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        legend.position = "right",
        legend.title = element_blank(), 
        # legend.text = element_text(size = 10),
        legend.key.size = unit(0.5, "line"),
        plot.margin = margin(t = 2, r = 7, b = 2, l = 7, unit = "pt")) +
      scale_fill_manual(values = colors_tmp) +
      scale_y_discrete(expand = c(0,0)) +
      scale_x_discrete(expand = c(0,0)) +
      guides(fill = guide_legend(nrow = nrow_legend, reverse = TRUE))
    
    
  })
  
  if(return_plotlist){
    plotlist
  }else{
    plot_grid(plotlist = plotlist, ncol = 1, align = "v", axis = "lr") 
  }
  
  
}














