



# geneset_var = "Geneset"; observed_var = "Observed"; adjp_var = "adj.P.Val";  
# title = ""; title_size = 10; axis_text_y_size = 10; axis_text_y_width = 70; color_point = 'darkslateblue'; size_range = c(2, 10)



#' Dot plot with ORA results for a single contrast 
#' 
#' @param x TopTable with selected ORA results obtained by running 'wrapper_dispaly_significant_ora'.
#' @export
wrapper_plot_ORA_dotplot_single <- function(x, geneset_var = "Geneset", observed_var = "Observed", adjp_var = "adj.P.Val",  genes_prefix = NA, color_point_var = NULL,
  trim_limits = 0.01, color_point = 'darkslateblue',
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', 
  size_range = c(2, 10),
  title = "", title_size = 10, title_width = 200, axis_text_y_size = 8, axis_text_y_width = 150){
  
  
  stopifnot(length(geneset_var) == 1)
  
  ## Wrap the title so it can be nicely displayed in the plots
  if(title_width > 0){
    title <- stringr::str_wrap(title, width = title_width)
  }
  
  
  
  ## Add gene lists to the gene set names 
  ## TODO Display genes with geom_text in a separate plot 
  
  if(!is.na(genes_prefix)){
    
    
    genes_vars <- setdiff(grep(paste0("^", genes_prefix), colnames(x), value = TRUE), "Geneset")
    if(genes_vars > 1){
      for(i in seq_along(genes_vars)){
        x[, genes_vars[i]] <- ifelse(x[, genes_vars[i]] == "", "", paste0(gsub(paste0(genes_prefix, ":"), "", genes_vars[i]), ": ", x[, genes_vars[i]]))
        
      }
    }
    
    gene_list_suffix <- paste0(" [", unlist(apply(x[, genes_vars, drop = FALSE], 1, paste0, collapse = " ")), "]") 
    gene_list_suffix <- stringr::str_wrap(gene_list_suffix, width = axis_text_y_width)
    
    x[, geneset_var] <- paste0(x[, geneset_var], "\n", gene_list_suffix)
    
    
  }
  
  
  x[, geneset_var] <- factor(x[, geneset_var], levels = rev(x[, geneset_var]))
  
  
  ## To avoid p-values equal to zero
  min_non_zero_adjp <- min(x[x[, adjp_var] > 0, adjp_var], na.rm = TRUE)
  x[x[, adjp_var] == 0, adjp_var] <- min(1e-15, min_non_zero_adjp)
  
  
  x$log_adjp <- -log10(x[, adjp_var])
  
  ## Derive the number of DE genes that overlap with the gene set
  
  if(!is.null(observed_var)){
    x$DE_in_set <- as.numeric(x[, observed_var])
    x$DE_in_set[x$DE_in_set == 0] <- NA
  }
  
  
  if(!is.null(color_point_var)){
    
    ### Trim limits
    
    if(is.null(trim_limits)){
      max_abs_value <- ceiling(max(abs(range(x[, color_point_var], na.rm = TRUE))))
    }else if(trim_limits >= 1){
      max_abs_value <- trim_limits
    }else{
      ### Use quantiles 
      max_abs_value <- ceiling(max(abs(quantile(x[, color_point_var], probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
    }
    
    
    limits <- c(-max_abs_value, max_abs_value)
    
  }
  
  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  if(!is.null(observed_var)){
    
    if(!is.null(color_point_var)){
      
      ggp <- ggplot(x, aes_string(x = "log_adjp", y = geneset_var, size = "DE_in_set", color = color_point_var)) +
        geom_point() +
        scale_size(name = "No. DE in set", range = size_range) +
        scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish)
      
    }else{
      
      ggp <- ggplot(x, aes_string(x = "log_adjp", y = geneset_var, size = "DE_in_set")) +
        geom_point(color = color_point) +
        scale_size(name = "No. DE in set", range = size_range) 
      
    }
    
  }else{
    
    if(!is.null(color_point_var)){
      
      ggp <- ggplot(x, aes_string(x = "log_adjp", y = geneset_var, color = color_point_var)) +
        geom_point(size = size_range[2]) +
        scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish)
      
    }else{
      
      ggp <- ggplot(x, aes_string(x = "log_adjp", y = geneset_var)) +
        geom_point(color = color_point, size = size_range[2])
      
    }
    
  }
  
  
  xintercept <- -log10(c(0.05, 0.1))
  xlab <- paste0("-log10(", adjp_var, ")")
  xlim <- c(0, max(c(x[, "log_adjp"], -log10(0.05)), na.rm = TRUE))
  
  
  ggp <- ggp +
    geom_vline(xintercept = xintercept, linetype = 2, color = "lightgrey") +
    ggtitle(title) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      plot.title = element_text(size = title_size, hjust = 1),
      legend.position = "right") +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_y_discrete(position = "left")
  
  
  ggp
  
  
  
}








# x <- topTable_significant_ora
# 
# geneset_var = "Geneset"
# 
# observed_prefix = "Observed"; adjp_prefix = "adj.P.Val";  sep = "_";
# title = ""; title_size = 10; title_width = 70; axis_text_y_size = 10; axis_text_y_width = 70; colors_point = NULL; size_range = c(2, 10)
# point_alpha = 0.8





#' Dot plot with ORA results for a single contrast 
#' 
#' @param x TopTable with ORA results.
#' @export
wrapper_plot_ORA_dotplot_multiple <- function(x, geneset_var = "Geneset", observed_prefix = "Observed", adjp_prefix = "adj.P.Val", sep = "_", directions = c("both", "up", "down"),
  title = "", title_size = 10, title_width = 140, axis_text_y_size = 8, axis_text_y_width = 70, 
  colors_point = NULL, size_range = c(2, 10), point_alpha = 0.8){
  
  
  stopifnot(length(geneset_var) == 1)
  stopifnot(all(directions %in% c("up", "down", "both")))
  
  directions_labels <- c("up" = "Up-regulation", "down" = "Down-regulation", "both" = "Up-Down-regulation")
  directions_labels <- directions_labels[directions]
  
  
  ## Wrap the title so it can be nicely displayed in the plots
  title <- stringr::str_wrap(title, width = title_width)
  
  
  ## Wrap the gene set names so they can be nicely displayed in the plots
  x[, geneset_var] <- stringr::str_wrap(x[, geneset_var], width = axis_text_y_width)
  
  x[, geneset_var] <- factor(x[, geneset_var], levels = rev(unique(x[, geneset_var])))
  
  
  # -----------------------------------------------------------------
  # Retrieve contrasts and data
  # -----------------------------------------------------------------
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts_and_directions <- gsub(paste0("^", adjp_prefix, sep), "", grep(paste0("^", adjp_prefix, sep), colnames(x), value = TRUE))
  
  
  
  ## We add '$' because we want to match expression at the end of the string
  contrasts <- gsub(paste0(sep, "up$"), "", contrasts_and_directions)
  contrasts <- gsub(paste0(sep, "down$"), "", contrasts)
  contrasts <- gsub(paste0(sep, "both$"), "", contrasts)
  contrasts <- unique(contrasts)
  
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix)
  
  data_observed <- wrapper_extract_from_topTable(x, extract_prefix = observed_prefix)
  
  
  
  data_adjp <- pivot_longer(data.frame(x[, geneset_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE), 
    cols = contrasts_and_directions, names_to = "contrasts_and_directions", values_to = adjp_prefix)
  
  data_observed <- pivot_longer(data.frame(x[, geneset_var, drop = FALSE], data_observed, stringsAsFactors = FALSE), 
    cols = contrasts_and_directions, names_to = "contrasts_and_directions", values_to = observed_prefix)
  
  
  
  data <- data_adjp %>% 
    left_join(data_observed, by = c(geneset_var, "contrasts_and_directions")) %>% 
    as.data.frame()
  
  
  ### Retrive contrasts
  data$contrasts <- gsub(paste0(sep, "up$"), "", data$contrasts_and_directions)
  data$contrasts <- gsub(paste0(sep, "down$"), "", data$contrasts)
  data$contrasts <- gsub(paste0(sep, "both$"), "", data$contrasts)
  
  data$contrasts <- factor(data$contrasts, levels = contrasts)
  
  ### Retrive directions
  
  data$directions <- sapply(seq_len(nrow(data)), function(i){
    gsub(paste0("^", data$contrasts[i], sep), "", data$contrasts_and_directions[i])
  })
  
  
  data$directions <- factor(data$directions, levels = directions, labels = directions_labels)
  
  data <- data[complete.cases(data$directions), , drop = FALSE]
  
  
  
  
  adjp_var <- adjp_prefix
  observed_var <- observed_prefix
  
  
  ## To avoid p-values equal to zero
  min_non_zero_adjp <- min(data[data[, adjp_var] > 0, adjp_var], na.rm = TRUE)
  data[data[, adjp_var] == 0, adjp_var] <- min(1e-15, min_non_zero_adjp)
  
  
  data$log_adjp <- -log10(data[, adjp_var])
  
  
  ## Derive the number of DE genes that overlap with the gene set
  
  if(!is.null(observed_var)){
    data$DE_in_set <- as.numeric(data[, observed_var])
    data$DE_in_set[data$DE_in_set == 0] <- NA
  }
  
  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  colors_point <- format_colors(levels = contrasts, colors = colors_point)
  
  
  if(!is.null(observed_var)){
    
    ggp <- ggplot(data, aes_string(x = "log_adjp", y = geneset_var, color = "contrasts", size = "DE_in_set")) +
      geom_point(alpha = point_alpha) +
      scale_size(name = "No. DE in set", range = size_range) 
    
  }else{
    
    ggp <- ggplot(data, aes_string(x = "log_adjp", y = geneset_var, color = "contrasts")) +
      geom_point(alpha = point_alpha, size = size_range[2]) 
    
  }
  
  
  xintercept <- -log10(c(0.05, 0.1))
  xlab <- paste0("-log10(", adjp_var, ")")
  xlim <- c(0, max(c(data[, "log_adjp"], -log10(0.05)), na.rm = TRUE))
  
  
  ggp <- ggp +
    geom_vline(xintercept = xintercept, linetype = 2, color = "lightgrey") +
    ggtitle(title) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      plot.title = element_text(size = title_size, hjust = 1),
      strip.background = element_rect(colour = "white", fill = "white")) +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_color_manual(name = "Contrast", values = colors_point) +
    facet_wrap(~ directions, nrow = 1)
  
  
  ggp
  
  
  
  
  
}
























