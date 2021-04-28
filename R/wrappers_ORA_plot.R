



# geneset_var = "GenesetID"; observed_var = "Observed"; adjp_var = "adj.P.Val";  
# title = ""; title_size = 10; axis_text_y_size = 10; axis_text_y_width = 80; color_point = 'darkslateblue'; size_range = c(2, 10)



#' Dot plot with ORA results for a single contrast 
#' 
#' @param x TopTable with selected ORA results obtained by running 'wrapper_ora_dispaly_significant'.
#' @export
wrapper_ora_dotplot_single <- function(x, geneset_var = "GenesetID", observed_var = "Observed", 
  adjp_var = "adj.P.Val", color_point_var = NULL,
  trim_limits = 0.01, color_point = 'darkslateblue',
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', 
  size_range = c(2, 10),
  title = "", title_size = 10, title_width = 100, axis_text_y_size = 8, axis_text_y_width = 80){
  
  
  stopifnot(length(geneset_var) == 1)
  
  ## Wrap the title so it can be nicely displayed in the plots
  if(title_width > 0){
    title <- stringr::str_wrap(title, width = title_width)
  }
  
  
  ## Wrap the gene set names so they can be nicely displayed in the plots
  x[, geneset_var] <- stringr::str_wrap(x[, geneset_var], width = axis_text_y_width)
  
  x[, geneset_var] <- factor(x[, geneset_var], levels = rev(x[, geneset_var]))
  
  
  ## TODO Display genes with geom_text in a separate plot 
  
  
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
      
      ggp <- ggplot(x, aes(x = .data[["log_adjp"]], y = .data[[geneset_var]], size = .data[["DE_in_set"]], color = .data[[color_point_var]])) +
        geom_point() +
        scale_size_area(name = "No. DE in set", max_size = size_range[2]) +
        scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish)
      
    }else{
      
      ggp <- ggplot(x, aes(x = .data[["log_adjp"]], y = .data[[geneset_var]], size = .data[["DE_in_set"]])) +
        geom_point(color = color_point) +
        scale_size_area(name = "No. DE in set", max_size = size_range[2]) 
      
    }
    
  }else{
    
    if(!is.null(color_point_var)){
      
      ggp <- ggplot(x, aes(x = .data[["log_adjp"]], y = .data[[geneset_var]], color = .data[[color_point_var]])) +
        geom_point(size = size_range[2]) +
        scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish)
      
    }else{
      
      ggp <- ggplot(x, aes(x = .data[["log_adjp"]], y = .data[[geneset_var]])) +
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
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      plot.title = element_text(size = title_size, hjust = 0.5),
      legend.position = "right") +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_y_discrete(position = "left")
  
  
  ggp
  
  
  
}








# x <- topTable_significant_ora
# 
# geneset_var = "GenesetID"
# 
# observed_prefix = "Observed"; adjp_prefix = "adj.P.Val";  sep = "_";
# title = ""; title_size = 10; title_width = ; axis_text_y_size = 10; axis_text_y_width = 80; colors_point = NULL; size_range = c(2, 10)
# point_alpha = 0.8





#' Dot plot with ORA results for multiple contrasts 
#' 
#' @param x TopTable with ORA results.
#' @export
wrapper_ora_dotplot_multiple <- function(x, geneset_var = "GenesetID", observed_prefix = "Observed", adjp_prefix = "adj.P.Val", sep = "_", directions = c("both", "up", "down"),
  title = "", title_size = 10, title_width = 100, axis_text_y_size = 8, axis_text_y_width = 80, 
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
    dplyr::left_join(data_observed, by = c(geneset_var, "contrasts_and_directions")) %>% 
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
  
  data <- data[stats::complete.cases(data$directions), , drop = FALSE]
  
  
  
  
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
  
  colors_point <- format_colors(contrasts, colors = colors_point)
  
  
  if(!is.null(observed_var)){
    
    ggp <- ggplot(data, aes(x = .data[["log_adjp"]], y = .data[[geneset_var]], color = .data[["contrasts"]], size = .data[["DE_in_set"]])) +
      geom_point(alpha = point_alpha) +
      scale_size_area(name = "No. DE in set", max_size = size_range[2]) 
    
  }else{
    
    ggp <- ggplot(data, aes(x = .data[["log_adjp"]], y = .data[[geneset_var]], color = .data[["contrasts"]])) +
      geom_point(alpha = point_alpha, size = size_range[2]) 
    
  }
  
  
  xintercept <- -log10(c(0.05, 0.1))
  xlab <- paste0("-log10(", adjp_var, ")")
  xlim <- c(0, max(c(data[, "log_adjp"], -log10(0.05)), na.rm = TRUE))
  
  
  ggp <- ggp +
    geom_vline(xintercept = xintercept, linetype = 2, color = "lightgrey") +
    ggtitle(title) +
    xlab(xlab) +
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      plot.title = element_text(size = title_size, hjust = 0.5),
      strip.background = element_rect(colour = "white", fill = "white")) +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_color_manual(name = "Contrast", values = colors_point) +
    facet_wrap(~ directions, nrow = 1)
  
  
  ggp
  
  
  
  
  
}
























