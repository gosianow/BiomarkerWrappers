

#' Bar plot with NES from GSEA
#' 
#' So far it works only for one contrast.
#' 
#' @param x TopTable with results from GSEA
#' @export
wrapper_NES_barplot <- function(x, geneset_var = "GenesetID", 
  enrichment_score_prefix = "NES", adjp_prefix = "adj.P.Val",
  pval = 0.05, sep = "_", title = "", legend_position = "right",
  axis_text_y_size = NULL, axis_text_y_width = 80, title_size = NULL){
  
  
  stopifnot(length(geneset_var) == 1)
  
  
  data_NES <- wrapper_extract_from_topTable(x, extract_prefix = enrichment_score_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  contrasts <- colnames(data_NES)
  contrasts
  
  if(length(contrasts) > 1){
    stop("So far this plot works only for one contrast!")
  }
  
  
  data_NES <- pivot_longer(data.frame(x[, geneset_var, drop = FALSE], data_NES, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = enrichment_score_prefix)
  
  data_adjp <- pivot_longer(data.frame(x[, geneset_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = adjp_prefix)
  
  

  data <- data_NES %>% 
    dplyr::left_join(data_adjp, by = c(geneset_var, "contrast")) %>% 
    data.frame()
  
  
  data$significance <- factor(ifelse(data[, adjp_prefix] <= pval, paste0("<=", pval), paste0(">", pval)), levels = paste0(c("<=", ">"), pval))
  
  values_color <- c("darkorchid4", "gold")
  names(values_color) <- levels(data$significance)
  
  
  data$contrast <- factor(data$contrast, levels = contrasts)
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  rownames_wrap <- stringr::str_wrap(data[, geneset_var], width = axis_text_y_width)
  
  rownames_wrap <- limma::strsplit2(rownames_wrap, split = "\\\n")[, 1]
  
  data[, geneset_var] <- factor(rownames_wrap, levels = rev(unique(rownames_wrap)))
  

  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  
  
  ggp <- ggplot(data, aes(x = .data[[enrichment_score_prefix]], y = .data[[geneset_var]], fill = .data[["significance"]])) +
    geom_bar(stat = "identity") +
    ggtitle(title) +
    theme(plot.title = element_text(size = title_size),
      axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      legend.position = legend_position) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_fill_manual(name = adjp_prefix, values = values_color, drop = FALSE) 
  
  
  
  return(ggp)
  
  
  
}





#' Bar plot with NES from GSEA
#' 
#' So far it works only for one contrast.
#' 
#' @param x TopTable with results from GSEA
#' @export
wrapper_NES_barplot2 <- function(x, geneset_var = "GenesetID", 
  enrichment_score_prefix = "NES", adjp_prefix = "adj.P.Val",
  pval = 0.05, sep = "_", title = "", legend_position = "right",
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', 
  trim_values = 3, trim_prop = NULL, trim_range = NULL, ceiling = FALSE,
  axis_text_y_size = NULL, axis_text_y_width = 80, title_size = NULL){
  
  
  stopifnot(length(geneset_var) == 1)
  
  
  data_NES <- wrapper_extract_from_topTable(x, extract_prefix = enrichment_score_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  contrasts <- colnames(data_NES)
  contrasts
  
  if(length(contrasts) > 1){
    stop("So far this plot works only for one contrast!")
  }
  
  
  data_NES <- pivot_longer(data.frame(x[, geneset_var, drop = FALSE], data_NES, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = enrichment_score_prefix)
  
  data_adjp <- pivot_longer(data.frame(x[, geneset_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = adjp_prefix)
  
  
  data <- data_NES %>% 
    dplyr::left_join(data_adjp, by = c(geneset_var, "contrast")) %>% 
    data.frame()
  
  data$statistic <- -log10(data[, adjp_prefix]) * sign(data[, enrichment_score_prefix])
  
  
  data$contrast <- factor(data$contrast, levels = contrasts)
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  rownames_wrap <- stringr::str_wrap(data[, geneset_var], width = axis_text_y_width)
  
  rownames_wrap <- limma::strsplit2(rownames_wrap, split = "\\\n")[, 1]
  
  data[, geneset_var] <- factor(rownames_wrap, levels = rev(unique(rownames_wrap)))
  
  
  if(is.null(trim_values)){
    trim_values <- compute_trim_values(x = data[, enrichment_score_prefix], centered = TRUE, trim_prop = trim_prop, trim_range = trim_range, ceiling = ceiling)
  }else{
    max_abs_value <- max(abs(trim_values))
    trim_values <- c(-max_abs_value, max_abs_value)
  }
  
  
  limits <- trim_values
  
  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  
  xintercept <- -log10(pval) * unique(sign(data[, "statistic"]))
  xlab <- paste0("-log10(", adjp_prefix, ") * sign(", enrichment_score_prefix, ")")
  
  
  
  ggp <- ggplot(data, aes(x = .data[["statistic"]], y = .data[[geneset_var]], fill = .data[[enrichment_score_prefix]])) +
    geom_bar(stat = "identity") +
    geom_vline(xintercept = xintercept, linetype = 2, color = "lightgrey") +
    ggtitle(title) +
    xlab(xlab) +
    theme(plot.title = element_text(size = title_size),
      axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      legend.position = legend_position) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_fill_gradient2(name = enrichment_score_prefix, low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish)
  
  
  
  return(ggp)
  
  
  
}



































