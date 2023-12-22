


#' Dot plot with logFC and p-values for multiple contrasts
#' 
#' @param x TopTable
#' @export
wrapper_logFC_dotplot <- function(x, gene_var = "Hgnc_Symbol", 
  lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val",  
  sep = "_", pval = 0.05, title = "", 
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', 
  trim_values = 3, trim_prop = NULL, trim_range = NULL, ceiling = FALSE, 
  radius_range = c(3, 10), legend_position = "right", 
  axis_text_x_angle = 90, axis_text_x_vjust = 0.5, axis_text_x_hjust = 1, 
  axis_text_y_size = NULL, axis_text_y_width = 80, title_size = NULL){
  
  
  
  stopifnot(length(gene_var) == 1)
  
  
  data_lfc <- wrapper_extract_from_topTable(x, extract_prefix = lfc_prefix, sep = sep)
  
  data_pval <- wrapper_extract_from_topTable(x, extract_prefix = pval_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  contrasts <- colnames(data_lfc)
  contrasts
  
  
  data_lfc <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_lfc, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = lfc_prefix)
  
  data_pval <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_pval, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = pval_prefix)
  
  data_adjp <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = all_of(contrasts), names_to = "contrast", values_to = adjp_prefix)
  
  
  if(pval_prefix == adjp_prefix){
    data <- data_lfc %>% 
      dplyr::left_join(data_pval, by = c(gene_var, "contrast")) %>% 
      data.frame()
  }else{
    data <- data_lfc %>% 
      dplyr::left_join(data_pval, by = c(gene_var, "contrast")) %>% 
      dplyr::left_join(data_adjp, by = c(gene_var, "contrast")) %>% 
      data.frame()
  }
  
  
  data$significance <- factor(ifelse(data[, adjp_prefix] <= pval, paste0("<=", pval), paste0(">", pval)), levels = paste0(c("<=", ">"), pval))
  
  values_shape <- c(4, 32)
  names(values_shape) <- levels(data$significance)
  
  
  data$contrast <- factor(data$contrast, levels = contrasts)
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  rownames_wrap <- stringr::str_wrap(data[, gene_var], width = axis_text_y_width)
  
  rownames_wrap <- limma::strsplit2(rownames_wrap, split = "\\\n")[, 1]
  
  data[, gene_var] <- factor(rownames_wrap, levels = rev(unique(rownames_wrap)))
  
  
  ### Use -log10(p-value) to define the size of the bubbles 
  
  data$log.P.Val <- -log10(data[, pval_prefix])
  
  
  if(lfc_prefix %in% c("logFC", "NES", "statistic", "sign.P.Val")){
    radius_labels = rev(c(1, 0.1, 0.05, 0.01, 1e-04, 1e-06, 1e-08, 1e-10))
    radius_breaks = -log10(radius_labels)
    radius_limits = range(radius_breaks)
    radius_labels <- formatC(radius_labels, format = "g", digits = 1)
  }else{
    radius_labels = rev(c(1, 0.1, 0.05, 0.01, 0.001))
    radius_breaks = -log10(radius_labels)
    radius_limits = range(radius_breaks)
    radius_labels <- formatC(radius_labels, format = "f", drop0trailing = TRUE, digits = 10)
  }
  
  ### Squish the limits because oob = scales::squish is not possible for scale_radius
  data$log.P.Val[data$log.P.Val > max(radius_limits)] <- max(radius_limits)
 
   
  if(is.null(trim_values)){
    trim_values <- compute_trim_values(x = data[, lfc_prefix], centered = TRUE, trim_prop = trim_prop, trim_range = trim_range, ceiling = ceiling)
  }else{
    max_abs_value <- max(abs(trim_values))
    trim_values <- c(-max_abs_value, max_abs_value)
  }
  
  
  limits <- trim_values
  
  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  # This is great. squish in this context converts clamps all values to be within the min and max of the limits argument. i.e., if value < min(limits) then value = min(limits) else if value > max(limits) then value = max(limits).
  # scale_colour_gradient2(limits = c(-1.5, 1.5), oob = scales::squish)
  
  
  ggp <- ggplot(data, aes(x = .data[["contrast"]], y = .data[[gene_var]], size = .data[["log.P.Val"]], color = .data[[lfc_prefix]])) +
    geom_point() +
    geom_point(aes(size = .data[["log.P.Val"]] - 0.5, shape = .data[["significance"]]), color = "black", show.legend = TRUE) +
    ggtitle(title) +
    theme(plot.title = element_text(size = title_size),
      axis.line = element_blank(), 
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      axis.text.y = element_text(size = axis_text_y_size),
      legend.position = legend_position) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_shape_manual(name = adjp_prefix, values = values_shape, drop = FALSE) +
    scale_colour_gradient2(name = lfc_prefix, low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish) +
    scale_radius(name = pval_prefix, range = radius_range, breaks = radius_breaks, labels = radius_labels, limits = radius_limits) 
  
  
  ggp
  
  
  
}






#' Heatmap with logFC for multiple contrasts
#' 
#' @param x TopTable
#' @export
wrapper_logFC_heatmap <- function(x, gene_var = "Hgnc_Symbol", 
  lfc_prefix = "logFC", adjp_prefix = "adj.P.Val", sep = "_", 
  title = "",
  palette = c('#42399B', "white", '#D70131'), rev = FALSE,
  trim_values = 3, trim_prop = NULL, trim_range = NULL, ceiling = FALSE,
  font_size = 10,
  rect_gp = grid::gpar(col = "grey"), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_dend_reorder = FALSE, column_dend_reorder = FALSE,
  show_row_names = TRUE, show_column_names = TRUE, 
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9), 
  row_names_max_width = unit(20, "cm"), column_names_max_height = unit(20, "cm"),
  row_names_width = 80, draw = TRUE, return = "ht", ...){
  
  
  stopifnot(length(gene_var) == 1)
  
  name <- lfc_prefix 
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  rownames_wrap <- stringr::str_wrap(x[, gene_var], width = row_names_width)
  
  rownames_wrap <- limma::strsplit2(rownames_wrap, split = "\\\n")[, 1]
  
  rownames(x) <- rownames_wrap
  
  # ---------------------------------------------------------------------------
  # Prepare matrix
  # ---------------------------------------------------------------------------
  
  
  data_lfc <- wrapper_extract_from_topTable(x, extract_prefix = lfc_prefix, sep = sep)
  
  matrix <- as.matrix(data_lfc)
  
  
  # ---------------------------------------------------------------------------
  # Prepare data_significance
  # ---------------------------------------------------------------------------
  
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  data_significance <- data.frame(apply(data_adjp, 2, function(x){
    
    x[is.na(x)] <- 1
    
    pval_asterisk <- ifelse(x <= 0.0001, "****", ifelse(x <= 0.001, "***", ifelse(x <= 0.01, "**", ifelse(x <= 0.05, "*", ifelse(x <= 0.1, ".", "")))))
    
    pval_asterisk
    
  }), stringsAsFactors = FALSE)
  
  
  
  # ---------------------------------------------------------------------------
  # Colors for the heatmap expression
  # ---------------------------------------------------------------------------
  
  
  colors_matrix <- format_colors_num(x = matrix, centered = TRUE, palette = palette, rev = rev, trim_values = trim_values, trim_prop = trim_prop, trim_range = trim_range, ceiling = ceiling)
  
  
  
  # ---------------------------------------------------------------------------
  # cell_fun
  # ---------------------------------------------------------------------------
  
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#customize-the-heatmap-body
  
  
  cell_fun <- function(j, i, x, y, width, height, fill){
    
    
    # grid::grid.text(data_significance[i, j], x, y = y, gp = grid::gpar(col = "black", fontsize = font_size))
    
    ### Center the displayed items
    
    if(data_significance[i, j] == "."){
      grid::grid.text(data_significance[i, j], x, y = y + (0.3 * unit(font_size, units = "points")), gp = grid::gpar(col = "black", fontsize = font_size))
      
    }else{
      grid::grid.text(data_significance[i, j], x, y = y - (0.2 * unit(font_size, units = "points")), gp = grid::gpar(col = "black", fontsize = font_size))
      
    }
    
    
  }
  
  
  
  # cell_fun <- function(j, i, x, y, width, height, fill){ 
  #   
  #   if(data_significance[i, j] == "***"){
  #     grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "black"))
  #   }
  #   
  #   if(data_significance[i, j] == "**"){
  #     grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "dimgrey"))
  #   }
  #   
  #   if(data_significance[i, j] == "*"){
  #     grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "darkgrey"))
  #   }
  #   
  #   if(data_significance[i, j] == "."){
  #     grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "transparent"))
  #   }
  #   
  # }
  
  
  # ---------------------------------------------------------------------------
  # Heatmap
  # ---------------------------------------------------------------------------
  
  
  ht <- ComplexHeatmap::Heatmap(matrix, 
    name = name,
    col = colors_matrix, 
    na_col = "grey90",
    rect_gp = rect_gp,
    heatmap_legend_param = list(color_bar = "continuous"), 
    
    cell_fun = cell_fun,
    
    cluster_rows = cluster_rows, 
    cluster_columns = cluster_columns, 
    
    row_dend_reorder = row_dend_reorder, 
    column_dend_reorder = column_dend_reorder,
    
    show_row_names = show_row_names, 
    show_column_names = show_column_names, 
    
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    
    row_names_max_width = row_names_max_width,
    column_names_max_height = column_names_max_height,
    
    ...)
  
  
  if(draw){
    ComplexHeatmap::draw(ht, column_title = title)
  }else{
    
    if(return == "ht"){
      return(ht)  
    }else{
      return(matrix)  
    }
    
  }
  
  
}


































