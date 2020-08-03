




# Cluster rows using hierarchical clustering

# d <- dist(matrix, method = "euclidean")
# cluster_rows <- hclust(d, method = "ward.D2")



# cols_lfc <- grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE)
# 
# mydata <- t(x[, cols_lfc, drop = FALSE])
# colnames(mydata) <- x[, gene_var]
# 
# consensus_clustering <- M3C::M3C(mydata, method = 2, silent = TRUE, seed = 123)
# 
# nr_clusters <- max(consensus_clustering$assignments)
# 
# gene_order <- rownames(consensus_clustering$realdataresults[[nr_clusters]]$ordered_annotation)
# 
# ### There is a bug in M3C and it removes 'X' from the names of genes. The code below is a way around.
# rownames(x) <- names(consensus_clustering$realdataresults[[nr_clusters]]$assignments)
# x <- x[gene_order, , drop = FALSE]
# rownames(x) <- NULL










# x = topTable_significant_genes;
# gene_var = "Endpoint";
# 
# lfc_prefix = "logFC"; pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val"; sep = "_";
# pval = 0.05;
# title = "";
# axis_text_x_angle = 30; axis_text_x_vjust = 1; axis_text_x_hjust = 1;
# axis_text_y_size = 10;
# color_low = '#42399B'; color_mid = "white"; color_high = '#D70131'; radius_range = c(10, 3)





#' Dot plot with logFC and p-values for mutiple contrasts
#' 
#' @param x TopTable
#' @export
wrapper_logFC_dotplot <- function(x, gene_var = "Hgnc_Symbol", lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val",  
  sep = "_", pval = 0.05, title = "", 
  axis_text_x_angle = 30, axis_text_x_vjust = 1, axis_text_x_hjust = 1, 
  axis_text_y_size = 10, title_size = 12,
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', 
  trim_limits = 2, radius_range = c(10, 3),
  legend_position = "right"){
  
  
  
  stopifnot(length(gene_var) == 1)
  
  
  data_lfc <- wrapper_extract_from_topTable(x, extract_prefix = lfc_prefix, sep = sep)
  
  data_pval <- wrapper_extract_from_topTable(x, extract_prefix = pval_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  contrasts <- colnames(data_lfc)
  
  
  data_lfc <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_lfc, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = contrasts, names_to = "contrast", values_to = lfc_prefix)
  
  data_pval <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_pval, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = contrasts, names_to = "contrast", values_to = pval_prefix)
  
  data_adjp <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE, check.names = FALSE), 
    cols = contrasts, names_to = "contrast", values_to = adjp_prefix)
  
  
  if(pval_prefix == adjp_prefix){
    data <- data_lfc %>% 
      dplyr::left_join(data_pval, by = c(gene_var, "contrast")) %>% 
      as.data.frame()
  }else{
    data <- data_lfc %>% 
      dplyr::left_join(data_pval, by = c(gene_var, "contrast")) %>% 
      dplyr::left_join(data_adjp, by = c(gene_var, "contrast")) %>% 
      as.data.frame()
  }
  
  
  data$significance <- factor(ifelse(data[, adjp_prefix] < pval, paste0("<", pval), paste0(">=", pval)), levels = paste0(c("<", ">="), pval))
  values_shape <- c(4, 32)
  names(values_shape) <- levels(data$significance)
  
  
  data$contrast <- factor(data$contrast, levels = contrasts)
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  data[, gene_var] <- factor(stringr::str_wrap(data[, gene_var], width = 70), levels = stringr::str_wrap(rev(x[, gene_var]), width = 70))
  
  
  pval_cut <- c(-1, 0.001, 0.01, 0.05, 0.1, 0.2, 1)
  data$pval_cut <- as.numeric(cut(data[, pval_prefix], breaks = pval_cut, labels = pval_cut[-1]))
  
  
  
  if(is.null(trim_limits)){
    max_abs_value <- ceiling(max(abs(range(data[, lfc_prefix], na.rm = TRUE))))
  }else if(trim_limits >= 1){
    max_abs_value <- trim_limits
  }else{
    ### Use quantiles 
    max_abs_value <- ceiling(max(abs(quantile(data[, lfc_prefix], probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
  }
  
  limits <- c(-max_abs_value, max_abs_value)
  
  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  
  ggp <- ggplot(data, aes(x = .data[["contrast"]], y = .data[[gene_var]], size = .data[["pval_cut"]], color = .data[[lfc_prefix]])) +
    geom_point() +
    geom_point(aes(size = .data[["pval_cut"]] + 1, shape = .data[["significance"]]), color = "black", show.legend = TRUE) +
    ggtitle(title) +
    theme_cowplot(12) +
    theme(plot.title = element_text(size = title_size),
      axis.line = element_blank(), 
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      axis.text.y = element_text(size = axis_text_y_size),
      legend.position = legend_position) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_radius(name = pval_prefix, range = radius_range, breaks = 1:(length(pval_cut) - 1), labels = pval_cut[-1], limits = c(1, length(pval_cut) - 1)) +
    scale_colour_gradient2(name = lfc_prefix, low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish) + 
    scale_shape_manual(name = adjp_prefix, values = values_shape, drop = FALSE)
  
  
  # This is great. squish in this context converts clamps all values to be within the min and max of the limits argument. i.e., if value < min(limits) then value = min(limits) else if value > max(limits) then value = max(limits).
  # scale_colour_gradient2(limits = c(-1.5, 1.5), oob = scales::squish)
  
  
  return(ggp)
  
  
  
}









# x = topTable_significant_genes; 
# gene_var = "Endpoint"; 
# 
# lfc_prefix = "logFC"; pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val";  
# sep = "_"; pval = 0.05; title = "";
# color_low = '#D70131'; color_mid = "white"; color_high = '#42399B';
# point_size = 3;
# cluster_rows = FALSE; 
# show_row_names = TRUE; show_column_names = TRUE; 
# row_names_gp = grid::gpar(fontsize = 8); column_names_gp = grid::gpar(fontsize = 8);
# row_split = NULL; column_split = NULL




#' Heatmap with logFC for mutiple contrasts
#' 
#' @param x TopTable
#' @export
wrapper_logFC_heatmap <- function(x, gene_var = "Hgnc_Symbol", 
  lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val",  
  sep = "_", draw = TRUE, title = "",
  color_low = '#42399B', color_mid = "white", color_high = '#D70131', trim_limits = NULL,
  rect_gp = grid::gpar(col = "grey"), point_size = 3,
  cluster_rows = FALSE, 
  row_split = NULL, column_split = NULL,
  show_row_names = TRUE, show_column_names = TRUE, 
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9), ...){
  
  
  stopifnot(length(gene_var) == 1)
  
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  rownames(x) <- stringr::str_wrap(x[, gene_var], width = 70)
  
  
  data_lfc <- wrapper_extract_from_topTable(x, extract_prefix = lfc_prefix, sep = sep)
  
  data_pval <- wrapper_extract_from_topTable(x, extract_prefix = pval_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  data_significance <- data.frame(apply(data_adjp, 2, function(x){
    
    pval_asterisk <- ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", ifelse(x < 0.1, ".", ""))))
    
    pval_asterisk
    
  }), stringsAsFactors = FALSE)
  
  
  
  matrix <- as.matrix(data_lfc)
  
  
  # ---------------------------------------------------------------------------
  # Colors for the heatmap expression
  # ---------------------------------------------------------------------------
  
  
  if(is.null(trim_limits)){
    max_abs_value <- ceiling(max(abs(range(matrix, na.rm = TRUE))))
  }else if(trim_limits >= 1){
    max_abs_value <- trim_limits
  }else{
    ### Use quantiles 
    max_abs_value <- ceiling(max(abs(quantile(matrix, probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
  }
  
  colors_matrix <- circlize::colorRamp2(c(-max_abs_value, 0, max_abs_value), c(color_low, color_mid, color_high))
  
  breaks <- seq(-max_abs_value, max_abs_value, 1)
  
  # ---------------------------------------------------------------------------
  # Display the significance levels
  # ---------------------------------------------------------------------------
  
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#customize-the-heatmap-body
  
  
  # cell_fun <- function(j, i, x, y, width, height, fill){ 
  #   grid::grid.text(data_significance[i, j], x, y, gp = grid::gpar(col = label_color, fontsize = label_size))
  # }
  
  
  
  cell_fun <- function(j, i, x, y, width, height, fill){ 
    
    if(data_significance[i, j] == "***"){
      grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "black"))
    }
    
    if(data_significance[i, j] == "**"){
      grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "dimgrey"))
    }
    
    if(data_significance[i, j] == "*"){
      grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "darkgrey"))
    }
    
    if(data_significance[i, j] == "."){
      grid::grid.points(x, y, pch = 21, size = grid::unit(point_size, "mm"), gp = grid::gpar(fill = "transparent"))
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Heatmap with the logFC
  # ---------------------------------------------------------------------------
  
  
  ht <- ComplexHeatmap::Heatmap(matrix, 
    name = lfc_prefix,
    col = colors_matrix, 
    na_col = "grey90",
    rect_gp = rect_gp,
    heatmap_legend_param = list(color_bar = "continuous", at = breaks), 
    cell_fun = cell_fun,
    
    cluster_columns = FALSE, 
    cluster_rows = cluster_rows, 
    
    row_dend_reorder = FALSE, 
    column_dend_reorder = FALSE,
    
    show_row_names = show_row_names, 
    show_column_names = show_column_names, 
    
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    
    row_split = row_split,
    column_split = column_split,
    
    ...)
  
  
  if(draw){
    ComplexHeatmap::draw(ht, column_title = title)
  }else{
    return(ht)  
  }
  
  
}





#' Heatmap with logFC or z-score normalized gene expression
#' 
#' 
#' @param x Matrix with expression to plot
#' @export
wrapper_gene_expression_heatmap <- function(x,  
  name = "sample\nlogFC", draw = TRUE, title = "", 
  color_low = 'cornflowerblue', color_mid = "royalblue4", color_high = '#ff717e', trim_limits = NULL,
  cluster_rows = FALSE, cluster_columns = FALSE,
  left_annotation = NULL, top_annotation = NULL, 
  row_split = NULL, column_split = NULL,
  row_title = NULL, column_title = NULL, 
  show_row_names = TRUE, show_column_names = FALSE, 
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9), ...){
  
  ### Colors for z-score
  # color_low = 'cornflowerblue', color_mid = "black", color_high = 'orangered'
  
  matrix <- as.matrix(x)
  
  
  # ---------------------------------------------------------------------------
  # Colors for the heatmap expression
  # ---------------------------------------------------------------------------
  
  
  if(is.null(trim_limits)){
    max_abs_value <- ceiling(max(abs(range(matrix, na.rm = TRUE))))
  }else if(trim_limits >= 1){
    max_abs_value <- trim_limits
  }else{
    ### Use quantiles 
    max_abs_value <- ceiling(max(abs(quantile(matrix, probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
  }
  
  
  colors_matrix <- circlize::colorRamp2(c(-max_abs_value, 0, max_abs_value), c(color_low, color_mid, color_high))
  
  breaks <- seq(-max_abs_value, max_abs_value, 1)
  
  
  # ---------------------------------------------------------------------------
  # Heatmap with the logFC
  # ---------------------------------------------------------------------------
  
  
  ht <- ComplexHeatmap::Heatmap(matrix, 
    name = name,
    col = colors_matrix, 
    na_col = "grey90",
    heatmap_legend_param = list(color_bar = "continuous", at = breaks), 
    
    cluster_columns = cluster_columns, 
    cluster_rows = cluster_rows, 
    
    row_dend_reorder = FALSE, 
    column_dend_reorder = FALSE,
    
    left_annotation = left_annotation,
    top_annotation = top_annotation,
    
    row_split = row_split,
    column_split = column_split,
    
    row_title = row_title,
    column_title = column_title,
    
    show_row_names = show_row_names, 
    show_column_names = show_column_names, 
    
    row_names_gp = row_names_gp,
    column_names_gp = column_names_gp,
    
    ...)
  
  
  if(draw){
    draw(ht, column_title = title)
  }else{
    return(ht)  
  }
  
  
}




































