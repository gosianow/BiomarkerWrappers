




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
# color_low = '#D70131'; color_mid = "white"; color_high = '#42399B'; radius_range = c(10, 3)





#' Dot plot with logFC and p-values for mutiple contrasts
#' 
#' @param x TopTable
wrapper_plot_logFC_dotplot <- function(x, gene_var = "Hgnc_Symbol", lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val",  
  sep = "_", pval = 0.05, title = "", 
  axis_text_x_angle = 30, axis_text_x_vjust = 1, axis_text_x_hjust = 1, 
  axis_text_y_size = 10, 
  color_low = '#D70131', color_mid = "white", color_high = '#42399B', 
  radius_range = c(10, 3)){
  
  
  
  stopifnot(length(gene_var) == 1)
  
  
  
  
  data_lfc <- wrapper_extract_from_topTable(x, extract_prefix = lfc_prefix, sep = sep)
  
  data_pval <- wrapper_extract_from_topTable(x, extract_prefix = pval_prefix, sep = sep)
  
  data_adjp <- wrapper_extract_from_topTable(x, extract_prefix = adjp_prefix, sep = sep)
  
  
  contrasts <- colnames(data_lfc)
  
  
  data_lfc <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_lfc, stringsAsFactors = FALSE), 
    cols = contrasts, names_to = "contrast", values_to = lfc_prefix)
  
  data_pval <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_pval, stringsAsFactors = FALSE), 
    cols = contrasts, names_to = "contrast", values_to = pval_prefix)
  
  data_adjp <- pivot_longer(data.frame(x[, gene_var, drop = FALSE], data_adjp, stringsAsFactors = FALSE), 
    cols = contrasts, names_to = "contrast", values_to = adjp_prefix)
  
  
  data <- data_lfc %>% 
    left_join(data_pval, by = c(gene_var, "contrast")) %>% 
    left_join(data_adjp, by = c(gene_var, "contrast")) %>% 
    as.data.frame()
  
  
  data$significance <- ifelse(data[, adjp_prefix] < pval, paste0("<", pval), paste0(">=", pval))
  values_shape <- c(4, 32)
  
  data$contrast <- factor(data$contrast, levels = contrasts)
  
  ## Shorten the gene set names so they can be nicely displayed in the plots
  data[, gene_var] <- factor(stringr::str_wrap(data[, gene_var], width = 70), levels = stringr::str_wrap(rev(x[, gene_var]), width = 70))
  
  
  pval_cut <- c(-1, 0.001, 0.01, 0.05, 0.1, 0.2, 1)
  data$pval_cut <- as.numeric(cut(data[, pval_prefix], breaks = pval_cut, labels = pval_cut[-1]))
  
  
  # ---------------------------------------------------------------------------
  # ggplot
  # ---------------------------------------------------------------------------
  
  
  ggp <- ggplot(data, aes_string(x = "contrast", y = gene_var, size = "pval_cut", color = lfc_prefix)) +
    geom_point() +
    geom_point(aes_string(size = "pval_cut", shape = "significance"), color = "black", show.legend = TRUE) +
    ggtitle(title) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = axis_text_x_angle, vjust = axis_text_x_vjust, hjust = axis_text_x_hjust),
      axis.text.y = element_text(size = axis_text_y_size)) +
    panel_border(colour = "black", linetype = 1, size = 1, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_radius(name = pval_prefix, range = radius_range, breaks = 1:(length(pval_cut) - 1), labels = pval_cut[-1], limits = c(1, length(pval_cut) - 1)) +
    scale_colour_gradient2(name = lfc_prefix, low = color_low, mid = color_mid, high = color_high) + 
    scale_shape_manual(name = adjp_prefix, values = values_shape)
  
  
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
# row_names_gp = gpar(fontsize = 8); column_names_gp = gpar(fontsize = 8);
# row_split = NULL; column_split = NULL




#' Heatmap with logFC and p-values for mutiple contrasts
#' 
#' @param x TopTable
wrapper_plot_logFC_heatmap <- function(x, gene_var = "Hgnc_Symbol", 
  lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val",  
  sep = "_", pval = 0.05, title = "",
  color_low = '#D70131', color_mid = "white", color_high = '#42399B',
  point_size = 3,
  cluster_rows = FALSE, 
  show_row_names = TRUE, show_column_names = TRUE, 
  row_names_gp = gpar(fontsize = 9), column_names_gp = gpar(fontsize = 9),
  row_split = NULL, column_split = NULL, ...){
  
  
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
  
  max_abs_value <- ceiling(max(abs(range(matrix, na.rm = TRUE))))
  
  colors_matrix <- circlize::colorRamp2(c(-max_abs_value, 0, max_abs_value), c(color_low, color_mid, color_high))
  
  breaks <- seq(-max_abs_value, max_abs_value, 1)
  
  # ---------------------------------------------------------------------------
  # Display the significance levels
  # ---------------------------------------------------------------------------
  
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html#customize-the-heatmap-body
  
  
  # cell_fun <- function(j, i, x, y, width, height, fill){ 
  #   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
  #   grid.text(data_significance[i, j], x, y, gp = gpar(col = label_color, fontsize = label_size))
  # }
  
  
  
  cell_fun <- function(j, i, x, y, width, height, fill){ 
    
    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
    
    
    if(data_significance[i, j] == "***"){
      grid.points(x, y, pch = 21, size = unit(point_size, "mm"), gp = gpar(fill = "black"))
    }
    
    if(data_significance[i, j] == "**"){
      grid.points(x, y, pch = 21, size = unit(point_size, "mm"), gp = gpar(fill = "dimgrey"))
    }
    
    if(data_significance[i, j] == "*"){
      grid.points(x, y, pch = 21, size = unit(point_size, "mm"), gp = gpar(fill = "darkgrey"))
    }
    
    if(data_significance[i, j] == "."){
      grid.points(x, y, pch = 21, size = unit(point_size, "mm"), gp = gpar(fill = "transparent"))
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Heatmap with the logFC
  # ---------------------------------------------------------------------------
  
  
  ht <- ComplexHeatmap::Heatmap(matrix, 
    name = "logFC",
    col = colors_matrix, 
    na_col = "grey90",
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
  
  
  draw(ht)
  
  
  
}






































