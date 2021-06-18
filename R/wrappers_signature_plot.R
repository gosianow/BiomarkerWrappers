


#' Heatmap with signature and expression of genes defining that signature
#' 
#' @param x Matrix with log transformed and normalized gene expression. Genes in rows and samples in columns.
#' @param signature Vector with precalculated signature.
#' @param title Title.
#' @export
wrapper_signature_heatmap <- function(x, signature, title, plot_expr = FALSE, show_row_names = FALSE){
  
  
  zscore <- apply(x, 1, scale, center = TRUE, scale = TRUE)
  colnames(zscore) <- rownames(x)
  rownames(zscore) <- colnames(x)
  
  x <- t(x)
  
  ### Order samples by signature 
  
  indx <- order(signature, decreasing = TRUE)
  
  x <- x[indx, , drop = FALSE]
  zscore <- zscore[indx, , drop = FALSE]
  signature <- signature[indx]
  
  
  cluster_rows = FALSE
  
  # Cluster genes using hierarchical clustering
  if(ncol(zscore) > 2){
    d <- dist(t(zscore), method = "euclidean")
    cluster_columns <- hclust(d, method = "ward.D2")
  }else{
    cluster_columns <- FALSE
  }
  
  
  # Colors for the heatmap
  
  color_zscore <- format_colors_num(zscore, trim_values = 2.5)
  
  
  ### Signature - barplot
  ha_bar <- rowAnnotation("Signature" = row_anno_barplot(
    x = signature, border = FALSE, axis = TRUE, 
    axis_param = list(gp = gpar(fontsize = 5)), gp = gpar(fill = "#696969", col = "#696969"), 
    bar_width = 0.9), width = unit(2, "cm"), show_annotation_name = TRUE, 
    annotation_name_rot = 0, annotation_name_offset = unit(5, "mm"), 
    annotation_name_gp = gpar(fontsize = 7))
  
  
  ### Heatmap with the z-score of signature genes
  ht2 <- Heatmap(zscore, name = "Z-score", column_title = title, 
    col = color_zscore, cluster_columns = cluster_columns, cluster_rows = cluster_rows, 
    row_dend_reorder = FALSE, heatmap_legend_param = list(color_bar = "continuous"), 
    show_row_names = show_row_names,
    row_dend_width = unit(2, "cm"), 
    rect_gp = gpar(col = NA), column_names_gp = gpar(fontsize = 8))
  
  
  if(!plot_expr){
    
    draw_out <- ht2 + ha_bar
    
  }else{
    
    # Colors for the heatmap
    
    color_expr <- format_colors_num(x, centered = FALSE)
    
    ### Heatmap with the expression of signature genes
    ht1 <- Heatmap(x, name = "Expr.",
      col = color_expr, cluster_columns = cluster_columns, cluster_rows = cluster_rows,
      row_dend_reorder = FALSE, heatmap_legend_param = list(color_bar = "continuous"),
      show_row_names = FALSE, 
      row_dend_width = unit(2, "cm"),
      rect_gp = gpar(col = NA), column_names_gp = gpar(fontsize = 8))
    
    draw_out <- ht1 + ht2 + ha_bar
    
    
  }
  
  
  ComplexHeatmap::draw(draw_out, auto_adjust = FALSE)
  
  
}














