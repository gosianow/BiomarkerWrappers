







#' Heatmap with logFC or z-score normalized gene expression
#' 
#' 
#' @param x Matrix with expression to plot
#' @param method Plot data 'asis' or transform it with 'z-score'.
#' @export
wrapper_gene_expression_heatmap <- function(x,  
  method = "z-score", scale = TRUE,
  title = "", name = "z-score", 
  centered = TRUE, palette = NULL, rev = FALSE,
  trim_values = 2.5, trim_prop = NULL, trim_range = NULL, ceiling = FALSE,
  rect_gp = grid::gpar(col = NA), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_dend_reorder = FALSE, column_dend_reorder = FALSE,
  row_split = NULL, column_split = NULL,
  show_row_names = TRUE, show_column_names = FALSE, 
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9), 
  row_names_max_width = unit(20, "cm"), column_names_max_height = unit(20, "cm"),
  row_names_width = 80, draw = TRUE, return = "ht", ...){
  
  
  stopifnot(method %in% c("asis", "z-score"))
  
  
  ## Shorten the row names so they can be nicely displayed in the plots
  rownames(x) <- stringr::str_wrap(rownames(x), width = row_names_width)
  
  
  # ---------------------------------------------------------------------------
  # Prepare matrix
  # ---------------------------------------------------------------------------
  
  
  matrix <- as.matrix(x)
  
  if(method == "z-score"){
    
    zscore <- t(apply(matrix, 1, scale, center = TRUE, scale = scale))
    colnames(zscore) <- colnames(matrix)
    
    matrix <- zscore
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Colors for the heatmap expression
  # ---------------------------------------------------------------------------
  
  
  colors_matrix <- format_colors_num(x = matrix, centered = centered, palette = palette, rev = rev, trim_values = trim_values, trim_prop = trim_prop, trim_range = trim_range, ceiling = ceiling)
  
  
  # ---------------------------------------------------------------------------
  # Heatmap with the logFC
  # ---------------------------------------------------------------------------
  
  
  ht <- ComplexHeatmap::Heatmap(matrix, 
    name = name,
    col = colors_matrix, 
    na_col = "grey90",
    rect_gp = rect_gp, 
    heatmap_legend_param = list(color_bar = "continuous"), 
    
    cluster_rows = cluster_rows, 
    cluster_columns = cluster_columns, 
    
    row_dend_reorder = row_dend_reorder, 
    column_dend_reorder = column_dend_reorder,
    
    row_split = row_split,
    column_split = column_split,
    
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




































