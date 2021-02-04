



#' Dotplot of scRNA-seq expression
#' 
#' @param x Log transformed and normalized expression
#' @param group Vector with grouping
#' @export
wrapper_scDot_plot <- function(x, group){
  
  
  x <- data.frame(group = group, t(x))
  group_var <- "group"
  
  
  ggdata_expr <- x %>% 
    group_by_at(group_var) %>% 
    summarize_all(list(~ mean(.))) %>% 
    pivot_longer(!all_of(group_var), names_to = "Gene", values_to = "Mean_expr") %>% 
    data.frame
  
  
  ggdata_expr_no_zeros <- x %>% 
    group_by_at(group_var) %>% 
    summarize_all(list(~ mean(.[. != 0]))) %>% 
    pivot_longer(!all_of(group_var), names_to = "Gene", values_to = "Mean_expr_no_zeros") %>% 
    data.frame
  
  
  ggdata_zeros <- x %>% 
    group_by_at(group_var) %>% 
    summarize_all(list(~ mean(. != 0))) %>% 
    pivot_longer(!all_of(group_var), names_to = "Gene", values_to = "Detected") %>% 
    data.frame
  
  
  ggdata <- full_join(ggdata_expr, ggdata_zeros, by = c(group_var, "Gene")) %>% 
    full_join(ggdata_expr_no_zeros, by = c(group_var, "Gene"))
  
  
  ggdata$Gene <- factor(ggdata$Gene, levels = rev(levels(factor(ggdata$Gene))))
  
  
  
  ### ------------------------------------------------------------------------
  ### Dot plot
  ### ------------------------------------------------------------------------
  
  
  
  ggplot(ggdata, aes(x = .data[[group_var]], y = .data[["Gene"]], size = .data[["Detected"]], color = .data[["Mean_expr"]])) +
    geom_point() +
    scale_size_area(max_size = 15, breaks = seq(0, 1, by = 0.25), limits = c(0, 1)) +
    scale_colour_gradientn(colours = hcl.colors(20, palette = "viridis"), limits = c(0, max(ggdata[, "Mean_expr"])), oob = scales::squish) +
    theme(axis.line = element_blank()) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25)
  
  
  
  
}





#' Dotplot heatmap of scRNA-seq expression
#' 
#' @param x Log transformed and normalized expression
#' @param group Vector with grouping
#' @export
wrapper_scDot_heatmap <- function(x, group, 
  method = "asis", scale = FALSE,
  cell_fun_method = "circle",
  name = "Mean expr.", title = "", draw = TRUE,
  colors = NULL, trim_limits = NULL,
  rect_gp = grid::gpar(col = "black"),
  cluster_rows = TRUE, cluster_columns = FALSE,
  row_dend_reorder = FALSE, column_dend_reorder = FALSE, 
  clustering_distance_rows = "euclidean",
  clustering_method_rows = "complete",
  clustering_distance_columns = "euclidean",
  clustering_method_columns = "complete",
  show_row_names = TRUE, show_column_names = TRUE, 
  row_names_width = 80, ...){
  
  
  x <- data.frame(group = group, t(x))
  group_var <- "group"
  
  
  ggdata_expr <- x %>% 
    group_by_at(group_var) %>% 
    summarize_all(list(~ mean(.))) %>% 
    data.frame
  
  rownames(ggdata_expr) <- ggdata_expr[, group_var]
  
  ggdata_expr[, group_var] <- NULL
  
  mat_expr <- as.matrix(t(ggdata_expr))
  
  
  if(method == "z-score"){
    
    zscore <- t(apply(mat_expr, 1, scale, center = TRUE, scale = scale))
    colnames(zscore) <- colnames(mat_expr)
    
    mat_expr <- zscore
    
  }
  
  
  ggdata_zeros <- x %>% 
    group_by_at(group_var) %>% 
    summarize_all(list(~ mean(. != 0))) %>% 
    data.frame
  
  rownames(ggdata_zeros) <- ggdata_zeros[, group_var]
  
  ggdata_zeros[, group_var] <- NULL
  
  mat_zeros <- as.matrix(t(ggdata_zeros))
  
  
  
  ## Shorten the row names so they can be nicely displayed in the plots
  rownames(mat_expr) <- stringr::str_wrap(rownames(mat_expr), width = row_names_width)
  
  # ---------------------------------------------------------------------------
  # Colors for the heatmap expression
  # ---------------------------------------------------------------------------
  
  
  
  if(is.null(trim_limits)){
    max_abs_value <- ceiling(max(abs(range(mat_expr, na.rm = TRUE))))
  }else if(trim_limits >= 1){
    max_abs_value <- trim_limits
  }else{
    ### Use quantiles 
    max_abs_value <- ceiling(max(abs(quantile(mat_expr, probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
  }
  
  
  
  if(is.null(colors)){
    
    if(method == "z-score"){
      col_fun <- format_colors_num(mat_expr, trim_values = max_abs_value)
    }else{
      col_fun <- format_colors_num(mat_expr, centered = FALSE)
    }
    
  }else if(is.function(colors)){
    
    col_fun <- colors
    
  }else{
    
    if(method == "z-score"){
      col_fun <- format_colors_num(mat_expr, trim_values = max_abs_value, palette = colors)
    }else{
      col_fun <- format_colors_num(mat_expr, centered = FALSE, palette = colors)
    }
    
  }
  
  
  
  
  
  ### ------------------------------------------------------------------------
  ### Heatmap
  ### ------------------------------------------------------------------------
  
  
  
  
  
  
  if(cell_fun_method == "none"){
    
    cell_fun <- NULL
    
  }
  
  
  if(cell_fun_method == "circle"){
    
    rect_gp <- gpar(type = "none")
    
    cell_fun <- function(j, i, x, y, width, height, fill){
      
      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
      
      ### Scale the radius
      
      # grid.circle(x = x, y = y, r = min(unit.c(width, height))/2 * (mat_zeros[i, j] + 0.01), 
      #   gp = gpar(fill = col_fun(mat_expr[i, j]), col = NA))
      
      ### Scale the area
      
      grid.circle(x = x, y = y, r = min(unit.c(width, height))/2 * sqrt(mat_zeros[i, j]), 
        gp = gpar(fill = col_fun(mat_expr[i, j]), col = NA))
      
    }
    
  }
  
  
  if(cell_fun_method == "rect"){
    
    rect_gp <- gpar(type = "none")
    
    cell_fun <- function(j, i, x, y, width, height, fill){
      
      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
      
      grid.rect(x = x - (width * (1 - mat_zeros[i, j]) / 2), y = y, width = width * mat_zeros[i, j], height = height, 
        gp = gpar(col = NA, fill = col_fun(mat_expr[i, j])))
      
    }
    
  }
  
  
  
  
  
  ht <- Heatmap(mat_expr, name = name, col = col_fun, 
    rect_gp = rect_gp, cell_fun = cell_fun, 
    
    cluster_rows = cluster_rows, 
    cluster_columns = cluster_columns,
    
    row_dend_reorder = row_dend_reorder, 
    column_dend_reorder = column_dend_reorder,
    
    clustering_distance_rows = clustering_distance_rows,
    clustering_method_rows = clustering_method_rows,
    clustering_distance_columns = clustering_distance_columns,
    clustering_method_columns = clustering_method_columns,
    
    show_row_names = show_row_names, 
    show_column_names = show_column_names, 
    
    ...)
  
  
  
  if(draw){
    draw(ht, column_title = title)
  }else{
    return(ht)  
  }
  
  
}
















