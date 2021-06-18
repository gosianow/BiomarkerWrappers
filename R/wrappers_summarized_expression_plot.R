



#' Dotplot of gene expression summarized per group
#' 
#' @param x Log transformed and normalized expression
#' @param group Vector with grouping
#' @export
wrapper_summarized_expression_dotplot <- function(x, group){
  
  
  x <- data.frame(group = group, t(x))
  group_var <- "group"
  
  
  ggdata_expr <- x %>% 
    dplyr::group_by_at(group_var) %>% 
    dplyr::summarize_all(list(~ mean(.))) %>% 
    tidyr::pivot_longer(!all_of(group_var), names_to = "Gene", values_to = "Mean_expr") %>% 
    data.frame
  
  
  ggdata_expr_no_zeros <- x %>% 
    dplyr::group_by_at(group_var) %>% 
    dplyr::summarize_all(list(~ mean(.[. != 0]))) %>% 
    tidyr::pivot_longer(!all_of(group_var), names_to = "Gene", values_to = "Mean_expr_no_zeros") %>% 
    data.frame
  
  
  ggdata_zeros <- x %>% 
    dplyr::group_by_at(group_var) %>% 
    dplyr::summarize_all(list(~ mean(. != 0))) %>% 
    tidyr::pivot_longer(!all_of(group_var), names_to = "Gene", values_to = "Detected") %>% 
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





#' Heatmap of gene expression summarized per group
#' 
#' @param x Matrix with gene expression.
#' @param group Vector of factors with grouping.
#' @param adjp Matrix with adjusted p-values for selected groups. Column names must correspond to the groups.
#' @param method Summarize and plot data 'asis' or transform it with 'z-score'.
#' @export
wrapper_summarized_expression_heatmap <- function(x, group, adjp = NULL,
  method = "z-score", scale = TRUE,
  title = "", name = NULL,
  centered = TRUE, palette = NULL, rev = FALSE,
  trim_values = NULL, trim_prop = NULL, trim_range = NULL, ceiling = FALSE,
  cell_fun_method = "none", font_size = 10,
  rect_gp = grid::gpar(col = "white"), 
  cluster_rows = FALSE, cluster_columns = FALSE, 
  row_dend_reorder = FALSE, column_dend_reorder = FALSE,
  show_row_names = TRUE, show_column_names = TRUE, 
  top_annotation = NULL,
  row_names_gp = grid::gpar(fontsize = 9), column_names_gp = grid::gpar(fontsize = 9), 
  row_names_max_width = unit(20, "cm"), column_names_max_height = unit(20, "cm"),
  row_names_width = 80, draw = TRUE, return = "ht", ...){
  
  
  stopifnot(is.factor(group))
  
  ### Remove unused levels
  group <- factor(group)
  
  
  stopifnot(method %in% c("asis", "z-score"))
  stopifnot(cell_fun_method %in% c("none", "circle", "rect"))
  
  
  if(is.null(name)){
    name <- ifelse(method == "z-score" || centered, "Mean\nz-score", "Mean\nexpr.")
  }
  
  if(is.null(trim_values)){
    if(method == "z-score" || centered){
      trim_values <- 1.5
    }else{
      trim_values <- NULL
    }
  }
  
  
  if(!is.null(adjp)){
    
    levels_group <- levels(group)
    
    stopifnot(length(intersect(levels_group, colnames(adjp))) > 0)
    
    stopifnot(all(rownames(adjp) == rownames(x)))
    
    adjp <- lapply(seq_along(levels_group), function(i){
      # i = 1
      
      if(levels_group[i] %in% colnames(adjp)){
        out <- adjp[, levels_group[i]]
      }else{
        out <- rep(NA, nrow(x))  
      }
      return(out)
      
    })
    adjp <- as.matrix(data.frame(adjp))
    rownames(adjp) <- rownames(x)
    colnames(adjp) <- levels_group
    
    ## Shorten the row names so they can be nicely displayed in the plots
    rownames(adjp) <- stringr::str_wrap(rownames(adjp), width = row_names_width)
    
  }
  
  ## Shorten the row names so they can be nicely displayed in the plots
  rownames(x) <- stringr::str_wrap(rownames(x), width = row_names_width)
  
  
  # ---------------------------------------------------------------------------
  # Prepare matrix
  # ---------------------------------------------------------------------------
  
  
  x <- as.matrix(x)
  
  if(method == "z-score"){
    
    zscore <- t(apply(x, 1, scale, center = TRUE, scale = scale))
    colnames(zscore) <- colnames(x)
    
    x <- zscore
    
    centered <- TRUE
    
  }
  
  
  x <- data.frame(group = group, t(x), check.names = FALSE)
  group_var <- "group"
  
  
  ggdata <- x %>% 
    dplyr::group_by_at(group_var) %>% 
    dplyr::summarize_all(list(~ mean(.))) %>% 
    data.frame(check.names = FALSE)
  
  
  rownames(ggdata) <- ggdata[, group_var]
  
  ggdata[, group_var] <- NULL
  
  matrix <- as.matrix(t(ggdata))
  
  
  # ---------------------------------------------------------------------------
  # Prepare top annotation that shows number of samples in groups
  # ---------------------------------------------------------------------------
  
  
  if(is.logical(top_annotation)){
    
    if(top_annotation){
      
      tbl <- as.numeric(table(group))
      
      top_annotation <- HeatmapAnnotation(
        "Sample number" = row_anno_barplot(x = tbl, border = TRUE, axis = TRUE, axis_param = list(gp = gpar(fontsize = 5)), gp = gpar(fill = "#696969", col = "#696969")), 
        which = "column",
        show_annotation_name = FALSE)
      
      
    }else{
      top_annotation <- NULL
    }
    
  }
  
  
  
  # ---------------------------------------------------------------------------
  # Prepare data_zeros
  # ---------------------------------------------------------------------------
  
  
  if(cell_fun_method %in% c("circle", "rect")){
    
    ggdata <- x %>% 
      dplyr::group_by_at(group_var) %>% 
      dplyr::summarize_all(list(~ mean(. != 0))) %>% 
      data.frame(check.names = FALSE)
    
    rownames(ggdata) <- ggdata[, group_var]
    
    ggdata[, group_var] <- NULL
    
    data_zeros <- as.matrix(t(ggdata))
    
  }
  
  
  
  # ---------------------------------------------------------------------------
  # Prepare data_significance
  # ---------------------------------------------------------------------------
  
  
  if(!is.null(adjp)){
    
    
    data_significance <- data.frame(apply(adjp, 2, function(x){
      
      x[is.na(x)] <- 1
      
      pval_asterisk <- ifelse(x < 0.001, "***", ifelse(x < 0.01, "**", ifelse(x < 0.05, "*", ifelse(x < 0.1, ".", ""))))
      
      pval_asterisk
      
    }), stringsAsFactors = FALSE)
    
    
  }
  
  
  # ---------------------------------------------------------------------------
  # Colors for the heatmap expression
  # ---------------------------------------------------------------------------
  
  
  colors_matrix <- format_colors_num(x = matrix, centered = centered, palette = palette, rev = rev, trim_values = trim_values, trim_prop = trim_prop, trim_range = trim_range, ceiling = ceiling)
  
  
  
  # ---------------------------------------------------------------------------
  # cell_fun
  # ---------------------------------------------------------------------------
  
  
  if(cell_fun_method == "none"){
    
    if(is.null(adjp)){
      
      cell_fun <- NULL
      
    }else{
      
      
      cell_fun <- function(j, i, x, y, width, height, fill){
        
        
        # grid::grid.text(data_significance[i, j], x, y = y, gp = grid::gpar(col = "black", fontsize = font_size))
        
        ### Center the displayed items
        
        if(data_significance[i, j] == "."){
          grid::grid.text(data_significance[i, j], x, y = y + (0.3 * unit(font_size, units = "points")), gp = grid::gpar(col = "black", fontsize = font_size))
          
        }else{
          grid::grid.text(data_significance[i, j], x, y = y - (0.2 * unit(font_size, units = "points")), gp = grid::gpar(col = "black", fontsize = font_size))
          
        }
        
        
      }
      
      
      
    }
    
    
  }else if(cell_fun_method == "circle"){
    
    rect_gp <- gpar(type = "none")
    
    cell_fun <- function(j, i, x, y, width, height, fill){
      
      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
      
      ### Scale the radius
      
      # grid.circle(x = x, y = y, r = min(unit.c(width, height))/2 * (data_zeros[i, j] + 0.01), 
      #   gp = gpar(fill = colors_matrix(matrix[i, j]), col = NA))
      
      ### Scale the area
      
      grid.circle(x = x, y = y, r = min(unit.c(width, height))/2 * sqrt(data_zeros[i, j]), 
        gp = gpar(fill = colors_matrix(matrix[i, j]), col = NA))
      
    }
    
  }else if(cell_fun_method == "rect"){
    
    rect_gp <- gpar(type = "none")
    
    cell_fun <- function(j, i, x, y, width, height, fill){
      
      grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
      
      grid.rect(x = x - (width * (1 - data_zeros[i, j]) / 2), y = y, width = width * data_zeros[i, j], height = height, 
        gp = gpar(col = NA, fill = colors_matrix(matrix[i, j])))
      
    }
    
  }
  
  
  
  
  ### ------------------------------------------------------------------------
  ### Heatmap
  ### ------------------------------------------------------------------------
  
  
  
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
    
    top_annotation = top_annotation,
    
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
















