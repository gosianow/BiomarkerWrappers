#' @include wrappers_logFC_plot.R 
NULL





###############################################################################
# Co-occurrence plots
###############################################################################



#' Compute odds ratios for pairs of covariates 
#' 
#' @param x Data frame with covariates where absence is encoded as 0 and presence is encoded as 1.
#' @export
or_matrix <- function(x){
  
  out_or <- matrix(Inf, nrow = ncol(x), ncol = ncol(x))
  rownames(out_or) <- colnames(x)
  colnames(out_or) <- colnames(x)
  
  out_pval <- matrix(Inf, nrow = ncol(x), ncol = ncol(x))
  rownames(out_pval) <- colnames(x)
  colnames(out_pval) <- colnames(x)
  
  for(i in 1:(ncol(x) - 1)){
    
    for(j in (i+1):ncol(x)){
      # i = 2; j = 3
      
      tbl <- table(x[, i], x[, j])
      
      ## Fisher's exact test
      prop_test_res <- NULL
      try(prop_test_res <- fisher.test(tbl), silent = TRUE)
      if(is.null(prop_test_res)){
        try(prop_test_res <- fisher.test(tbl, simulate.p.value = TRUE))
      }
      
      if(is.null(prop_test_res)){
        out_or[i, j] <- NA
        out_or[j, i] <- NA
        out_pval[i, j] <- NA
        out_pval[j, i] <- NA
      }else{
        out_or[i, j] <- prop_test_res$estimate
        out_or[j, i] <- prop_test_res$estimate
        out_pval[i, j] <- prop_test_res$p.value
        out_pval[j, i] <- prop_test_res$p.value
      }
      
      
    }# j
    
  }# i
  
  
  return(list(or = out_or, pval = out_pval))
  
  
}


#' Compute Jaccard similarity for pairs of covariates 
#' 
#' @param x Data frame with covariates where absence is encoded as 0 and presence is encoded as 1.
#' @export
jaccard_matrix <- function(x){
  
  out <- matrix(1, nrow = ncol(x), ncol = ncol(x))
  rownames(out) <- colnames(x)
  colnames(out) <- colnames(x)
  
  for(i in 1:(ncol(x) - 1)){
    
    for(j in (i+1):ncol(x)){
      # i = 2; j = 3
      
      tbl <- table(x[, i], x[, j])
      
      jaccard <- tbl[2, 2] / (tbl[1, 2] + tbl[2, 1] + tbl[2, 2])
      
      out[i, j] <- jaccard
      out[j, i] <- jaccard
      
      
    }# j
    
  }# i
  
  return(out)
  
  
}



#' Heatmap of Jaccard similarity between pairs of covariates
#' 
#' @param x Data frame with covariates where absence is encoded as 0 and presence is encoded as 1.
#' @export
wrapper_cooccurence_heatmap <- function(x, method = "Jaccard", column_title = "", colors = NULL, print_text = FALSE, text_fontsize = 5, text_col = "black", order = TRUE, draw = TRUE, title = ""){
  
  stopifnot(method %in% c("Jaccard", "LogOR"))
  
  
  ## -----------------------------
  ## Setup for LogOR
  ## -----------------------------
  
  
  if(method == "LogOR"){
    
    out <- or_matrix(x)
    
    ## OR
    xor <- out$or
    ## LogOR
    xlor <- log10(xor)
    
    xx <- xlor
    
    # xx[xx %in% c(-Inf, Inf)] <- NA
    
    name <- "LogOR"
    
    max_or <- round(max(abs(xx[! xx %in% c(-Inf, Inf)]), na.rm = TRUE) + 0.2, 1)
    half_max_or <- round(max_or/2, 1)
    
    legend_breaks = c(-max_or, -half_max_or, 0, half_max_or, max_or)
    
    if(is.null(colors)){
      colors <- c("steelblue", "ghostwhite", "violetred2")
    }
    
    color_heat <- circlize::colorRamp2(legend_breaks, colorRampPalette(colors)(length(legend_breaks)))
    
    
  }
  
  
  ## -----------------------------
  ## Setup for Jaccard similarity
  ## -----------------------------
  
  if(method == "Jaccard"){
    
    ## Jaccard similarity
    xj <- jaccard_matrix(x)
    
    xx <- xj
    
    name <- "Jaccard similarity"
    
    
    legend_breaks = seq(from = 0, to = 1, by = 0.2)
    
    if(is.null(colors)){
      # colors <- c("#4682b4", "white", "#dc143c")
      # colors <- c("royalblue4", "ghostwhite", "violetred2")
      colors <- c("royalblue", "darkgoldenrod1", "darkgoldenrod1", "darkorange2")
    }
    
    
    # color_heat <- colorRampPalette(colors)(length(legend_breaks))
    color_heat <- circlize::colorRamp2(legend_breaks, colorRampPalette(colors)(length(legend_breaks)))
    
  }
  
  
  ## -----------------------------
  
  
  if(order){
    d <- dist(t(x), method = "binary")
    cluster <- hclust(d, method = "ward.D2")
    genes <- colnames(x)
    gene_order <- genes[cluster$order]
    matrix <- xx[gene_order, gene_order]
  }else{
    matrix <- xx
  }
  
  
  cluster_rows <- FALSE
  cluster_columns <- FALSE
  
  
  ## -----------------------------
  
  if(print_text){
    cell_fun <- function(j, i, x, y, width, height, fill){ 
      # if(j < i)
      grid.text(round(matrix[i, j], 2), x, y, gp = gpar(col = text_col, fontsize = text_fontsize))
    }
  }else{
    cell_fun <- NULL
  }
  
  
  ## -----------------------------
  
  
  ht <- ComplexHeatmap::Heatmap(matrix, name = name, column_title = column_title, 
    col = color_heat,
    na_col = "grey90",
    cell_fun = cell_fun,
    show_row_names = TRUE, 
    show_column_names = TRUE,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_dend = TRUE, 
    show_column_dend = TRUE, 
    row_dend_reorder = FALSE, 
    column_dend_reorder = FALSE,
    row_names_side = "left", 
    column_names_side = "top",
    heatmap_legend_param = list(at = legend_breaks, labels = legend_breaks, color_bar = "continuous"))
  
  if(draw){
    draw(ht, column_title = title)
  }else{
    return(ht)  
  }
  
  
}









#' Dotplot with cooccurence odds ratios
#' 
#' @param x Data frame with covariates where absence is encoded as 0 and presence is encoded as 1.
#' @export
wrapper_cooccurence_dotplot <- function(x, order = TRUE, title = ""){
  
  
  ### Order using Jaccard similarity
  if(order){
    d <- dist(t(x), method = "binary")
    cluster <- hclust(d, method = "ward.D2")
    x <- x[cluster$order]
  }
  
  
  out <- or_matrix(x)
  
  ## LogOR
  xlor <- log10(out$or)
  
  xlor[xlor == Inf] <- 100
  xlor[xlor == -Inf] <- -100
  
  
  colnames(xlor) <- paste0("LogOR_", gsub("_", " ", colnames(xlor)))
  rownames(xlor) <- gsub("_", " ", rownames(xlor))
  
  xpval <- out$pval
  colnames(xpval) <- paste0("P.Value_", gsub("_", " ", colnames(xpval)))
  rownames(xpval) <- gsub("_", " ", rownames(xpval))
  
  
  xx <- data.frame(Gene = rownames(xlor), as.data.frame.matrix(xlor), as.data.frame.matrix(xpval), stringsAsFactors = FALSE, check.names = FALSE, row.names = NULL)
  
  
  ggp <- wrapper_logFC_dotplot(x = xx, gene_var = "Gene", lfc_prefix = "LogOR", pval_prefix = "P.Value", adjp_prefix = "P.Value", title = title, axis_text_x_angle = 90, axis_text_x_vjust = 0.5, axis_text_x_hjust = 0, axis_text_y_size = 11, radius_range = c(12, 5))
  
  
  ggp <- ggp +
    scale_x_discrete(position = "top") 
  
  ggp 
  
  
  return(ggp)
  
}

















