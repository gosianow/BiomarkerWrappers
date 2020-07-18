




# x <- topTable
# genesets; contrast; gene_var = "EntrezIDs"; 
# statistic_prefix = "t"; sep = "_";
# trim_limits = 0.01;
# color_low = '#42399B'; color_mid = "darkgrey"; color_high = '#D70131';
# title = ""; title_size = 10; title_width = 100; axis_text_y_size = 8; axis_text_y_width = 70




#' Plot GSEA statistics ranks
#' 
#' @param x TopTable with DGE results, for example, from limma.
#' @export
wrapper_plot_GSEA <- function(x, contrast, genesets, gene_var = "EntrezIDs", statistic_prefix = "t", sep = "_", 
  gsea_results = NULL, geneset_var = "Geneset", adjp_var = "adj.P.Val", color_point_var = "NES",
  trim_limits = 0.01,
  color_low = '#42399B', color_mid = "darkgrey", color_high = '#D70131',
  title = "", title_size = 10, title_width = 100, axis_text_y_size = 8, axis_text_y_width = 70){
  
  
  ## Wrap the title so it can be nicely displayed in the plots
  title <- stringr::str_wrap(title, width = title_width)
  
  
  # -------------------------------------------------------------------------
  # Prepare data
  # -------------------------------------------------------------------------
  
  statistic <- x[, paste0(statistic_prefix, sep, contrast)]
  names(statistic) <- x[, gene_var]
  
  statistic <- sort(statistic, decreasing = TRUE)
  
  
  data <- data.frame(universe = names(statistic), statistic = statistic, rank = rank(-statistic), 
    stringsAsFactors = FALSE)
  
  data$direction <- factor(ifelse(data$statistic > 0, "up", "down"), levels = c("up", "down"))
  
  
  
  ### Trim limits
  
  if(is.null(trim_limits)){
    max_abs_value <- ceiling(max(abs(range(data[, "statistic"], na.rm = TRUE))))
  }else if(trim_limits >= 1){
    max_abs_value <- trim_limits
  }else{
    ### Use quantiles 
    max_abs_value <- ceiling(max(abs(quantile(data[, "statistic"], probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
  }
  
  
  limits <- c(-max_abs_value, max_abs_value)
  
  
  # data$statistic[data$statistic > max_abs_value] <- max_abs_value
  # data$statistic[data$statistic < -max_abs_value] <- -max_abs_value
  
  ### Normalize the hight of lines to 1
  data$statistic_adj <- data$statistic / (2 * max_abs_value)
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  genesets <- lapply(genesets, unique)
  
  # genesets[[1]]
  # statistic[genesets[[1]]]
  
  
  ### Keep genes in the gene sets that are in the universe
  
  ## We want to keep the order that is in the univarse so it has to be the first argument in intersect
  genesets <- lapply(genesets, function(x){intersect(data$universe, x)})
  
  # genesets[[1]]
  # statistic[genesets[[1]]]
  
  
  if(length(genesets) == 0){
    message("There are no common genes between the universe and the genesets.")
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # Data for ggplot
  # -------------------------------------------------------------------------
  
  
  ggdata <- lapply(seq_along(genesets), function(i){
    # i = 1
    
    mm <- match(genesets[[i]], data$universe)
    
    ggdata <- data[mm, , drop = FALSE]
    
    ggdata$Geneset <- names(genesets[i])
    
    return(ggdata)
    
  })
  
  
  ggdata <- plyr::rbind.fill(ggdata)
  
  
  ## Wrap the gene set names so they can be nicely displayed in the plots
  ggdata$Geneset <- factor(ggdata$Geneset, levels = rev(names(genesets)), labels = stringr::str_wrap(rev(names(genesets)), width = axis_text_y_width))
  
  ggdata$Geneset_num <- as.numeric(ggdata$Geneset)
  
  ggdata$statistic_adj2 <- ggdata$statistic_adj + ggdata$Geneset_num + ifelse(ggdata$direction == "up", 0.1, -0.1)
  
  ggdata$statistic_fixed <- ggdata$Geneset_num + ifelse(ggdata$direction == "up", 0.5, -0.5)
  
  
  # -------------------------------------------------------------------------
  # ggplot
  # -------------------------------------------------------------------------
  
  
  xlim <- c(0, max(data$rank))
  xlab <- "Rank"
  
  
  # ggp <-  ggplot(ggdata, aes(x = .data$rank, xend = .data$rank, y = .data$Geneset, yend = .data$statistic_fixed, color = .data$direction))
  
  ggp <-  ggplot(ggdata, aes(x = 0, y = .data$Geneset, group = .data$Geneset)) +
    geom_line() +
    geom_segment(aes(x = .data$rank, xend = .data$rank, y = .data$Geneset_num - 0.45, yend = .data$Geneset_num + 0.45, color = .data$statistic)) +
    ggtitle(title) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      # axis.title.x = element_blank(), 
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      plot.title = element_text(size = title_size, hjust = 1),
      legend.position = "none") +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "y", minor = "none", size.major = 0.25) +
    scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish) + 
    scale_x_continuous(expand = c(0, 0))
  
  
  # ggp
  
  
  
  
  ggp2 <- ggplot(data, aes(x = .data$rank, xend = .data$rank, y = 0, yend = .data$statistic, color = .data$statistic)) +
    geom_segment() +
    xlab(xlab) +
    ylab(statistic_prefix) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      # axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      plot.title = element_text(size = axis_text_y_size, hjust = 0.5, face = "plain"),
      legend.position = "none") +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 0.5, remove = FALSE) +
    background_grid(major = "y", minor = "none", size.major = 0.25) +
    scale_colour_gradient2(low = color_low, mid = color_mid, high = color_high, midpoint = 0, limits = limits, oob = scales::squish) + 
    scale_x_continuous(expand = c(0, 0))
  
  
  # ggp2
  
  out1 <- cowplot::plot_grid(ggp, ggp2, nrow = 2, ncol = 1, rel_heights = c(max(c(1, length(genesets) / 4)), 1), align = "v", axis = 'lr')
  
  
  
  if(!is.null(gsea_results)){
    
    if(all(names(genesets) == gsea_results[, geneset_var])){
      
      ### Use a trick with the title. Otherwise, the plots are not aligned :/
      
      ggp3 <- wrapper_plot_ORA_dotplot_single(gsea_results, geneset_var = geneset_var, observed_var = NULL, adjp_var = adjp_var, color_point_var = color_point_var, title = paste0(" ", paste0(rep("\n", stringr::str_count(title, "\n")), collapse = " ")," "), title_width = 0, title_size = title_size, size_range = c(2, 5)) +
        theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
      
      
      out2 <- cowplot::plot_grid(ggp3, nrow = 2, ncol = 1, rel_heights = c(max(c(1, length(genesets) / 4)), 1))
      
      
      out <- cowplot::plot_grid(out1, out2, nrow = 1, ncol = 2, rel_widths = c(4, 1), align = "h", axis = "tb")
      
      out
      
      
    }else{
      
      out1
      
    }
    
    
  }else{
    
    out1
    
  }
  
  
  
  
}





























