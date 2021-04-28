




#' Plot GSEA statistics ranks
#' 
#' @param statistic Named vector with statistic. Names correspond to genes.
#' @param genesets Named list with gene sets to plot. Names correspond to gene sets.
#' @param adjp Vector with adjusted p-values for the gene sets.
#' @param enrichment_score Vector with enrichment scores for the gene sets.
#' @export
wrapper_gsea_plot_core <- function(statistic, genesets, adjp = NULL, enrichment_score = NULL,
  statistic_name = "t", geneset_name = "GenesetID", adjp_name = "adj.P.Val", enrichment_score_name = "NES",
  trim_limits = 0.01,
  color_low = '#42399B', color_mid = "darkgrey", color_high = '#D70131',
  title = "", title_size = 10, title_width = 100, axis_text_y_size = 8, axis_text_y_width = 80){
  
  
  stopifnot(!is.null(names(statistic)))
  
  ## Wrap the title so it can be nicely displayed in the plots
  title <- stringr::str_wrap(title, width = title_width)
  
  
  # -------------------------------------------------------------------------
  # Prepare data
  # -------------------------------------------------------------------------
  
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
  
  
  ### Normalize the hight of lines to 1 (-0.5;0.5)
  
  # data$statistic[data$statistic > max_abs_value] <- max_abs_value
  # data$statistic[data$statistic < -max_abs_value] <- -max_abs_value
  # data$statistic_adj <- data$statistic / (2 * max_abs_value)
  
  
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
    stop("There are no common genes between the universe and the genesets.")
  }
  
  
  # -------------------------------------------------------------------------
  # Data for ggplot
  # -------------------------------------------------------------------------
  
  
  ggdata <- lapply(seq_along(genesets), function(i){
    # i = 1
    
    mm <- match(genesets[[i]], data$universe)
    
    ggdata <- data[mm, , drop = FALSE]
    
    ggdata$GenesetID <- names(genesets[i])
    
    return(ggdata)
    
  })
  
  
  ggdata <- plyr::rbind.fill(ggdata)
  
  
  ## Wrap the gene set names so they can be nicely displayed in the plots
  ggdata$GenesetID <- factor(ggdata$GenesetID, levels = rev(names(genesets)), labels = stringr::str_wrap(rev(names(genesets)), width = axis_text_y_width))
  
  ggdata$GenesetID_num <- as.numeric(ggdata$GenesetID)
  
  # ggdata$statistic_adj2 <- ggdata$statistic_adj + ggdata$GenesetID_num + ifelse(ggdata$direction == "up", 0.1, -0.1)
  
  ggdata$statistic_fixed <- ggdata$GenesetID_num + ifelse(ggdata$direction == "up", 0.5, -0.5)
  
  
  # -------------------------------------------------------------------------
  # ggplot
  # -------------------------------------------------------------------------
  
  
  xlim <- c(0, max(data$rank))
  xlab <- "Rank"
  
  
  # ggp <-  ggplot(ggdata, aes(x = .data$rank, xend = .data$rank, y = .data$GenesetID, yend = .data$statistic_fixed, color = .data$direction))
  
  ggp <-  ggplot(ggdata, aes(x = 0, y = .data$GenesetID, group = .data$GenesetID)) +
    geom_line() +
    geom_segment(aes(x = .data$rank, xend = .data$rank, y = .data$GenesetID_num - 0.45, yend = .data$GenesetID_num + 0.45, color = .data$statistic)) +
    ggtitle(title) +
    xlab(xlab) +
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      # axis.title.x = element_blank(), 
      # axis.text.x = element_blank(),
      # axis.ticks.x = element_blank(),
      plot.title = element_text(size = title_size, hjust = 0.5),
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
    ylab(statistic_name) +
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
  
  out1 <- cowplot::plot_grid(ggp, ggp2, nrow = 2, ncol = 1, rel_heights = c((length(genesets) + 4) / 5, 1), align = "v", axis = 'lr')
  
  
  
  if(!is.null(adjp) && !is.null(enrichment_score)){
    
    stopifnot(length(adjp) == length(genesets))
    stopifnot(length(enrichment_score) == length(genesets))
    
    
    gsea_results <- data.frame(geneset = names(genesets), adjp = adjp, enrichment_score = enrichment_score, stringsAsFactors = FALSE)
    
    colnames(gsea_results) <- c(geneset_name, adjp_name, enrichment_score_name)
    geneset_var <- geneset_name
    adjp_var <- adjp_name
    enrichment_score_var <- enrichment_score_name
    
    
    if(all(names(genesets) == gsea_results[, geneset_var])){
      
      ### Use a trick with the title. Otherwise, the plots are not aligned :/
      
      ggp3 <- wrapper_ora_dotplot_single(gsea_results, geneset_var = geneset_var, observed_var = NULL, adjp_var = adjp_var, color_point_var = enrichment_score_var, title = paste0(" ", paste0(rep("\n", stringr::str_count(title, "\n")), collapse = " ")," "), title_width = 0, title_size = title_size, size_range = c(2, 5)) +
        theme(axis.text.y = element_blank(), 
          axis.ticks.y = element_blank())
      
      
      out2 <- cowplot::plot_grid(ggp3, nrow = 2, ncol = 1, rel_heights = c((length(genesets) + 4) / 5, 1))
      
      
      out <- cowplot::plot_grid(out1, out2, nrow = 1, ncol = 2, rel_widths = c(4, 1), align = "h", axis = "tb")
      
      out
      
      
    }else{
      
      out1
      
    }
    
    
  }else{
    
    out1
    
  }
  
  
  
  
}



















#' Plot GSEA statistics ranks
#' 
#' @param x TopTable with statistic, for example, DGE results from limma.
#' @param gsea_results TopTable with selected GSEA results obtained by running 'wrapper_dispaly_significant_gsea'.
#' @export
wrapper_gsea_plot <- function(x, contrast, genesets, gene_var = "EntrezIDs", statistic_prefix = "t", sep = "_", 
  gsea_results = NULL, geneset_var = "GenesetID", adjp_var = "adj.P.Val", enrichment_score_var = "NES",
  trim_limits = 0.01,
  color_low = '#42399B', color_mid = "darkgrey", color_high = '#D70131',
  title = "", title_size = 10, title_width = 100, axis_text_y_size = 8, axis_text_y_width = 80){
  
  
  if(contrast == ""){
    sep <- ""
  }
  
  
  # -------------------------------------------------------------------------
  # Prepare data
  # -------------------------------------------------------------------------
  
  statistic <- x[, paste0(statistic_prefix, sep, contrast)]
  names(statistic) <- x[, gene_var]
  
  
  if(!is.null(gsea_results)){
    stopifnot(all(names(genesets) == gsea_results[, geneset_var]))
    adjp <- gsea_results[, adjp_var]
    enrichment_score <- gsea_results[, enrichment_score_var]
  }else{
    adjp <- NULL  
    enrichment_score <- NULL
  }
  
  
  
  wrapper_gsea_plot_core(statistic = statistic, genesets = genesets, adjp = adjp, enrichment_score = enrichment_score,
    statistic_name = statistic_prefix, geneset_name = geneset_var, adjp_name = adjp_var, enrichment_score_name = enrichment_score_var,
    trim_limits = trim_limits,
    color_low = color_low, color_mid = color_mid, color_high = color_high,
    title = title, title_size = title_size, title_width = title_width, axis_text_y_size = axis_text_y_size, axis_text_y_width = axis_text_y_width)
  
  
  
  
}





























