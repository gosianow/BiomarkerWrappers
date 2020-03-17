




# x <- topTable
# genesets; contrast; gene_var = "EntrezIDs"; statistic_prefix = "t"; sep = "_"; 
# trim_limits = 0.02; min_GS_size = 10; max_GS_size = 500;
# title = ""; title_size = 10; title_width = 100; axis_text_y_size = 8; axis_text_y_width = 70




#' Plot GSEA statistics ranks
#' 
#' @param x TopTable.
wrapper_plot_GSEA <- function(x, genesets, contrast, gene_var = "EntrezIDs", statistic_prefix = "t", sep = "_", 
  trim_limits = 0.02, min_GS_size = 10, max_GS_size = 500,
  title = "", title_size = 10, title_width = 100, axis_text_y_size = 8, axis_text_y_width = 70){
  
  
  ## Wrap the title so it can be nicely displayed in the plots
  title <- stringr::str_wrap(title, width = title_width)
  
  
  # -------------------------------------------------------------------------
  # Prepare data
  # -------------------------------------------------------------------------
  
  statistic <- topTable[, paste0(statistic_prefix, sep, contrast)]
  names(statistic) <- topTable[, gene_var]
  
  statistic <- sort(statistic, decreasing = TRUE)
  
  
  data <- data.frame(universe = names(statistic), statistic = statistic, rank = rank(-statistic), 
    stringsAsFactors = FALSE)
  
  data$direction <- factor(ifelse(data$statistic > 0, "up", "down"))
  
  
  
  ### Trim limits
  
  if(is.null(trim_limits)){
    max_abs_value <- ceiling(max(abs(range(data[, "statistic"], na.rm = TRUE))))
  }else if(trim_limits >= 1){
    max_abs_value <- trim_limits
  }else{
    ### Use quantiles 
    max_abs_value <- ceiling(max(abs(quantile(data[, "statistic"], probs = c(trim_limits, 1 - trim_limits), na.rm = TRUE))))
  }
  
  
  data$statistic[data$statistic > max_abs_value] <- max_abs_value
  data$statistic[data$statistic < -max_abs_value] <- -max_abs_value
  
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
  
  
  ### Exclude too small or too large gene sets
  
  size_genesets <- sapply(genesets, length)
  
  genesets <- genesets[size_genesets >= min_GS_size & size_genesets <= max_GS_size]
  
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
  
  ggdata$statistic_adj2 <- ggdata$statistic_adj + ggdata$Geneset_num + ifelse(ggdata$direction == "up", 0.2, -0.2)
  
  ggdata$statistic_fixed <- ggdata$Geneset_num + ifelse(ggdata$direction == "up", 0.5, -0.5)
    
    
  # -------------------------------------------------------------------------
  # ggplot
  # -------------------------------------------------------------------------
  
  
  xlim <- c(0, max(data$rank))
  xlab <- "Rank"
  
  
  
  
  ggp <-  ggplot(ggdata, aes(x = rank, xend = rank, y = Geneset, yend = statistic_adj2, color = direction)) +
    geom_segment() +
    ggtitle(title) +
    xlab(xlab) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      axis.title.y = element_blank(), 
      axis.text.y = element_text(size = axis_text_y_size),
      plot.title = element_text(size = title_size, hjust = 1),
      legend.position = "none") +
    coord_cartesian(xlim = xlim) +
    panel_border(colour = "black", linetype = 1, size = 1, remove = FALSE) +
    background_grid(major = "y", minor = "none", size.major = 0.25) +
    scale_color_manual(values = c('#42399B', '#D70131')) +
    scale_x_continuous(expand = c(0, 0))
  
  
  ggp
  
  
  
  
}






























