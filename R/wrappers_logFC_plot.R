



# gene_var = "Endpoint"
# lfc_prefix = "logFC"
# pval_prefix = "P.Value"
# adjp_prefix = "adj.P.Val"
# sep = "_"
# pval = 0.05
# title = ""
# axis.text.x.angle = 20
# axis.text.x.vjust = 1
# axis.text.x.hjust = 1
# axis.text.y.size = 10
# order_genes = FALSE
# color_low = '#D70131'
# color_high = '#42399B'
# radius_range = c(10, 3)




#' Dot plot with logFC and p-values for mutiple contrasts
#' 
#' @param x TopTable
wrapper_plot_logFC <- function(x, gene_var = "Hgnc_Symbol", lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", sep = "_", pval = 0.05, title = "", axis.text.x.angle = 30, axis.text.x.vjust = 1, axis.text.x.hjust = 1, axis.text.y.size = 10, order_genes = FALSE, color_low = '#D70131', color_high = '#42399B', radius_range = c(10, 3)){
  
  
  # axis.text.x.angle = 90, axis.text.x.vjust = 0.5, axis.text.x.hjust = 1
  # axis.text.x.angle = 45, axis.text.x.vjust = 1, axis.text.x.hjust = 1
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", lfc_prefix, sep), "", grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE))
  
  
  if(order_genes){
    for(i in length(contrasts):1){
      # i = 1
      x <- x[order(x[, paste0(lfc_prefix, sep, contrasts[i])], decreasing = TRUE), , drop = FALSE]
    }
  }
  
  
  cols_lfc <- grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE)
  data_lfc <- pivot_longer(x[, c(gene_var, cols_lfc), drop = FALSE], cols = cols_lfc, names_to = "contrast", values_to = lfc_prefix)
  data_lfc$contrast <- gsub(paste0("^", lfc_prefix, sep), "", data_lfc$contrast)
  
  
  cols_pval <- grep(paste0("^", pval_prefix, sep), colnames(x), value = TRUE)
  data_pval <- pivot_longer(x[, c(gene_var, cols_pval), drop = FALSE], cols = cols_pval, names_to = "contrast", values_to = pval_prefix)
  data_pval$contrast <- gsub(paste0("^", pval_prefix, sep), "", data_pval$contrast)
  
  
  cols_adjp <- grep(paste0("^", adjp_prefix, sep), colnames(x), value = TRUE)
  data_adjp <- pivot_longer(x[, c(gene_var, cols_adjp), drop = FALSE], cols = cols_adjp, names_to = "contrast", values_to = adjp_prefix)
  data_adjp$contrast <- gsub(paste0("^", adjp_prefix, sep), "", data_adjp$contrast)
  
  
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
  
  
  
  
  ggp <- ggplot(data, aes_string(x = "contrast", y = gene_var, size = "pval_cut", color = lfc_prefix)) +
    geom_point() +
    geom_point(aes_string(size = "pval_cut", shape = "significance"), color = "black", show.legend = TRUE) +
    ggtitle(title) +
    theme_cowplot(12) +
    theme(axis.line = element_blank(), 
      axis.title = element_blank(), 
      axis.text.x = element_text(angle = axis.text.x.angle, vjust = axis.text.x.vjust, hjust = axis.text.x.hjust),
      axis.text.y = element_text(size = axis.text.y.size)) +
    panel_border(colour = "black", linetype = 1, size = 1, remove = FALSE) +
    background_grid(major = "xy", minor = "none", size.major = 0.25) +
    scale_radius(name = pval_prefix, range = radius_range, breaks = 1:(length(pval_cut) - 1), labels = pval_cut[-1], limits = c(1, length(pval_cut) - 1)) +
    scale_colour_gradient2(name = lfc_prefix, low = color_low, high = color_high) + 
    scale_shape_manual(name = adjp_prefix, values = values_shape)
  
  
  # This is great. squish in this context converts clamps all values to be within the min and max of the limits argument. i.e., if value < min(limits) then value = min(limits) else if value > max(limits) then value = max(limits).
  # scale_colour_gradient2(limits = c(-1.5, 1.5), oob = scales::squish)
  
  
  return(ggp)
  
  
  
}












































