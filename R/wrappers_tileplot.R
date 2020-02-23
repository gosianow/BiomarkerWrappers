
#' Tile plot showing co-occurrence and mutual exclusivity of variant types for selected genes
#' 
#' @param data_variants_spread Data frame where each row corresponds to one sample.
#' @param var_samples Column name defining samples.
#' @param vars_cat Column names of categorical variables defining variants for genes of interest. Order of the variables defined order in which they will be plotted starting from the top of the figure.
#' @param names_cat Names corresponding to the categorical variables that should be displayed on the plot, for example, gene names.
#' @param vars_num Column names of numerical variables, for example, one corresponding to TMB.
#' @param names_num Names corresponding to the numerical variables that should be displayed on the plot.
#' @param cut_points List of the same length as vars_num specifying vectors of cutting points for each of the numerical variables in vars_num.
#' @param cut_colors List of the same length as vars_num specifying vectors of colors for the intervals obtained by cutting data using cut_points.
#' @param levels_variants Vector defining the order in which variants should be plotted.
#' @param colors_variants Vector defining colors for the variants.
#' @param title String defining the plot title.
wrapper_tileplot <- function(data_variants_spread, var_samples, vars_cat, names_cat, vars_num = NULL, names_num = NULL, cut_points = NULL, cut_colors = NULL, levels_variants, colors_variants, title = ""){
  
  ## Order samples by values in vars_num
  if(!is.null(vars_num)){
    for(i in length(vars_num):1){
      oo <- order(data_variants_spread[, vars_num[i]], decreasing = TRUE)
      data_variants_spread <- data_variants_spread[oo, ]
    }
  }
  
  ## Order samples by variants in vars_cat
  for(i in length(vars_cat):1){
    data_variants_spread[, vars_cat[i]] <- factor(data_variants_spread[, vars_cat[i]], levels = levels_variants)
    oo <- order(data_variants_spread[, vars_cat[i]], decreasing = TRUE)
    data_variants_spread <- data_variants_spread[oo, ]
  }
  
  levels_samples <- data_variants_spread[, var_samples]
  
  ### Generate a tile plot for categorical variables
  
  ## Transform data into the long format for plotting with ggplot
  data_sub <- data_variants_spread[, c(var_samples, vars_cat)]
  colnames(data_sub) <- c(var_samples, names_cat)
  data_variants_gather <- gather(data_sub, key = "gene", value = "variant_type", names_cat)
  
  data_variants_gather$gene <- factor(data_variants_gather$gene, levels = rev(names_cat))
  data_variants_gather[, var_samples] <- factor(data_variants_gather[, var_samples], levels = levels_samples)
  data_variants_gather$variant_type <- factor(data_variants_gather$variant_type, levels = levels_variants)
  

  ## Generate the tile plot
  ggp <- ggplot(data_variants_gather, aes_string(x = var_samples, y = "gene", fill = "variant_type")) +
    geom_tile(height = 0.9) +
    ggtitle(title) +
    labs(caption = paste0("Number of samples: ", nrow(data_variants_spread))) +
    theme_classic() + 
    theme(legend.position = "bottom",
      axis.title = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_text(color = "black", size = 14),
      axis.ticks = element_blank(),
      axis.line = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA)) +
    scale_fill_manual(name = "", values = colors_variants) +
    scale_x_discrete(expand = c(0,0))
  
  
  ### Generate bar plots for numerical variables
  
  if(!is.null(vars_num)){
    
    ## Generate the bar plots
    ggpl <- list()
    
    for(i in 1:length(vars_num)){
      # i = 1
      
      data_sub <- data_variants_spread[, c(var_samples, vars_num[i])]
      data_sub[, var_samples] <- factor(data_sub[, var_samples], levels = levels_samples)
      ## Cut the numeric data
      data_sub[, paste0(vars_num[i], "_cut")] <- cut(data_sub[, vars_num[i]], breaks = cut_points[[i]], include.lowest = TRUE, right = FALSE, ordered_result = TRUE)
      
      
      ggpl[[i]] <- ggplot(data_sub, aes_string(x = var_samples, y = vars_num[i], fill = paste0(vars_num[i], "_cut"))) +
        geom_col() +
        ggtitle(title) +
        ylab(names_num[i]) +
        theme_classic() + 
        theme(legend.position = "top",
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(color = "black", size = 8),
          axis.ticks.x = element_blank(),
          axis.line = element_blank(),
          panel.border = element_rect(colour = "black", fill = NA)) +
        scale_x_discrete(expand = c(0,0)) +
        scale_fill_manual(name = names_num[i], values = cut_colors[[i]])
      
    }
    
    ggpl[[length(vars_num) + 1]] <- ggp
    
    ggp <- plot_grid(plotlist = ggpl, ncol = 1, align = "v", axis = "rlbt", rel_heights = c(rep(2, length(vars_num)), length(vars_cat)))
    
    
  }
  
  
  return(ggp)
  
}








