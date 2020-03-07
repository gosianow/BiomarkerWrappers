










#' Calculate logFC for expression that will be plotted in a heatmap 
#' 
#' @param x eSet object
wrapper_calculate_sample_logFC <- function(x, comparison_var, subgroup_var = NULL){
  
  
  if(is.null(subgroup_var)){
    pData(x)$dummy_subgroup_var <- factor("dummy_subgroup_var")
    subgroup_var <- "dummy_subgroup_var"
  }
  
  subgroup_levels <- levels(pData(x)[, subgroup_var])
  
  comparison_levels <- levels(pData(x)[, comparison_var])
  
  reference_level <- comparison_levels[1]
  
  tbl <- table(pData(x)[, comparison_var])
  
  stopifnot(tbl[reference_level] > 0)
  
  
  expr_heatmap <- lapply(1:length(subgroup_levels), function(i){
    # i = 1
    
    x_sub <- x[, pData(x)[, subgroup_var] == subgroup_levels[i]]
    
    expr_sub <- Biobase::exprs(x_sub)
    
    expr_lfc <- expr_sub - rowMeans(expr_sub[, pData(x_sub)[, comparison_var] == reference_level], na.rm = TRUE)
    
    expr_lfc
    
    
  }) 
  
  
  expr_heatmap <- do.call("cbind", expr_heatmap)
  
  
  ### Bring the same sample order as in x
  
  expr_heatmap <- expr_heatmap[, colnames(Biobase::exprs(x)), drop = FALSE]
  
  ### Return an eSet
  
  out <- ExpressionSet(assayData = expr_heatmap, 
    phenoData = phenoData(x), 
    featureData = featureData(x))
  

  out
  
}




duplicated2 <- function(x, value = FALSE, na.rm = FALSE){
  
  if(is.factor(x)){
    x <- as.character(x)
  }
  
  ### Make the NAs to count as unique
  if(na.rm){
    x[is.na(x)] <- paste0("X...NA...", 1:sum(is.na(x)))
  }
  
  ifduplicated <- duplicated(x, fromLast = FALSE) | duplicated(x, fromLast = TRUE)
  
  if(value){
    out <- x[ifduplicated]
  }else{
    out <- ifduplicated
  }
  
  return(out)
  
}



length_nonNA <- function(x){
  sum(!is.na(x))
}




#' Read gmt file and return a list of genes
wrapper_read_gmt <- function(filename){
  x <- GSEABase::getGmt(filename)
  x <- GSEABase::geneIds(x)
  return(x)
}


#' Write a list of genes into a gmt file
wrapper_write_gmt <- function(x, filename){
  
  gsc <- GSEABase::GeneSetCollection(lapply(1:length(x), function(i){
    # i = 1
    GSEABase::GeneSet(setName = names(x[i]), geneIds = x[[i]])
    
  }))
  
  GSEABase::toGmt(gsc, filename)
  
  invisible(NULL)
  
}











wrapper_print_plot_grid <- function(plotlist, nsplit = 2, ncol = 2, nrow = 1){
  
  indx <- 1:length(plotlist)
  indx_split <- split(indx, ceiling(seq_along(indx) / nsplit))
  
  for(i in 1:length(indx_split)){
    print(plot_grid(plotlist = plotlist[indx_split[[i]]], ncol = ncol, nrow = nrow))
  }
  
  invisible(NULL)
  
}






# x <- data$FCGR1A


cut_core_quartiles <- function(x, labels = c("(0%, 25%]", "(25%, 50%]", "(50%, 75%]", "(75%, 100%]")){
  
  stopifnot(length(labels) == 4)
  
  out <- ggplot2::cut_number(x, n = 4)
  
  out <- factor(out, labels = labels)
  
  return(out)
  
}


cut_core_median <- function(x, labels = c("<=MED", ">MED")){
  
  stopifnot(length(labels) == 2)
  
  out <- ggplot2::cut_number(x, n = 2)
  
  out <- factor(out, labels = labels)
  
  return(out)
  
}


cut_core_2groups <- function(x, probs = 0.5, cutoff = NULL, labels = c("low", "high")){
  
  stopifnot(length(probs) == 1)
  stopifnot(length(labels) == 2)
  
  if(is.null(cutoff)){
    cutoff <- quantile(x, probs = probs, na.rm = TRUE)
  }
  
  stopifnot(length(cutoff) == 1)
  
  out <- factor(ifelse(x <= cutoff, labels[1], labels[2]), levels = labels)
  
  return(out)
  
}



cut_core_quartiles_strat <- function(x, strata, labels = c("(0%, 25%]", "(25%, 50%]", "(50%, 75%]", "(75%, 100%]")){
  
  stopifnot(is.factor(strata))
  
  strata_levels <- levels(strata)
  
  out <- factor(rep(NA, length(x)), levels = labels)
  
  for(i in 1:length(strata_levels)){
    # i = 1
    
    indx_sub <- which(strata == strata_levels[i])
    
    out[indx_sub] <- cut_core_quartiles(x[indx_sub], labels = labels)
    
  }
  
  return(out)
  
}



cut_core_median_strat <- function(x, strata, labels = c("<=MED", ">MED")){
  
  stopifnot(is.factor(strata))
  
  strata_levels <- levels(strata)
  
  out <- factor(rep(NA, length(x)), levels = labels)
  
  for(i in 1:length(strata_levels)){
    # i = 1
    
    indx_sub <- which(strata == strata_levels[i])
    
    out[indx_sub] <- cut_core_median(x[indx_sub], labels = labels)
    
  }
  
  return(out)
  
}



cut_core_2groups_strat <- function(x, strata, probs = 0.5, cutoff = NULL, labels = c("low", "high")){
  
  stopifnot(is.factor(strata))
  
  strata_levels <- levels(strata)
  
  out <- factor(rep(NA, length(x)), levels = labels)
  
  for(i in 1:length(strata_levels)){
    # i = 1
    
    indx_sub <- which(strata == strata_levels[i])
    
    out[indx_sub] <- cut_core_2groups(x[indx_sub], probs = probs, cutoff = cutoff, labels = labels) 
    
  }
  
  return(out)
  
}










indicate_blocks <- function(d, block_vars, return = "block"){
  
  stopifnot(length(return) == 1)
  stopifnot(return %in% c("block", "line", "empty"))
  
  
  all_data <- apply(d[, block_vars, drop = FALSE], 1, paste, collapse = ".")
  
  all_rle <- c(0, cumsum(rle(all_data)$lengths))
  
  
  if(return == "block"){
    
    if(length(all_rle) > 2){
      
      max_rle_indx <- ifelse(length(all_rle) %% 2 == 0, length(all_rle), length(all_rle) - 1)
      
      out <- unlist(lapply(seq(1, max_rle_indx, by = 2), function(i){
        seq(all_rle[i] + 1, all_rle[i + 1], by = 1)
      }))
      
    }else{
      out <- NULL
    }
    
    
  }else if(return == "line"){
    
    out <- all_rle[-1]
    
    
  }else if(return == "empty"){
    
    ## Indicate first lines that are unique
    out <- c(1, all_rle[-1] + 1)
    out <- out[-length(out)]
    
    ## But make empty all besides them
    
    out <- setdiff(1:nrow(d), out)
    
    
  }
  
  
  return(out)
  
}













#' Format p-values
#' 
#' @param x Vector of p-values to be formatted.
#' @param digits Number of digits after decimial to display.
#' @param asterisk Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05 `.` 0.1.
format_pvalues <- function(x, digits = 4, asterisk = TRUE){
  
  # digits = 4
  # asterisk <- TRUE
  # x <- c(0.2, 0.05, 0.034534, 1.366332e-05, 1.366332e-04, NA)
  
  if(sum(is.na(x)) == length(x)){
    return(rep("", length(x)))
  }
  
  min_pval <- 1/10^digits
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = TRUE)
  output[x < min_pval] <- paste0("<", formatC(min_pval, format = "f", digits = digits))
  output[is.na(x)] <- ""
  
  
  if(asterisk){
    
    pval_asterisk <- ifelse(x < 0.001, " ***", ifelse(x < 0.01, " **", ifelse(x < 0.05, " *", ifelse(x < 0.1, " .", ""))))
    
    pval_asterisk[is.na(pval_asterisk)] <- ""
    
    output <- paste0(output, pval_asterisk)
    
  }
  
  
  return(output)
  
}




#' Format Odds Ratios
#' 
#' @param x Vector of levels.
#' @param digits Number of decimal places.
format_or <- function(x, digits = 2){
  
  if(sum(is.na(x)) == length(x)){
    return(rep("", length(x)))
  }
  
  min_val <- 1/10^digits
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = TRUE)
  output[x < min_val] <- paste0("<", formatC(min_val, format = "f", digits = digits))
  output[is.na(x)] <- ""
  
  return(output)
  
}


#' Format difference 
#' 
#' @param x Vector of levels.
#' @param digits Number of decimal places.
format_difference <- function(x, digits = 2){
  
  if(sum(is.na(x)) == length(x)){
    return(rep("", length(x)))
  }
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = TRUE)
  output[is.na(x)] <- ""
  
  return(output)
  
}



#' Format CIs (confidence intervals)
#' 
#' @param CI_lower Vector of lower CIs.
#' @param CI_upper Vector of upper CIs.
format_CIs <- function(CI_lower, CI_upper, digits = 2){
  
  output <- paste0("(", round(CI_lower, digits), "-", round(CI_upper, digits), ")")
  
  return(output)
  
}



#' Format versus
#' 
#' @param level Vector of levels.
#' @param reference Vector of references.
format_vs <- function(level, reference){
  
  output <- paste0(level, " vs ", reference)
  
  output[level == "" & reference == ""] <- ""
  output[level == 0 & reference == 0] <- ""
  
  return(output)
  
}


#' Paste counts and proportions corresponding to one subgroup
#' 
#' @param counts Data frame with counts.
#' @param props Data frame with proportions.
#' @param digits Number of decimal places when rounding proportions.
format_counts_and_props <- function(counts, props, digits = 2, prefix_counts = "counts_"){
  
  ### Match parrern at the beginning
  pattern <- paste0("^", prefix_counts)
  
  output_names <- gsub(pattern, "", colnames(counts))
  
  stopifnot(all(dim(counts) == dim(props)))
  
  
  output <- lapply(1:nrow(counts), function(i){
    # i = 1
    
    out <- paste0(ifelse(is.na(counts[i, ]), "", counts[i, ]), ifelse(is.na(props[i, ]), "", paste0(" (", round(props[i, ], digits), "%)")))
    
    # Remove white spaces from the beginning and the end of a string
    
    out <- stringr::str_trim(out, side = "both")
    
    return(out)
    
  })
  
  output <- data.frame(do.call("rbind", output), stringsAsFactors = FALSE)
  
  colnames(output) <- output_names
  
  return(output)
  
}


#' Format summary statistics such as N, mean, median, min, max
#' 
#' @param summ Data frame.
format_summ <- function(summ, per = "row", digits = 2){
  
  stopifnot(per %in% c("row", "col"))
  
  if(per == "row"){
    
    output <- lapply(1:nrow(summ), function(i){
      # i = 1
      
      out <- ifelse(is.na(summ[i, ]), "", formatC(round(summ[i, ], digits), format = "f", digits = digits, drop0trailing = TRUE))
      
      
      return(out)
      
    })
    
    output <- data.frame(do.call("rbind", output), stringsAsFactors = FALSE)
    
  }else{
    
    output <- lapply(1:ncol(summ), function(i){
      # i = 1
      
      out <- ifelse(is.na(summ[, i]), "", formatC(round(summ[, i], digits), format = "f", digits = digits, drop0trailing = TRUE))
      
      
      return(out)
      
    })
    
    output <- data.frame(do.call("cbind", output), stringsAsFactors = FALSE)
    
  }
  
  
  colnames(output) <- colnames(summ)
  
  return(output)
  
}



#' Create variable names
#' 
#' Make sure that unique variable names exists for all the variables in the data frame. Creates from scratch or adds missing variable names for variables in data.
#' 
#' @param data Data frame.
#' @param variable_names Named vector of variable names corresponding to variables in data. This vector does not have to contain names for all the variables in data. If names for some variables are missing, they will be created. If NULL, variable names are created by subtracting underscore from the column names of data.
#' @return Named vector of unique variable names for all variables from data.
format_variable_names <- function(data, variable_names = NULL){
  
  
  new_variable_names <- gsub("_", " ", colnames(data))
  names(new_variable_names) <- colnames(data)
  
  if(!is.null(variable_names)){
    
    stopifnot(!is.null(names(variable_names)))
    
    mm <- match(names(new_variable_names), names(variable_names))
    new_variable_names[!is.na(mm)] <- variable_names[na.omit(mm)]
    
  }
  
  ### There cannot be duplicated names. If so, fix it.
  
  if(sum(duplicated(new_variable_names)) > 0){
    
    dupl_names <- unique(new_variable_names[duplicated(new_variable_names)])
    
    fixed_names <- lapply(1:length(dupl_names), function(i){
      # i = 1
      
      dupl_indx <- which(new_variable_names == dupl_names[i])
      
      new_names <- paste(dupl_names[i], 1:length(dupl_indx))
      
      out <- data.frame(dupl_indx = dupl_indx, new_names = new_names)
      
      
    })
    
    fixed_names <- plyr::rbind.fill(fixed_names)
    
    new_variable_names[fixed_names[, "dupl_indx"]] <- fixed_names[, "new_names"]
    
  }
  
  
  return(new_variable_names)
  
  
}


# levels <- levels(data_goya[, "Ann_Arbor_Stage"])


#' Check or create colors
#' 
#' Make sure that unique named colors are created for each level.
#' 
#' @param levels Vector of factor levels for which we want to specify colors.
#' @param colors Vector of colors longer or equal the number of levels. Can be named or non-named. If NULL, colors are created.
#' @param palette Vector of at least two colors used to create a color palette. 
#' @return Named vector of unique colors for all levels.
format_colors <- function(levels, colors = NULL, palette = NULL){
  
  if(is.null(colors)){
    
    if(is.null(palette)){
      
      stopifnot(length(levels) <= 12)
      
      colors <- RColorBrewer::brewer.pal(12, "Paired")[1:length(levels)]
      names(colors) <- levels
      
      
    }else{
      
      stopifnot(length(palette) >= 2)
      
      colors <- colorRampPalette(palette)(length(levels))
      names(colors) <- levels
      
      
    }
    
    
    # barplot(rep(1, length(colors)), col = colors)
    
    
  }else{
    
    stopifnot(length(colors) >= length(levels))
    
    if(is.null(names(colors))){
      colors <- colors[1:length(levels)]
      names(colors) <- levels
    }else{
      stopifnot(all(levels %in% names(colors)))
      colors <- colors[levels]
    }
    
    ### Colors have to be unique. Otherwise, ggsurvplot does not work.
    stopifnot(sum(duplicated(colors)) == 0)
    
  }
  
  
  return(colors)
  
  
}



calculate_break_time <- function(x, n_breaks = 10){
  
  
  max_tte_tmp <- max(x, na.rm = TRUE) / n_breaks
  
  decimial_nr <- round(log10(max_tte_tmp))
  
  break.time.by <- round(max_tte_tmp / 10^decimial_nr) * 10^decimial_nr
  
  return(break.time.by)
  
}















