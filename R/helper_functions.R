






indicate_blocks <- function(d, block_vars, return = "seq"){
  
  
  all_data <- apply(d[, block_vars, drop = FALSE], 1, paste, collapse = ".")
  
  all_rle <- c(0, cumsum(rle(all_data)$lengths))
  
  
  if(length(all_rle) > 2){
    
    
    if(return == "seq"){
      
      max_rle_indx <- ifelse(length(all_rle) %% 2 == 0, length(all_rle), length(all_rle) - 1)
      
      out <- unlist(lapply(seq(1, max_rle_indx, by = 2), function(i){
        seq(all_rle[i] + 1, all_rle[i + 1], by = 1)
      }))
      
      
    }else if(return == "line"){
      
      out <- all_rle[-1]
      
    }
    
    
    
  }else{
    out <- NULL
  }
  
  
  return(out)
  
}













#' Format p-values
#' 
#' @param x Vector of p-values to be formatted.
#' @param digits Number of digits after decimial to display.
#' @param asterisk Signif. codes:  0 `***` 0.001 `**` 0.01 `*` 0.05.
format_pvalues <- function(x, digits = 4, asterisk = TRUE){
  
  # digits = 4
  # asterisk <- TRUE
  # x <- c(0.2, 0.05, 0.034534, 1.366332e-05, 1.366332e-04, NA)
  
  
  min_pval <- 1/10^digits
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = TRUE)
  output[x < min_pval] <- paste0("<", formatC(min_pval, format = "f", digits = digits))
  output[is.na(x)] <- ""
  
  output
  
  
  if(asterisk){
    
    pval_asterisk <- ifelse(x < 0.001, " ***", ifelse(x < 0.01, " **", ifelse(x < 0.05, " *", "")))
    pval_asterisk[is.na(pval_asterisk)] <- ""
    
    output <- paste0(output, pval_asterisk)
    
  }
  
  output
  
  
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



















