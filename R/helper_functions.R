


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



















