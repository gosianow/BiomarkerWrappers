


#' Create variable names
#' 
#' Creates from scratch or adds missing variable names for variables in data.
#' 
#' @param data Data frame.
#' @param variable_names Named vector of variable names corresponding to variables in data. This vector does not have to contain names for all the variables in data. If names for some variables are missing, they will be created. If NULL, variable names are created by subtracting underscore from the column names of data.
format_variable_names <- function(data, variable_names = NULL){
  
  
  new_variable_names <- gsub("_", " ", colnames(data))
  names(new_variable_names) <- colnames(data)
  
  if(!is.null(variable_names)){
    
    stioifnot(!is.null(names(variable_names)))
    
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




#' Check or create colors
#' 
#' Check or create colors. The returned vector must be named.
#' 
#' @param levels Vector of factor levels for which we want to specify colors.
#' @param colors Vector of colors of length of levels. If NULL, colors are created.
format_colors <- function(levels, colors = NULL){
  
  if(is.null(colors)){
    
    colors <- colorRampPalette(brewer.pal(12, "Paired"))(12)[1:length(levels)]
    names(colors) <- levels
    
    
  }else{
    
    stopifnot(length(colors) == length(levels))
    stopifnot(names(colors) == levels)
    
  }
  
  
  return(colors)
  
  
}



















