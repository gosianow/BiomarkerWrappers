



# data <- data_goya
# 
# x_var <- "Treatment_Arm"
# y_var <- "Ann_Arbor_Stage"
# 
# variable_names = NULL
# caption = NULL




#' Table with distribution summary for a categorical covariate
#' 
#' @param data Data frame.
wrapper_core_characteristics_cat <- function(data, x_var, y_var, variable_names = NULL, caption = NULL){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  stopifnot(nlevels(data[, x_var]) >= 2)
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.factor(data[, y_var]))
  stopifnot(nlevels(data[, y_var]) >= 2)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, y_var], data[, x_var])
  
  prop <- prop.table(tbl, margin = 2) * 100
  
  tbl_isna <- table(is.na(data[, y_var]), data[, x_var])
  
  empty_row <- matrix(NA, nrow = 1, ncol = nlevels(data[, x_var]))
  colnames(empty_row) <- levels(data[, x_var])
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  countdf <- as.data.frame.matrix(tbl)
  
  propdf <- as.data.frame.matrix(prop)
  
  countm <- rbind.fill.matrix(empty_row, tbl_isna, countdf)
  colnames(countm) <- paste0("counts_", colnames(countm))
  
  propm <- rbind.fill.matrix(empty_row, empty_row, empty_row, propdf)
  colnames(propm) <- paste0("proportions_", colnames(propm))
  
  
  res <- data.frame(category = c(y_var, "Total (non-NA)", "NAs", rownames(tbl)),
    countm,
    propm,
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Category = c(variable_names[y_var], "Total (non-NA)", "NAs", rownames(tbl)), 
    
    format_counts_and_props(counts = countm, props = propm),
    
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  num_start_cols <- 1
  
  
  header <- c(num_start_cols, nlevels(data[, x_var]))
  header <- as.integer(header)
  names(header) <- c(" ", variable_names[x_var])
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(!is.null(caption)){
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  

  
  bout <- BclassCharacteristics(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}




# data <- data_goya
# 
# x_var <- "Treatment_Arm"
# y_var <- "FCGR2B"
# 
# variable_names = NULL
# caption = NULL



#' Table with distribution summary for a numerical covariate
#' 
#' @param data Data frame.
wrapper_core_characteristics_num <- function(data, x_var, y_var, variable_names = NULL, caption = NULL){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  stopifnot(nlevels(data[, x_var]) >= 2)
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.numeric(data[, y_var]))
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics
  # --------------------------------------------------------------------------
  
  tbl_isna <- table(is.na(data[, y_var]), data[, x_var])
  
  empty_row <- matrix(NA, nrow = 1, ncol = nlevels(data[, x_var]))
  colnames(empty_row) <- levels(data[, x_var])
  
  
  Median = aggregate(data[, y_var], list(subgroup = data[, x_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = aggregate(data[, y_var], list(subgroup = data[, x_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = aggregate(data[, y_var], list(subgroup = data[, x_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = aggregate(data[, y_var], list(subgroup = data[, x_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- t(data.frame(Median, Mean, Min, Max))
  colnames(summdf) <- levels(data[, x_var])
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  summm <- rbind.fill.matrix(empty_row, tbl_isna, summdf)
  
  
  res <- data.frame(category = c(y_var, "Total (non-NA)", "NAs", rownames(summdf)),
    summm,
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Category = c(variable_names[y_var], "Total (non-NA)", "NAs", rownames(summdf)), 
    
    format_summ(summ = summm),
    
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  num_start_cols <- 1
  
  
  header <- c(num_start_cols, nlevels(data[, x_var]))
  header <- as.integer(header)
  names(header) <- c(" ", variable_names[x_var])
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(!is.null(caption)){
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  
  
  
  bout <- BclassCharacteristics(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}



# data <- data_goya
# 
# x_var <- "Treatment_Arm"
# y_vars <- c("FCGR2B", "Ann_Arbor_Stage")
# 
# variable_names = NULL
# caption = NULL



#' Table with distribution summary for categorical and numerical covariates
#' 
#' @param data Data frame.
wrapper_core_characteristics <- function(data, x_var, y_vars, variable_names = NULL, caption = NULL){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  stopifnot(nlevels(data[, x_var]) >= 2)
  
  stopifnot(length(y_vars) >= 1)
  
  vars_class <- sapply(data[, y_vars], class)
  stopifnot(all(vars_class %in% c("factor", "numeric", "integer")))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate summaries
  # --------------------------------------------------------------------------
  
  
  wrapper_res <- lapply(1:length(y_vars), function(i){
    # i = 1
    
    y_var <- y_vars[i]
    
    
    if(class(data[, y_var]) == "factor"){
      
      wrapper_res <- wrapper_core_characteristics_cat(data = data, x_var = x_var, y_var = y_var, variable_names = variable_names, caption = caption)
      
    }else{
      
      wrapper_res <- wrapper_core_characteristics_num(data = data, x_var = x_var, y_var = y_var, variable_names = variable_names, caption = caption)
      
    }
    
    
    return(wrapper_res)
    
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(!is.null(caption)){
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  
  
  header <- Bheader(wrapper_res[[1]])
  
  
  bout <- BclassCharacteristics(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}












