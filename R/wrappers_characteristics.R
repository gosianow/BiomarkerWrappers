





#' Table with distribution summary for a categorical covariate
#' 
#' @param data Data frame.
#' @export
wrapper_core_characteristics_cat <- function(data, covariate_var, strat_var = NULL, variable_names = NULL, caption = NULL, out_colname = "Value"){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(covariate_var, strat_var), drop = FALSE]
  
  stopifnot(length(covariate_var) == 1)
  stopifnot(is.factor(data[, covariate_var]))
  
  stopifnot(length(out_colname) == 1)
  
  if(!is.null(strat_var)){
    stopifnot(length(strat_var) == 1)
    stopifnot(is.factor(data[, strat_var]))
  }else{
    ### Add dummy variable to data
    stopifnot(!"strat_dummy" %in% colnames(data))
    data[, "strat_dummy"] <- factor(out_colname)
    strat_var <- "strat_dummy"
  }
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, covariate_var], data[, strat_var])
  
  prop <- prop.table(tbl, margin = 2) * 100
  
  tbl_isna <- table(factor(is.na(data[, covariate_var]), levels = c(FALSE, TRUE)), data[, strat_var])
  
  empty_row <- matrix(NA, nrow = 1, ncol = nlevels(data[, strat_var]))
  colnames(empty_row) <- levels(data[, strat_var])
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  countdf <- as.data.frame.matrix(tbl)
  
  propdf <- as.data.frame.matrix(prop)
  
  countm <- rbind.fill.matrix(empty_row, tbl_isna, countdf)
  colnames(countm) <- paste0("counts_", colnames(countm))
  
  stopifnot(ncol(countm) == ncol(tbl_isna))
  
  propm <- rbind.fill.matrix(empty_row, empty_row, empty_row, propdf)
  colnames(propm) <- paste0("proportions_", colnames(propm))
  
  stopifnot(ncol(propm) == ncol(tbl_isna))
  
  
  res <- data.frame(covariate = c(covariate_var, "Total (non-NA)", "NAs", rownames(tbl)),
    countm,
    propm,
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = c(variable_names[covariate_var], "Total (non-NA)", "NAs", rownames(tbl)), 
    
    format_counts_and_props(counts = countm, props = propm),
    
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  ### Do not display header when strat_var = NULL 
  
  if(strat_var == "strat_dummy"){
    header <- NULL
  }else{
    num_start_cols <- 1
    
    header <- c(num_start_cols, nlevels(data[, strat_var]))
    header <- as.integer(header)
    names(header) <- c(" ", variable_names[strat_var])
  }
  
  
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(!is.null(caption)){
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassCharacteristics(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}




#' Table with distribution summary for a numerical covariate
#' 
#' @param data Data frame.
#' @export
wrapper_core_characteristics_num <- function(data, covariate_var, strat_var = NULL, variable_names = NULL, caption = NULL, out_colname = "Value"){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(covariate_var, strat_var), drop = FALSE]
  
  stopifnot(length(covariate_var) == 1)
  stopifnot(is.numeric(data[, covariate_var]))
  
  stopifnot(length(out_colname) == 1)
  
  if(!is.null(strat_var)){
    stopifnot(length(strat_var) == 1)
    stopifnot(is.factor(data[, strat_var]))
  }else{
    ### Add dummy variable to data
    stopifnot(!"strat_dummy" %in% colnames(data))
    data[, "strat_dummy"] <- factor(out_colname)
    strat_var <- "strat_dummy"
  }
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics
  # --------------------------------------------------------------------------
  
  tbl_isna <- table(factor(is.na(data[, covariate_var]), levels = c(FALSE, TRUE)), data[, strat_var])
  
  empty_row <- matrix(NA, nrow = 1, ncol = nlevels(data[, strat_var]))
  colnames(empty_row) <- levels(data[, strat_var])
  
  
  Median = aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- t(data.frame(Median, Mean, Min, Max))
  colnames(summdf) <- levels(data[, strat_var])
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  summm <- rbind.fill.matrix(empty_row, tbl_isna, summdf)
  
  stopifnot(ncol(summm) == ncol(tbl_isna))
  
  
  res <- data.frame(covariate = c(covariate_var, "Total (non-NA)", "NAs", rownames(summdf)),
    summm,
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = c(variable_names[covariate_var], "Total (non-NA)", "NAs", rownames(summdf)), 
    
    format_summ(summ = summm),
    
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  ### Do not display header when strat_var = NULL 
  
  if(strat_var == "strat_dummy"){
    header <- NULL
  }else{
    num_start_cols <- 1
    
    header <- c(num_start_cols, nlevels(data[, strat_var]))
    header <- as.integer(header)
    names(header) <- c(" ", variable_names[strat_var])
  }
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(!is.null(caption)){
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassCharacteristics(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}




#' Table with distribution summary for categorical and numerical covariates
#' 
#' @param data Data frame.
#' @export
wrapper_core_characteristics <- function(data, covariate_vars, strat_var = NULL, variable_names = NULL, caption = NULL, out_colname = "Value"){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(covariate_vars, strat_var), drop = FALSE]
  
  
  stopifnot(length(covariate_vars) >= 1)
  
  vars_class <- sapply(data[, covariate_vars], class)
  stopifnot(all(vars_class %in% c("factor", "numeric", "integer")))
  
  stopifnot(length(out_colname) == 1)
  
  
  if(!is.null(strat_var)){
    stopifnot(length(strat_var) == 1)
    stopifnot(is.factor(data[, strat_var]))
  }else{
    ### Add dummy variable to data
    stopifnot(!"strat_dummy" %in% colnames(data))
    data[, "strat_dummy"] <- factor(out_colname)
    strat_var <- "strat_dummy"
  }
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate summaries
  # --------------------------------------------------------------------------
  
  
  wrapper_res <- lapply(1:length(covariate_vars), function(i){
    # i = 1
    
    covariate_var <- covariate_vars[i]
    
    
    if(class(data[, covariate_var]) == "factor"){
      
      wrapper_res <- wrapper_core_characteristics_cat(data = data, covariate_var = covariate_var, strat_var = strat_var, variable_names = variable_names, caption = caption, out_colname = out_colname)
      
    }else{
      
      wrapper_res <- wrapper_core_characteristics_num(data = data, covariate_var = covariate_var, strat_var = strat_var, variable_names = variable_names, caption = caption, out_colname = out_colname)
      
    }
    
    
    return(wrapper_res)
    
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(!is.null(caption)){
    
    ## Remove all underscores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  
  
  header <- bheader(wrapper_res[[1]])
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassCharacteristics(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}




#' Table with distribution summary for a list of covariates for ITT and BEP
#' 
#' @param data Data frame.
#' @param bep_vars Vector with column names for logical variables where TRUE indicates the biomarker evaluable population (BEP).
#' @export
wrapper_characteristics_bep <- function(data, covariate_vars, bep_vars = NULL, treatment_var = NULL, variable_names = NULL, caption = NULL, itt_name = "ITT"){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(covariate_vars, bep_vars, treatment_var), drop = FALSE]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  

  # --------------------------------------------------------------------------
  # Calculate characteristics for ITT
  # --------------------------------------------------------------------------
  
  
  characteristics_itt <- wrapper_core_characteristics(data = data, covariate_vars = covariate_vars, strat_var = NULL, variable_names = variable_names, caption = NULL, out_colname = itt_name)
  
  
  if(!is.null(treatment_var)){
    
    data$population_dummy <- factor(itt_name)
    
    ### Calculate population-treatment interaction
    data$population_treatment_interaction <- interaction(data$population_dummy, data[, treatment_var], drop = FALSE, lex.order = TRUE, sep = " : ")
    
    
    characteristics_itt_treatment <- wrapper_core_characteristics(data = data, covariate_vars = covariate_vars, strat_var = "population_treatment_interaction", variable_names = variable_names, caption = NULL)
    
    
    res <- cbind(bresults(characteristics_itt), bresults(characteristics_itt_treatment)[, -1, drop = FALSE])
    out <- cbind(boutput(characteristics_itt), boutput(characteristics_itt_treatment)[, -1, drop = FALSE])
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    characteristics_itt <- BclassCharacteristics(results = res, output = out)
    
  }
  
  
  # --------------------------------------------------------------------------
  # Calculate characteristics for BEPs
  # --------------------------------------------------------------------------
  
  
  if(!is.null(bep_vars)){
    
    characteristics_beps <- lapply(1:length(bep_vars), function(i){
      # i = 1
      
      data_bep <- data[data[, bep_vars[i]] %in% TRUE, , drop = FALSE]
      
      if(nrow(data_bep) == 0){
        return(NULL)
      }
      
      
      characteristics_bep <- wrapper_core_characteristics(data = data_bep, covariate_vars = covariate_vars, strat_var = NULL, variable_names = variable_names, caption = NULL, out_colname = variable_names[bep_vars[i]])
      
      
      if(!is.null(treatment_var)){
        
        data_bep$population_dummy <- factor(variable_names[bep_vars[i]])
        
        ### Calculate population-treatment interaction
        data_bep$population_treatment_interaction <- interaction(data_bep$population_dummy, data_bep[, treatment_var], drop = FALSE, lex.order = TRUE, sep = " : ")
        
        
        characteristics_bep_treatment <- wrapper_core_characteristics(data = data_bep, covariate_vars = covariate_vars, strat_var = "population_treatment_interaction", variable_names = variable_names, caption = NULL)
        
        
        res <- cbind(bresults(characteristics_bep), bresults(characteristics_bep_treatment)[, -1, drop = FALSE])
        out <- cbind(boutput(characteristics_bep), boutput(characteristics_bep_treatment)[, -1, drop = FALSE])
        
        rownames(res) <- NULL
        rownames(out) <- NULL
        
        characteristics_bep <- BclassCharacteristics(results = res, output = out)
        
      }
      
      
      characteristics_bep <- BclassCharacteristics(results = bresults(characteristics_bep)[, -1, drop = FALSE], output = boutput(characteristics_bep)[, -1, drop = FALSE])
      
      return(characteristics_bep)
      
    })
    
    
    res <- cbind(bresults(characteristics_itt), do.call("cbind", lapply(characteristics_beps, bresults)))
    out <- cbind(boutput(characteristics_itt), do.call("cbind", lapply(characteristics_beps, boutput)))
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    characteristics_itt <- BclassCharacteristics(results = res, output = out)
    
    
  }
  
  
  bcaption(characteristics_itt) <- caption
  
  
  return(characteristics_itt)
  
  
}




































