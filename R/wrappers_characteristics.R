





#' Table with distribution summary for a categorical covariate
#' 
#' @param data Data frame.
#' @export
wrapper_characteristics_core_cat <- function(data, covariate_var, strat_var = NULL, variable_names = NULL, caption = NULL, out_colname = "Value"){
  
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
  
  countm <- plyr::rbind.fill.matrix(empty_row, tbl_isna, countdf)
  colnames(countm) <- paste0("counts_", colnames(countm))
  
  stopifnot(ncol(countm) == ncol(tbl_isna))
  
  propm <- plyr::rbind.fill.matrix(empty_row, empty_row, empty_row, propdf)
  colnames(propm) <- paste0("proportions_", colnames(propm))
  
  stopifnot(ncol(propm) == ncol(tbl_isna))
  
  
  res <- data.frame(covariate = c(covariate_var, "N", "NAs", rownames(tbl)),
    countm,
    propm,
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = c(variable_names[covariate_var], "N", "NAs", rownames(tbl)), 
    
    format_counts_and_props_df(counts = countm, props = propm),
    
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
wrapper_characteristics_core_num <- function(data, covariate_var, strat_var = NULL, variable_names = NULL, caption = NULL, out_colname = "Value", display_statistics = c("Median", "Mean")){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  stopifnot(length(display_statistics) >= 1)
  stopifnot(display_statistics %in% c("Median", "Mean", "Min", "Max", "First.Quartile", "Third.Quartile"))
  
  
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
  
  
  Median = stats::aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = stats::aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = stats::aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = stats::aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  First.Quartile = stats::aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = quantile, probs = 0.25, na.rm = TRUE, drop = FALSE)[, 2]
  Third.Quartile = stats::aggregate(data[, covariate_var], list(subgroup = data[, strat_var]), FUN = quantile, probs = 0.75, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- data.frame(Median, Mean, Min, Max, First.Quartile, Third.Quartile)
  summdf <- summdf[, display_statistics, drop = FALSE]
  summdf <- t(summdf)
  colnames(summdf) <- levels(data[, strat_var])
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  summm <- plyr::rbind.fill.matrix(empty_row, tbl_isna, summdf)
  
  stopifnot(ncol(summm) == ncol(tbl_isna))
  
  
  res <- data.frame(covariate = c(covariate_var, "N", "NAs", rownames(summdf)),
    summm,
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = c(variable_names[covariate_var], "N", "NAs", rownames(summdf)), 
    
    format_summ_df(summ = summm, digits = c(rep(0, 3), rep(2, nrow(summdf))), per = "col"),
    
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
wrapper_characteristics_core <- function(data, covariate_vars, strat_var = NULL, variable_names = NULL, caption = NULL, out_colname = "Value", display_statistics = c("Median", "Mean")){
  
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
      
      wrapper_res <- wrapper_characteristics_core_cat(data = data, covariate_var = covariate_var, strat_var = strat_var, variable_names = variable_names, caption = caption, out_colname = out_colname)
      
    }else{
      
      wrapper_res <- wrapper_characteristics_core_num(data = data, covariate_var = covariate_var, strat_var = strat_var, variable_names = variable_names, caption = caption, out_colname = out_colname, display_statistics = display_statistics)
      
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
#' @param covariate_vars Covariates to summarise
#' @param bep_vars Vector with column names for logical variables where TRUE indicates the biomarker evaluable population (BEP).
#' @export
wrapper_characteristics_bep <- function(data, covariate_vars, bep_vars = NULL, treatment_var = NULL, population_var = "Population", strat_vars = c(population_var, treatment_var), strat1_var = NULL, variable_names = NULL, caption = NULL, itt_name = "ITT", display_statistics = c("Median", "Mean"), lex_order = TRUE, include_pooled_arms = TRUE){
  
  
  
  ### Keep only those variables that are used for the analysis
  
  data <- data[, c(covariate_vars, bep_vars, treatment_var, strat1_var), drop = FALSE]
  
  
  data_list <- list()
  
  data_itt <- data
  data_itt[, population_var] <- itt_name
  
  
  data_list[[itt_name]] <- data_itt
  
  for(i in seq_along(bep_vars)){
    # i = 1
    
    data_bep <- data[data[, bep_vars[i]], ]
    data_bep[, population_var] <- bep_vars[i]
    data_list[[bep_vars[i]]] <- data_bep
    
  }
  
  
  data_rbind <- plyr::rbind.fill(data_list)
  data_rbind[, population_var] <- factor(data_rbind[, population_var], levels = c(itt_name, bep_vars))
  

  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  
  variable_names <- format_variable_names(data = data_rbind, variable_names = variable_names)
  
  
  
  # --------------------------------------------------------------------------
  # Calculate characteristics for pooled arms
  # --------------------------------------------------------------------------
  
  
  strat_var <- paste0(strat_vars, collapse = " : ")
  
  
  if(!is.null(treatment_var) && treatment_var %in% strat_vars && include_pooled_arms){
    
    data_rbind_pooled <- data_rbind
    data_rbind_pooled[, treatment_var] <- ""
    
    levels_treatment <- levels(data_rbind[, treatment_var])
    
    data_rbind[, treatment_var] <- as.character(data_rbind[, treatment_var])
    
    
    data_rbind <- rbind.fill(data_rbind_pooled, data_rbind)
    
    data_rbind[, treatment_var] <- factor(data_rbind[, treatment_var], levels = c("", levels_treatment))
    
    
  }
  
  
  data_rbind[, strat_var] <- interaction(data_rbind[, strat_vars, drop = FALSE], drop = TRUE, lex.order = lex_order, sep = " : ")
  
  table(data_rbind[, strat_var])
  
  
  characteristics_bep <- wrapper_characteristics_core(data = data_rbind, covariate_vars = covariate_vars, strat_var = strat_var, variable_names = variable_names, caption = caption, display_statistics = display_statistics)
  
  
  bheader(characteristics_bep) <- NULL
  
  
  return(characteristics_bep)
  
  
}




































