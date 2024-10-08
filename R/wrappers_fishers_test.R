





#' Fisher's test
#' 
#' @param data Data frame.
#' @export
wrapper_fishers_test_core <- function(data, col_var, row_var, variable_names = NULL, caption = NULL, margin = 1, force_empty_cols = FALSE, print_pvalues = TRUE, print_OR = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(col_var) == 1)
  stopifnot(is.factor(data[, col_var]))
  stopifnot(nlevels(data[, col_var]) >= 2)
  
  stopifnot(length(row_var) == 1)
  stopifnot(is.factor(data[, row_var]))
  stopifnot(nlevels(data[, row_var]) >= 2)
  
  
  ### Keep non-missing data
  all_vars <- c(col_var, row_var)
  data <- data[stats::complete.cases(data[, all_vars]), all_vars, drop = FALSE]
  stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions and do testing
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, row_var], data[, col_var])

  prop <- prop.table(tbl, margin = margin) * 100
  
  
  if(sum(margin.table(tbl, margin = 1) >= 1) >= 2 && sum(margin.table(tbl, margin = 2) >= 1) >= 2){
    ## Fisher's exact test: get odds rations and CI for OR
    
    test_res <- NULL
    try(test_res <- fisher.test(tbl), silent = TRUE)
    if(is.null(test_res)){
      try(test_res <- fisher.test(tbl, simulate.p.value = TRUE), silent = TRUE)
    }
    
    if(is.null(test_res)){
      pvalue <- NA
      OR <- NA
    }else{
      pvalue <- test_res$p.value
      OR <- test_res$estimate
      if(is.null(OR)){
        OR <- NA
      }
    }
    
  }else{
    pvalue <- NA
    OR <- NA
  }
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  countdf <- as.data.frame.matrix(tbl)
  colnames(countdf) <- paste0("counts_", colnames(countdf))
  
  propdf <- as.data.frame.matrix(prop)
  colnames(propdf) <- paste0("proportions_", colnames(propdf))
  
  
  res <- data.frame(covariate = row_var,
    subgroup = rownames(tbl),
    countdf, propdf,
    OR = c(OR, rep(NA, nrow(tbl) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(tbl) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = res$subgroup, 
    format_counts_and_props_df(counts = countdf, props = propdf, digits = 1),
    OR = format_or(res$OR, digits = 2, non_empty = 1),
    `P-value` = format_pvalues(res$pvalue, non_empty = 1), 
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  if(!print_OR){
    out$OR <- NULL
  }
  
  ### If all OR are empty, do not display that column.
  if(all(out$OR %in% c("", "NA")) && !force_empty_cols){
    out$OR <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  header <- format_header(all_colnames = colnames(out), header_colnames = levels(data[, col_var]), header_name = variable_names[col_var])
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Fisher's exact test.")
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}





#' @rdname wrapper_fishers_test_core
#' @inheritParams wrapper_fishers_test_core
#' @param strat1_var Name of the first stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_fishers_test_core_strat <- function(data, col_var, row_var, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, margin = 1, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE, print_OR = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check on strat vars
  # --------------------------------------------------------------------------
  
  if(!is.null(strat1_var)){
    stopifnot(length(strat1_var) == 1)
    stopifnot(is.factor(data[, strat1_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat1_dummy"] <- factor("strat1_dummy")
    strat1_var <- "strat1_dummy"
  }
  
  if(!is.null(strat2_var)){
    stopifnot(length(strat2_var) == 1)
    stopifnot(is.factor(data[, strat2_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat2_dummy"] <- factor("strat2_dummy")
    strat2_var <- "strat2_dummy"
  }
  
  ### Strata cannot include col_var, row_var
  stopifnot(length(intersect(c(strat1_var, strat2_var), c(col_var, row_var))) == 0)
  
  if(!print_pvalues){
    print_adjpvalues <- FALSE
  }
  
  ### Keep non-missing data
  all_vars <- c(col_var, row_var, strat1_var, strat2_var)
  data <- data[stats::complete.cases(data[, all_vars]), all_vars, drop = FALSE]
  stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculations within strata
  # --------------------------------------------------------------------------
  
  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  wrapper_res <- lapply(1:length(strata2_levels), function(j){
    # j = 1
    
    data_strata2 <- data[data[, strat2_var] == strata2_levels[j] & !is.na(data[, strat2_var]), ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    wrapper_res <- lapply(1:length(strata1_levels), function(i){
      # i = 1
      
      # print(paste(j, i))
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      wrapper_res <- wrapper_fishers_test_core(data = data_strata1, col_var = col_var, row_var = row_var, variable_names = variable_names, caption = caption, margin = margin, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues, print_OR = print_OR)
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      
      rownames(res) <- NULL
      rownames(out) <- NULL
      
      wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res))
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if(print_adjpvalues){
    out$`Adj. P-value` <- format_pvalues(stats::p.adjust(res$pvalue, method = "BH"))
  }
  
  
  ### Set repeating Strata names to empty
  out[out$Covariate == "", variable_names[strat1_var]] <- ""
  out[out$Covariate == "", variable_names[strat2_var]] <- ""
  
  ### Remove dummy columns
  if(strat2_var == "strat2_dummy"){
    res$strat2_dummy <- NULL
    out$`strat2 dummy` <- NULL
  }
  if(strat1_var == "strat1_dummy"){
    res$strat1_dummy <- NULL
    out$`strat1 dummy` <- NULL
  }
  
  ## Update header
  hdr <- format_header(all_colnames = colnames(out), header_colnames = levels(data[, col_var]), header_name = variable_names[col_var])
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
  
}





#' Fisher's test
#' 
#' Run Fisher's test for multiple covariates.
#' 
#' @inheritParams wrapper_fishers_test_core_strat
#' @param row_vars Vector with names of categorical variables.
#' @export
wrapper_fishers_test <- function(data, col_var, row_vars, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, margin = 1, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE, print_OR = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  stopifnot(length(col_var) == 1)
  stopifnot(length(row_vars) >= 1)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  
  wrapper_res <- lapply(1:length(row_vars), function(i){
    # i = 1
    
    row_var <- row_vars[i]
    
    wrapper_res <- wrapper_fishers_test_core_strat(data, col_var = col_var, row_var = row_var, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, margin = margin, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues, print_OR = print_OR)
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  hdr <- bheader(wrapper_res[[1]])
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
  }
  
  
  if("OR" %in% colnames(out)){
    ### If all OR are empty, do not display that column.
    if(all(out$OR %in% c("", "NA")) && !force_empty_cols){
      out$OR <- NULL
      ### Update header
      if(hdr[3] == 1){
        hdr <- hdr[-3]
      }else{
        hdr[3] <- hdr[3] - 1
      }
    }
  }
  
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
}















