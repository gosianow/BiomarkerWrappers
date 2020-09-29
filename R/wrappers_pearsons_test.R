#' @include wrappers_logistic_regression.R
NULL




#' Pearson's Chi-squared test or Fisher's exact test or Cochran-Mantel-Haenszel Chi-Squared Test 
#' 
#' @param data Data frame.
#' @param response_var Name of categorical variable defining successes and failures where the first level corresponds to failure and the second level corresponds to success.
#' @param covariate_var Name of categorical variable defining subgroups.
#' @param strata_vars Stratification factors. If defined, then the Cochran-Mantel-Haenszel Chi-Squared Test is applied.
#' @param method Method 'pearson' or 'fisher' test.
#' @param variable_names Named vector with variable names. If not supplied, variable names are created by replacing in column names underscores with spaces.
#' @param caption Caption for the table with results.
#' @param force_empty_cols Logical. Whether to display output columns which are all empty.
#' @param print_total Logical. Whether to print total number of samples.
#' @param print_pvalues Logical. Whether to print p-values.
#'  
#' @export
wrapper_pearsons_test_core <- function(data, response_var, covariate_var, strata_vars = NULL, method = "pearson", variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_total = TRUE, print_pvalues = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(response_var) == 1)
  stopifnot(is.factor(data[, response_var]))
  stopifnot(nlevels(data[, response_var]) == 2)
  
  stopifnot(length(covariate_var) == 1)
  stopifnot(is.factor(data[, covariate_var]))
  stopifnot(nlevels(data[, covariate_var]) >= 2)
  
  
  stopifnot(method %in% c("pearson", "fisher"))
  
  if(!is.null(strata_vars)){
    strata_class <- sapply(data[, strata_vars], class)
    stopifnot(all(strata_class %in% c("factor")))
  }
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(response_var, covariate_var, strata_vars)]), ]
  
  stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  

  # --------------------------------------------------------------------------
  # Prepare data frame with results
  # --------------------------------------------------------------------------
  
  # -----------------------------------
  # Fit logistic regression to obtain Response rates
  # -----------------------------------
  
  covariate_vars <- covariate_var
  return_vars <- covariate_var
  
  wrapper_res <- wrapper_logistic_regression_core_simple(data = data, response_var = response_var, covariate_vars = covariate_vars, return_vars = return_vars, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_total = TRUE, print_pvalues = FALSE, print_adjpvalues = FALSE)
  
  
  res <- bresults(wrapper_res)
  out <- boutput(wrapper_res)
  
  ### Remove OR
  
  res <- res[, -grep("^OR", colnames(res))]
  out <- out[, -grep("^OR", colnames(out))]
  
  
  
  # --------------------------------------------------------------------------
  # Testing
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, covariate_var], data[, response_var])
  
  tbl_test <- tbl[, rev(seq_len(ncol(tbl)))]
  
  
  if(sum(margin.table(tbl, margin = 1) >= 1) >= 2 && sum(margin.table(tbl, margin = 2) >= 1) >= 2){
    
    
    if(method == "pearson"){
      
      ## Pearson's Chi-squared test
      ## Success category is defined by the first column
      prop_test_res <- prop.test(tbl_test, correct = TRUE)
      
      
      if(!is.null(prop_test_res$conf.int)){
        difference <- (prop_test_res$estimate[2] - prop_test_res$estimate[1]) * 100
        difference_CI95_upper <- -prop_test_res$conf.int[1] * 100
        difference_CI95_lower <- -prop_test_res$conf.int[2] * 100
      }else{
        difference <- NA
        difference_CI95_upper <- NA
        difference_CI95_lower <- NA
      }
      
      
      ### Calculate p-values
      if(is.null(strata_vars)){
        
        ## Using Pearson's Chi-squared test
        pvalue <- prop_test_res$p.value
        
      }else{
        
        ### Using Cochran-Mantel-Haenszel Chi-Squared Test
        mantelhaen_test_res <- NULL
        
        try(mantelhaen_test_res <- mantelhaen.test(x = data[, response_var], y = data[, covariate_var], z = interaction(data[, strata_vars]), correct = TRUE, exact = FALSE), silent = TRUE)
        
        if(is.null(mantelhaen_test_res)){
          pvalue <- NA
        }else{
          pvalue <- mantelhaen_test_res$p.value
        }
        
      }
      
      
    }else{
      
      ### Calculate OR and p-values 
      
      if(is.null(strata_vars)){
        
        ## Using Fisher's test
        fisher_test_res <- NULL
        try(fisher_test_res <- fisher.test(tbl), silent = TRUE)
        if(is.null(fisher_test_res)){
          try(fisher_test_res <- fisher.test(tbl, simulate.p.value = TRUE), silent = TRUE)
        }
        
        if(is.null(fisher_test_res)){
          pvalue <- NA
          OR <- NA
          OR_CI95_lower <- NA
          OR_CI95_upper <- NA
        }else{
          pvalue <- fisher_test_res$p.value
          OR <- fisher_test_res$estimate
          OR_CI95_lower <- fisher_test_res$conf.int[1]
          OR_CI95_upper <- fisher_test_res$conf.int[2]
          if(is.null(OR)){
            OR <- NA
            OR_CI95_lower <- NA
            OR_CI95_upper <- NA
          }
        }
        
        
      }else{
        
        ### Using Cochran-Mantel-Haenszel Chi-Squared Test
        mantelhaen_test_res <- NULL
        
        try(mantelhaen_test_res <- mantelhaen.test(x = data[, response_var], y = data[, covariate_var], z = interaction(data[, strata_vars]), correct = TRUE, exact = FALSE), silent = TRUE)
        
        if(is.null(mantelhaen_test_res)){
          OR <- NA
          OR_CI95_lower <- NA
          OR_CI95_upper <- NA
          pvalue <- NA
        }else{
          OR <- mantelhaen_test_res$estimate
          OR_CI95_lower <- mantelhaen_test_res$conf.int[1]
          OR_CI95_upper <- mantelhaen_test_res$conf.int[2]
          pvalue <- mantelhaen_test_res$p.value
          if(is.null(OR)){
            OR <- NA
            OR_CI95_lower <- NA
            OR_CI95_upper <- NA
          }
        }
        
      }
      
      
    }
    
    
    
    
  }else{
    
    difference <- NA
    difference_CI95_lower <- NA
    difference_CI95_upper <- NA
    
    OR <- NA
    OR_CI95_lower <- NA
    OR_CI95_upper <- NA
    
    pvalue <- NA
    
  }
  
  
  ### Add statistics
  
  if(method == "pearson"){
    res$difference = c(difference, rep(NA, nrow(tbl) - 1))
    res$difference_CI95_lower = c(difference_CI95_lower, rep(NA, nrow(tbl) - 1))
    res$difference_CI95_upper = c(difference_CI95_upper, rep(NA, nrow(tbl) - 1))
  }else{
    res$OR = c(OR, rep(NA, nrow(tbl) - 1))
    res$OR_CI95_lower = c(OR_CI95_lower, rep(NA, nrow(tbl) - 1))
    res$OR_CI95_upper = c(OR_CI95_upper, rep(NA, nrow(tbl) - 1))
  }
  
  
  ### Add pvalue 
  res$pvalue <- c(pvalue, rep(NA, nrow(res) - 1))
  
  ### Remove adj_pvalue
  res$adj_pvalue <- NULL
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  if(method == "pearson"){
    out$Difference = format_difference(res$difference, digits = 2, non_empty = 1)
    out$`Difference 95% CI` = format_CIs(res$difference_CI95_lower, res$difference_CI95_upper, digits = 2, non_empty = 1)
    ### If all Difference are empty, do not display that column.
    if(all(out$Difference %in% c("", "NA")) && !force_empty_cols){
      out$Difference <- NULL
      out$`Difference 95% CI` <- NULL
    }
    
  }else{
    out$OR = format_or(res$OR, digits = 2, non_empty = 1)
    out$`OR 95% CI` = format_CIs(res$OR_CI95_lower, res$OR_CI95_upper, digits = 2, non_empty = 1)
    
    ### If all OR are empty, do not display that column.
    if(all(out$OR %in% c("", "NA")) && !force_empty_cols){
      out$OR <- NULL
      out$`OR 95% CI` <- NULL
    }
    
  }
  
  
  out$`P-value` <- format_pvalues(res$pvalue, non_empty = 1)
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  
  if(!print_total){
    col_total <- grep("^Total", colnames(out), value = TRUE)
    for(i in seq_along(col_total)){
      out[, col_total[i]] <- NULL
    }
  }
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    if(method == "pearson" && is.null(strata_vars)){
      caption <- paste0("Unstratified analysis with Pearson's Chi-squared test.")
    }else if(method == "fisher" && is.null(strata_vars)){
      caption <- paste0("Unstratified analysis with Fisher's exact test.")
    }else{
      caption <- paste0("Stratified analysis with Cochran-Mantel-Haenszel Chi-squared test.", " Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ".")
    }
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  return(bout)
  
  
}




#' @rdname wrapper_pearsons_test_core
#' @inheritParams wrapper_pearsons_test_core
#' @param strat1_var Name of the first stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @param print_adjpvalues Logical. Whether to print adjusted p-values.
#' @export
wrapper_pearsons_test_core_strat <- function(data, response_var, covariate_var, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, method = "pearson", variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
  
  ### Strata cannot include response_var, covariate_var
  stopifnot(length(intersect(c(strata_vars, strat1_var, strat2_var), c(response_var, covariate_var))) == 0)
  
  if(print_adjpvalues){
    print_pvalues <- TRUE
  }
  
  ### Keep non-missing data
  data <- data[complete.cases(data[, c(response_var, covariate_var, strata_vars, strat1_var, strat2_var)]), , drop = FALSE]
  
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
      
      
      wrapper_res <- wrapper_pearsons_test_core(data = data_strata1, response_var = response_var, covariate_var = covariate_var, strata_vars = strata_vars, method = method, variable_names = variable_names, caption = caption, force_empty_cols = force_empty_cols, print_total = print_total, print_pvalues = print_pvalues)
      
      
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
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if(print_adjpvalues){
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
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
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]))
  
  return(wrapper_res)
  
  
}







#' Testing biomarker effect on response with Pearson's Chi-squared test or Cochran-Mantel-Haenszel Chi-Squared Test
#' 
#' @inheritParams wrapper_pearsons_test_core_strat
#' @param biomarker_vars Vector of biomarker names.
#' @export
wrapper_pearsons_test_biomarker <- function(data, response_var, biomarker_vars, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, method = "pearson", variable_names = NULL, caption = NULL, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(biomarker_vars), c(strata_vars, strat1_var, strat2_var))) == 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 1
    
    covariate_var <- biomarker_vars[i]
    
    
    wrapper_res <- wrapper_pearsons_test_core_strat(data = data, response_var = response_var, covariate_var = covariate_var, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, method = method, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_total = print_total, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
    return(wrapper_res)
    
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
  }
  
  
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference == "")){
    out$Difference <- NULL
    out$`Difference 95% CI` <- NULL
  }
  
  ### Rename 'Covariate' column name to 'Biomarker'
  
  colnames(res)[colnames(res) == "covariate"] <- "biomarker"
  colnames(out)[colnames(out) == "Covariate"] <- "Biomarker"
  
  
  ### Generate caption
  
  
  if(is.null(caption)){
    
    caption <- paste0("Biomarker effect on ", variable_names[response_var], ". ")
    
    if(method == "pearson" && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified Pearson's Chi-squared test.")
      
    }else if(method == "fisher" && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified Fisher's exact test.")
      
    }else{
      
      caption <- paste0(caption, "Stratified Cochran-Mantel-Haenszel Chi-squared test. Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}






#' Testing treatment effect on response with Pearson's Chi-squared test or Cochran-Mantel-Haenszel Chi-Squared Test
#' 
#' @inheritParams wrapper_pearsons_test_core_strat
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector with names of categorical biomarkers. When NULL, overall treatment effect is estimated. 
#' @export
wrapper_pearsons_test_treatment <- function(data, response_var, treatment_var, strata_vars = NULL, biomarker_vars = NULL, strat2_var = NULL, method = "pearson", variable_names = NULL, caption = NULL, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(treatment_var), c(strata_vars, biomarker_vars, strat2_var))) == 0)
  
  
  ### Allow biomarker_vars = NULL
  if(is.null(biomarker_vars)){
    ### Add dummy variable to data
    data[, "biomarker_dummy"] <- factor("biomarker_dummy")
    biomarker_vars <- "biomarker_dummy"
  }
  
  
  vars_class <- sapply(data[, c(treatment_var, biomarker_vars)], class)
  stopifnot(all(vars_class == "factor"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)

  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 1
    
    covariate_var <- treatment_var
    strat1_var <- biomarker_vars[i]
    
    wrapper_res <- wrapper_pearsons_test_core_strat(data = data, response_var = response_var, covariate_var = covariate_var, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, method = method, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_total = print_total, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
    res <- bresults(wrapper_res)
    out <- boutput(wrapper_res)
    
    
    ### Rename strata column name to 'Biomarker Subgroup'
    
    colnames(res)[colnames(res) == strat1_var] <- "biomarker_subgroup"
    colnames(out)[colnames(out) == variable_names[strat1_var]] <- "Biomarker Subgroup"
    
    ### Treatment is the same for all the biomarkers and biomarker info is missing. Thus, we remove covariate and add biomarker.
    
    res[, "biomarker"] <- strat1_var
    res[res[, "biomarker_subgroup"] == "", "biomarker"] <- ""
    out[, "Biomarker"] <- variable_names[strat1_var]
    out[out[, "Biomarker Subgroup"] == "", "Biomarker"] <- ""
    
    res[, "covariate"] <- NULL
    out[, "Covariate"] <- NULL
    
    
    res <- dplyr::select(res, c(strat2_var, "biomarker", "biomarker_subgroup"), everything())
    out <- dplyr::select(out, c(as.character(variable_names[strat2_var]), "Biomarker", "Biomarker Subgroup"), everything())
    
    
    bresults(wrapper_res) <- res
    boutput(wrapper_res) <- out
    
    
    return(wrapper_res)
    
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
  }
  
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference == "")){
    out$Difference <- NULL
    out$`Difference 95% CI` <- NULL
  }
  
  if(length(biomarker_vars) == 1){
    if(biomarker_vars == "biomarker_dummy"){
      res$biomarker <- NULL
      out$`Biomarker` <- NULL
      res$biomarker_subgroup <- NULL
      out$`Biomarker Subgroup` <- NULL
    }
  }
  
  
  ### Change the name to Treatment Subgroup
  
  colnames(res)[colnames(res) == "subgroup"] <- "treatment_subgroup"
  colnames(out)[colnames(out) == "Subgroup"] <- "Treatment Subgroup"
  
  ### Generate caption
  
  if(is.null(caption)){
    
    caption <- paste0("Treatment effect on ", variable_names[response_var], ". ")
    
    if(method == "pearson" && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified Pearson's Chi-squared test.")
      
    }else if(method == "fisher" && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified Fisher's exact test.")
      
    }else{
      
      caption <- paste0(caption, "Stratified Cochran-Mantel-Haenszel Chi-squared test. Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}















