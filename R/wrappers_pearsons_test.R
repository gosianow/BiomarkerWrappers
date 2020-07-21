






#' Pearson's Chi-squared test or Cochran-Mantel-Haenszel Chi-Squared Test 
#' 
#' @param data Data frame.
#' @param response_var Name of categorical variable defining successes and failures where the first level corresponds to failure and the second level corresponds to success.
#' @param covariate_var Name of categorical variable defining subgroups.
#' @param strata_vars Stratification factors. If defined, then the Cochran-Mantel-Haenszel Chi-Squared Test is applied. 
#' @export
wrapper_core_pearsons_test <- function(data, response_var, covariate_var, strata_vars = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_pvalues = TRUE){
  
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
  
  
  method <- "pearson"
  
  if(!is.null(strata_vars)){
    strata_class <- sapply(data[, strata_vars], class)
    stopifnot(all(strata_class %in% c("factor")))
    method <- "mantelhaen"
  }
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(response_var, covariate_var, strata_vars)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions and do testing
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, covariate_var], data[, response_var])
  
  margin <- 1
  
  prop <- prop.table(tbl, margin = margin) * 100
  
  tbl_test <- tbl[, rev(seq_len(ncol(tbl)))]
  
  success_level <- colnames(tbl_test)[1]
  
  
  if(sum(margin.table(tbl, margin = 1) >= 1) >= 2 && sum(margin.table(tbl, margin = 2) >= 1) >= 2){
    ## Pearson's Chi-squared test: get difference of proportions of succeses and CI for difference
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
    if(method == "pearson"){
      
      ## Using Pearson's Chi-squared test
      pvalue <- prop_test_res$p.value
      
    }else{
      
      ### Using Cochran-Mantel-Haenszel Chi-Squared Test
      mantelhaen_test_res <- mantelhaen.test(x = data[, response_var], y = data[, covariate_var], z = interaction(data[, strata_vars]), correct = TRUE, exact = FALSE)
      
      pvalue <- mantelhaen_test_res$p.value
      
    }
    
    
  }else{
    
    difference <- NA
    difference_CI95_lower <- NA
    difference_CI95_upper <- NA
    pvalue <- NA
    
  }
  
  
  ### Calculate CIs for proportions of success 
  
  response_CI <- lapply(1:nrow(tbl_test), function(k){
    # k = 1
    
    binom_test_res <- binom.test(as.numeric(tbl_test[k, ]))
    
    out <- data.frame(response_CI95_lower = binom_test_res$conf.int[1] * 100, response_CI95_upper = binom_test_res$conf.int[2] * 100)
    
  })
  
  response_CI <- rbind.fill(response_CI)
  colnames(response_CI) <- paste0("response_", c("CI95_lower", "CI95_upper"))
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  res <- data.frame(covariate = covariate_var,
    subgroup = rownames(tbl),
    n_total = sum(tbl),
    n = as.numeric(margin.table(tbl, margin = 1)),
    nresponse = as.numeric(tbl[, success_level]),
    propresponse = as.numeric(prop[, success_level]),
    response_CI,
    difference = c(difference, rep(NA, nrow(tbl) - 1)),
    difference_CI95_lower = c(difference_CI95_lower, rep(NA, nrow(tbl) - 1)),
    difference_CI95_upper = c(difference_CI95_upper, rep(NA, nrow(tbl) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(tbl) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = res$subgroup, 
    `Total N` = as.character(res$n_total),
    `N` = as.character(res$n),
    Response = format_counts_and_props_core(counts = res$nresponse, props = res$propresponse),
    `Response 95% CI` = format_CIs(res$response_CI95_lower, res$response_CI95_upper),
    Difference = format_difference(res$difference, digits = 2),
    `Difference 95% CI` = format_CIs(res$difference_CI95_lower, res$difference_CI95_upper, digits = 2),
    `P-value` = format_pvalues(res$pvalue), 
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference == "") && !force_empty_cols){
    out$Difference <- NULL
    out$`Difference 95% CI` <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    if(method == "pearson"){
      caption <- paste0("Unstratified analysis with Pearson's Chi-squared test.")
    }else{
      caption <- paste0("Stratified analysis with Cochran-Mantel-Haenszel Chi-squared test.", " Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ".")
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  return(bout)
  
  
}




#' @rdname wrapper_core_pearsons_test
#' @inheritParams wrapper_core_pearsons_test
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_core_pearsons_test_strat <- function(data, response_var, covariate_var, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
      
      
      wrapper_res <- wrapper_core_pearsons_test(data = data_strata1, response_var = response_var, covariate_var = covariate_var, strata_vars = strata_vars, variable_names = variable_names, caption = caption, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues)
      
      
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
#' @inheritParams wrapper_core_pearsons_test_strat
#' @param biomarker_vars Vector of biomaker names.
#' @export
wrapper_pearsons_test_biomarker <- function(data, response_var, biomarker_vars, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
    
    
    wrapper_res <- wrapper_core_pearsons_test_strat(data = data, response_var = response_var, covariate_var = covariate_var, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
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
    
    if(is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified Pearson's Chi-squared test. Subgroups defined by the biomarker.")
      
    }else{
      
      caption <- paste0(caption, "Stratified Cochran-Mantel-Haenszel Chi-squared test. Subgroups defined by the biomarker. Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}






#' Testinging treatment effect on response with Pearson's Chi-squared test or Cochran-Mantel-Haenszel Chi-Squared Test
#' 
#' @inheritParams wrapper_core_pearsons_test_strat
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector with names of categorical biomarkers. When NULL, overall treatment effect is estimated. 
#' @export
wrapper_pearsons_test_treatment <- function(data, response_var, treatment_var, strata_vars = NULL, biomarker_vars = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
    
    wrapper_res <- wrapper_core_pearsons_test_strat(data = data, response_var = response_var, covariate_var = covariate_var, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
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
    out <- dplyr::select(out, c(variable_names[strat2_var], "Biomarker", "Biomarker Subgroup"), everything())
    
    
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
    
    if(is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified Pearson's Chi-squared test. Subgroups defined by the treatment.")
      
    }else{
      
      caption <- paste0(caption, "Stratified Cochran-Mantel-Haenszel Chi-squared test. Subgroups defined by the treatment. Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}















