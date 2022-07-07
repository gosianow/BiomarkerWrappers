







#' Logistic regression with simple additive model
#' 
#' @param data Data frame.
#' @param response_var Name of the response variable. This variable must be a factor where 'success' is interpreted as the factor not having the first level.
#' @param covariate_vars Vector with names of covariate that should be included in the formula.
#' @param return_vars Vector with names of covariate that for which the statistics should be returned. If NULL, statistics for all covariates are returned.
#' @param variable_names Named vector with variable names. If not supplied, variable names are created by replacing in column names underscores with spaces.
#' @param caption Caption for the table with results.
#' @param force_empty_cols Logical. Whether to display output columns which are all empty.
#' @param print_total Logical. Whether to print total number of samples.
#' @param print_pvalues Logical. Whether to print p-values.
#' @param print_adjpvalues Logical. Whether to print adjusted p-values.
#' @details 
#' If for a factor covariate that should be returned the reference level has zero counts, results are set to NA because this levels is not used as a reference which means that it is not possible to estimate odds ratios that we want.
#' @export
wrapper_logistic_regression_core_simple <- function(data, response_var, covariate_vars, return_vars = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_total = TRUE, print_non_response = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE, print_OR = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(isValidAndUnreservedName(c(response_var, covariate_vars)))
  
  stopifnot(is.data.frame(data))
  
  ## Response variable must be factor with two levels
  stopifnot(length(response_var) == 1)
  stopifnot(is.factor(data[, response_var]))
  stopifnot(nlevels(data[, response_var]) == 2)
  
  
  covariate_class <- sapply(data[, covariate_vars], class)
  stopifnot(all(covariate_class %in% c("factor", "numeric", "integer")))
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(response_var, covariate_vars)]), ]
  
  stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  response_levels <- levels(data[, response_var])
  
  
  # --------------------------------------------------------------------------
  # Generate data frame with coefficient names and levels and information about reference groups
  # --------------------------------------------------------------------------
  
  coef_info <- lapply(1:length(covariate_vars), function(i){
    # i = 1
    
    if(covariate_class[i] %in% c("numeric", "integer")){
      
      
      out_dummy <- as.data.frame.matrix(matrix(NA, nrow = 1, ncol = 2 * length(response_levels) + 2))
      colnames(out_dummy) <- c(paste0("n_", response_levels), paste0("prop_", response_levels), paste0(response_levels[2], c("_CI95_lower", "_CI95_upper")))
      
      
      res <- data.frame(covariate = covariate_vars[i], covariate_class = covariate_class[i], 
        subgroup = "", reference = "", reference_indx = 0, n = NA, 
        out_dummy,
        OR_non_empty = TRUE,
        stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
      
      
      return(res)
      
      
    }else{
      
      ## Check if the first level has non zero counts so it can be used as a reference group in the regression
      ## When the first level has zero counts, then the next level with non-zero counts is used as a reference 
      tbl <- table(data[, covariate_vars[i]])
      
      reference_indx <- which(tbl > 0)[1]
      
      ## Calculate nresponse and propresponse
      tbl_response <- table(data[, covariate_vars[i]], data[, response_var])
      prop_response <- prop.table(tbl_response, margin = 1) * 100
      ## Replace NaN with NA
      prop_response[is.na(prop_response)] <- NA
      
      
      ndf <- as.data.frame.matrix(tbl_response)
      colnames(ndf) <- paste0("n_", colnames(tbl_response))
      
      propdf <- as.data.frame.matrix(prop_response)
      colnames(propdf) <- paste0("prop_", colnames(prop_response))
      
      
      
      res <- data.frame(covariate = covariate_vars[i], covariate_class = covariate_class[i], 
        subgroup = levels(data[, covariate_vars[i]]), reference = names(reference_indx), reference_indx = as.numeric(reference_indx), n = as.numeric(tbl), 
        ndf, propdf,
        stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
      
      res$OR_non_empty <- c(FALSE, rep(TRUE, nrow(res) - 1))
      
      
      # --------------------------------------------------------------------------
      # Calculate CIs for proportions of response/success 
      # --------------------------------------------------------------------------

      
      tbl_test <- tbl_response[, rev(seq_len(ncol(tbl_response))), drop = FALSE]
      
      
      response_CI <- lapply(1:nrow(tbl_test), function(k){
        # k = 1
        
        binom_test_res <- binom.test(as.numeric(tbl_test[k, ]))
        
        res <- data.frame(response_CI95_lower = binom_test_res$conf.int[1] * 100, response_CI95_upper = binom_test_res$conf.int[2] * 100)
        
      })
      
      response_CI <- rbind.fill(response_CI)
      colnames(response_CI) <- paste0(response_levels[2], c("_CI95_lower", "_CI95_upper"))
      
      
      res <- cbind(res, response_CI)
      
      return(res) 
      
    }
    
    
  })
  
  coef_info <- plyr::rbind.fill(coef_info)
  
  coef_info$n_total <- nrow(data)
  
  coef_info$coefficient <- paste0(coef_info$covariate, coef_info$subgroup)
  
  rownames(coef_info) <- coef_info$coefficient
  
  
  # --------------------------------------------------------------------------
  # Logistic regression
  # --------------------------------------------------------------------------
  
  
  ## Create the formula
  formula_covariates <- paste0(covariate_vars, collapse = " + ")
  f <- stats::as.formula(paste0(response_var, " ~ ", formula_covariates))
  
  
  ## Fit the logistic model
  regression_fit <- NULL
  
  try(regression_fit <- glm(formula = f, family = binomial(link = "logit"), data = data), silent = TRUE)
  
  if(is.null(regression_fit)){
    regression_summ <- NULL
  }else{
    regression_summ <- summary(regression_fit)
  }
  
  
  # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
  # h(mm)
  
  
  
  # --------------------------------------------------------------------------
  # Calculate confidence intervals
  # --------------------------------------------------------------------------
  
  
  ## Note: When MASS is loaded, the 95% CIs are calculated (using confint) based on profile likelihood
  ## To compute 95% Wald CIs (based on asymptotic normality), one needs to use confint.default()
  
  ## There can be an error from confint when all samples have the same response
  # confint_res <- NULL
  # try(confint_res <- exp(stats::confint(regression_fit)), silent = TRUE)
  # if(is.null(confint_res)){
  #   confint_res <- exp(stats::confint.default(regression_fit))
  # }
  
  ## I use the Wald CIs because they are in concordance with p-values i.e. they contain 1 when p-value is not significant
  
  
  if(!is.null(regression_fit)){
    confint_res <- exp(stats::confint.default(regression_fit))
  }
  
  
  # --------------------------------------------------------------------------
  # Parse the regression summary
  # --------------------------------------------------------------------------
  
  
  if(is.null(regression_summ)){
    
    conf_int <- data.frame(coefficient = coef_info$coefficient, OR_CI95_lower = NA, OR_CI95_upper = NA, stringsAsFactors = FALSE)
    coefficients <- data.frame(coefficient = coef_info$coefficient, OR = NA, pvalue = NA, stringsAsFactors = FALSE)
    
  }else{
    
    conf_int <- data.frame(coefficient = rownames(confint_res), confint_res[, c("2.5 %", "97.5 %"), drop = FALSE], stringsAsFactors = FALSE)
    colnames(conf_int) <- c("coefficient", "OR_CI95_lower", "OR_CI95_upper")
    
    
    coefficients <- data.frame(coefficient = rownames(regression_summ$coefficients), regression_summ$coefficients[, c("Estimate", "Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE)
    
    colnames(coefficients) <- c("coefficient", "OR", "pvalue")
    
    coefficients$OR <- exp(coefficients$OR)
    
  }
  

  # --------------------------------------------------------------------------
  # Append results from regression
  # --------------------------------------------------------------------------
  
  
  coef_info <- coef_info %>% 
    dplyr::left_join(conf_int, by = "coefficient") %>% 
    dplyr::left_join(coefficients, by = "coefficient")
  
  
  ## Calculate adjusted p-values using the Benjamini & Hochberg method
  coef_info$adj_pvalue <- stats::p.adjust(coef_info$pvalue, method = "BH")
  
  
  # --------------------------------------------------------------------------
  # Return results for the covariates of interest defined by `return_vars`
  # --------------------------------------------------------------------------
  
  if(is.null(return_vars)){
    
    ## Return results for all covariates
    res <- coef_info
    
  }else{
    
    ## Return results for covariates defined by `return_vars`
    stopifnot(any(return_vars %in% coef_info$covariate))
    
    res <- coef_info[coef_info$covariate %in% return_vars, , drop = FALSE]
    
    ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA because eventually this level is not used as a reference in the fitted model.
    if(all(res$reference_indx > 1)){
      
      res[, colnames(res) %in% c("OR", "OR_CI95_lower", "OR_CI95_upper", "pvalue", "adj_pvalue")] <- NA
      
    }
    
  }
  
  
  # --------------------------------------------------------------------------
  # Prepare the output data frame that will be displayed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = as.character(res$subgroup),
    `Total N` = as.character(res$n_total),
    `N` = format_difference(res$n, digits = 0),
    format_counts_and_props_df(counts = res[, paste0("n_", response_levels)], props = res[, paste0("prop_", response_levels)], digits = 1, prefix_counts = "n_"),
    
    as.data.frame(matrix(format_CIs(res[, paste0(response_levels[2], "_CI95_lower")], res[, paste0(response_levels[2], "_CI95_upper")], non_empty = res$covariate_class == "factor"), ncol = 1, dimnames = list(NULL, paste0(response_levels[2], " 95% CI")))),

    `OR` = format_or(res$OR, non_empty = res$OR_non_empty),
    `OR 95% CI` = format_CIs(res$OR_CI95_lower, res$OR_CI95_upper, non_empty = res$OR_non_empty),
    `P-value` = format_pvalues(res$pvalue, non_empty = res$OR_non_empty),
    `Adj. P-value` = format_pvalues(res$adj_pvalue, non_empty = res$OR_non_empty),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  ### When all covariates are numerical, some outputs are empty, and we remove them
  if(!force_empty_cols){
    for(i in seq_len(ncol(out))){
      if(all(out[, i] %in% "")){
        out[, i] <- NULL
      }
    }
  }
  
  if(!print_total){
    col_total <- grep("^Total", colnames(out), value = TRUE)
    for(i in seq_along(col_total)){
      out[, col_total[i]] <- NULL
    }
  }
  
  if(!print_non_response){
    out[, response_levels[1]] <- NULL
  }
  
  if(!print_pvalues){
    out$`P-value` <- NULL
    out$`Adj. P-value` <- NULL
  }
  
  if(!print_adjpvalues){
    out$`Adj. P-value` <- NULL
  }
  
  if(!print_OR){
    out$OR <- NULL
    out$`OR 95% CI` <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Covariate effect on ", variable_names[response_var], 
      ". Logistic regression model includes ", paste0(variable_names[covariate_vars], collapse = ", "), ".")
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}







#' @rdname wrapper_logistic_regression_core_simple
#' @inheritParams wrapper_logistic_regression_core_simple
#' @param strat1_var Name of the first stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_logistic_regression_core_simple_strat <- function(data, response_var, covariate_vars, return_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_total = TRUE, print_non_response = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE, print_OR = TRUE){
  
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
  
  ### Covariates cannot include strata
  stopifnot(length(intersect(c(strat1_var, strat2_var), covariate_vars)) == 0)
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(response_var, covariate_vars, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
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
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      wrapper_res <- wrapper_logistic_regression_core_simple(data = data_strata1, response_var = response_var, covariate_vars = covariate_vars, return_vars = return_vars, variable_names = variable_names, caption = caption, force_empty_cols = force_empty_cols, print_total = print_total, print_non_response = print_non_response, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues, print_OR = print_OR)
      
      
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      
      
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      bresults(wrapper_res) <- res
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      boutput(wrapper_res) <- out
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    ## This should work now too
    # wrapper_res <- do.call(rbind, wrapper_res)
    
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
  
  if("Adj. P-value" %in% colnames(out)){
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









#' Logistic regression estimating biomarker effect 
#' 
#' @inheritParams wrapper_logistic_regression_core_simple_strat
#' @param biomarker_vars Vector of biomarker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
#' @export
wrapper_logistic_regression_biomarker <- function(data, response_var, biomarker_vars, treatment_var = NULL, adjustment_vars = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_total = TRUE, print_non_response = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE, print_OR = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include biomarker variables
  stopifnot(length(intersect(biomarker_vars, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(biomarker_vars, adjustment_vars), c(treatment_var, strat2_var))) == 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 2
    
    covariate_vars <- c(biomarker_vars[i], adjustment_vars)
    return_vars <- biomarker_vars[i]
    
    
    wrapper_res <- wrapper_logistic_regression_core_simple_strat(data = data, response_var = response_var, covariate_vars = covariate_vars, return_vars = return_vars, strat1_var = treatment_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_total = print_total, print_non_response = print_non_response, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues, print_OR = print_OR)
    
    
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
  
  
  ### When all covariates are numerical, some outputs are empty, and we remove them
  for(i in seq_len(ncol(out))){
    if(all(out[, i] %in% "")){
      out[, i] <- NULL
    }
  }
  
  
  ### Rename 'Covariate' column name to 'Biomarker'
  
  colnames(res)[colnames(res) == "covariate"] <- "biomarker"
  colnames(out)[colnames(out) == "Covariate"] <- "Biomarker"
  
  
  ### Generate caption
  
  
  if(is.null(caption)){
    
    caption <- paste0("Biomarker effect on ", variable_names[response_var], ". ")
    
    
    if(is.null(adjustment_vars)){
      
      caption <- paste0(caption, "Unadjusted analysis. Logistic regression model includes only the biomarker.")
      
    }else{
      
      caption <- paste0(caption, "Adjusted analysis. Logistic regression model includes the biomarker and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}






#' Logistic regression estimating treatment effect within biomarker subgroups
#' 
#' @inheritParams wrapper_logistic_regression_core_simple_strat
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector of biomarker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
#' @export
wrapper_logistic_regression_treatment <- function(data, response_var, treatment_var, biomarker_vars = NULL, adjustment_vars = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_total = TRUE, print_non_response = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE, print_OR = TRUE){
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include treatment
  stopifnot(length(intersect(treatment_var, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(treatment_var, adjustment_vars), c(biomarker_vars, strat2_var))) == 0)
  
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
    
    covariate_vars <- c(treatment_var, adjustment_vars)
    return_vars <- treatment_var
    strat1_var <- biomarker_vars[i]
    
    
    wrapper_res <- wrapper_logistic_regression_core_simple_strat(data = data, response_var = response_var, covariate_vars = covariate_vars, return_vars = return_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_total = print_total, print_non_response = print_non_response, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues, print_OR = print_OR)
    
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
    
    if(is.null(adjustment_vars)){
      
      caption <- paste0(caption, "Unadjusted analysis. Logistic regression model includes only the treatment.")
      
    }else{
      
      caption <- paste0(caption, "Adjusted analysis. Logistic regression model includes the treatment and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}



###############################################################################
# Model with interaction
###############################################################################


# response_var = "Response_by_pCR"; interaction1_var = "Treatment_Arm"; interaction2_var = "AKT.S473.D9E"; covariate_vars = c("Ipat_Dx"); variable_names = NULL; caption = NULL; print_pvalues = TRUE; print_adjpvalues = TRUE



#' Logistic regression with additive model with interaction
#' 
#' @inheritParams wrapper_logistic_regression_core_simple
#' @param interaction1_var Name of the first interaction variable. 
#' @param interaction2_var Name of the second interaction variable.
#' @export
wrapper_logistic_regression_core_interaction <- function(data, response_var, interaction1_var, interaction2_var, covariate_vars = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  ## Response variable must be factor with two levels
  stopifnot(length(response_var) == 1)
  stopifnot(is.factor(data[, response_var]))
  stopifnot(nlevels(data[, response_var]) == 2)
  
  covariate_class <- sapply(data[, c(interaction1_var, interaction2_var, covariate_vars)], class)
  
  stopifnot(all(covariate_class %in% c("factor", "numeric", "integer")))
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(response_var, interaction1_var, interaction2_var, covariate_vars)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(print_adjpvalues){
    print_pvalues <- TRUE
  }

  
  # --------------------------------------------------------------------------
  # Logistic regression
  # --------------------------------------------------------------------------
  
  
  ## Create the formula
  formula_covariates <- paste0(paste0(covariate_vars, collapse = " + "), " + ", interaction1_var, " * ", interaction2_var)
  f <- stats::as.formula(paste0(response_var, " ~ ", formula_covariates))
  
  
  ## Fit the logistic model
  regression_fit <- glm(formula = f, family = binomial(link = "logit"), data = data)
  regression_summ <- summary(regression_fit)
  
  
  # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
  # h(mm)
  
  
  # --------------------------------------------------------------------------
  ### Parse the regression summary for the interaction terms
  # --------------------------------------------------------------------------
  
  ## Generate data frame with coefficient names and levels and information about reference groups
  
  if(class(data[, interaction1_var]) %in% c("numeric", "integer") && class(data[, interaction2_var]) %in% c("numeric", "integer")){
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = "", reference1 = "", reference1_indx = 0, covariate2 = interaction2_var, levels2 = "", reference2 = "", reference2_indx = 0, 
      stringsAsFactors = FALSE)
    
  }else if(class(data[, interaction1_var]) %in% c("numeric", "integer") && class(data[, interaction2_var]) %in% c("factor")){
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## When the first level has zero counts, then the next level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction2_var])
    
    reference_indx <- which(tbl > 0)[1]
    
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = "", reference1 = "", reference1_indx = 0, covariate2 = interaction2_var, levels2 = levels(data[, interaction2_var]), reference2 = names(reference_indx), reference2_indx = as.numeric(reference_indx),
      stringsAsFactors = FALSE)
    
  }else if(class(data[, interaction1_var]) %in% c("factor") && class(data[, interaction2_var]) %in% c("numeric", "integer")){
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## When the first level has zero counts, then the next level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction1_var])
    
    reference_indx <- which(tbl > 0)[1]
    
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = levels(data[, interaction1_var]), reference1 = names(reference_indx), reference1_indx = as.numeric(reference_indx), covariate2 = interaction2_var, levels2 = "", reference2 = "", reference2_indx = 0,
      stringsAsFactors = FALSE)
    
  }else{
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## When the first level has zero counts, then the next level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction1_var])
    
    reference1_indx <- which(tbl > 0)[1]
    
    tbl <- table(data[, interaction2_var])
    
    reference2_indx <- which(tbl > 0)[1]
    
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = rep(levels(data[, interaction1_var]), times = nlevels(data[, interaction2_var])), reference1 = names(reference1_indx), reference1_indx = as.numeric(reference1_indx), covariate2 = interaction2_var, levels2 = rep(levels(data[, interaction2_var]), each = nlevels(data[, interaction1_var])), reference2 = names(reference2_indx), reference2_indx = as.numeric(reference2_indx),
      stringsAsFactors = FALSE)
    
    
  }
  
  
  coef_info <- out
  
  coef_info$coefficient <- paste0(coef_info$covariate1, coef_info$levels1, ":", coef_info$covariate2, coef_info$levels2)
  
  rownames(coef_info) <- coef_info$coefficient
  
  
  # --------------------------------------------------------------------------
  ## Calculate confidence intervals
  # --------------------------------------------------------------------------
  
  
  ## Note: When MASS is loaded, the 95% CIs are calculated (using confint) based on profile likelihood
  ## To compute 95% Wald CIs (based on asymptotic normality), one needs to use confint.default()
  
  
  ## There can be an error from confint when all samples have the same response
  # confint_res <- NULL
  # try(confint_res <- exp(stats::confint(regression_fit)), silent = TRUE)
  # if(is.null(confint_res)){
  #   confint_res <- exp(stats::confint.default(regression_fit))
  # }
  
  ## I use the Wald CIs because they are in concordance with p-values i.e. they contain 1 when p-value is not significant
  confint_res <- exp(stats::confint.default(regression_fit))
  
  
  # --------------------------------------------------------------------------
  ## Append results from regression
  # --------------------------------------------------------------------------
  
  conf_int <- data.frame(coefficient = rownames(confint_res), confint_res[, c("2.5 %", "97.5 %"), drop = FALSE], stringsAsFactors = FALSE)
  colnames(conf_int) <- c("coefficient", "CI95_lower", "CI95_upper")
  
  coefficients <- data.frame(coefficient = rownames(regression_summ$coefficients), regression_summ$coefficients[, c("Estimate", "Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE)
  colnames(coefficients) <- c("coefficient", "OR", "pvalue")
  coefficients$OR <- exp(coefficients$OR)
  
  coef_info$n <- nrow(data) - length(regression_summ$na.action)
  
  coef_info <- coef_info %>% 
    dplyr::left_join(conf_int, by = "coefficient") %>% 
    dplyr::left_join(coefficients, by = "coefficient")
  
  
  ## Calculate adjusted p-values using the Benjamini & Hochberg method
  coef_info$adj_pvalue <- stats::p.adjust(coef_info$pvalue, method = "BH")
  
  
  # --------------------------------------------------------------------------
  ### Return results 
  # --------------------------------------------------------------------------
  
  
  res <- coef_info[coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
  
  ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA because eventually this level is not used as a reference in the fitted model.
  if(any(c(res$reference1_indx > 1, res$reference2_indx > 1))){
    
    res[, -which(colnames(res) %in% c("covariate1", "levels1", "reference1", "reference1_indx", "covariate2", "levels2", "reference2", "reference2_indx", "coefficient", "n"))] <- NA
    
  }
  
  
  
  # --------------------------------------------------------------------------
  ### Prepare the output data frame that will be displayed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate1 = variable_names[res$covariate1], 
    Effect1 = format_vs(res$levels1, res$reference1),
    Covariate2 = variable_names[res$covariate2], 
    Effect2 = format_vs(res$levels2, res$reference2),
    `Total n` = as.character(res$n),
    `OR` = as.character(round(res$OR, 2)),
    `OR 95% CI` = format_CIs(res$CI95_lower, res$CI95_upper),
    `P-value` = format_pvalues(res$pvalue),
    `Adj. P-value` = format_pvalues(res$adj_pvalue),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  if(!print_adjpvalues){
    out$`Adj. P-value` <- NULL
  }
  
  ### For all numerical covariate Effect1 and Effect2 etc. are empty and we do not display them.
  if(all(out$Effect1 == "")){
    out$Effect1 <- NULL
  }
  
  if(all(out$Effect2 == "")){
    out$Effect2 <- NULL
  }
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Effect of interaction between ", variable_names[interaction1_var], " and ", variable_names[interaction2_var], " on ", variable_names[response_var], ".")
    
    if(!is.null(covariate_vars)){
      caption <- paste0(caption, paste0("Logistic regression model additionally includes ", paste0(variable_names[covariate_vars], collapse = ", "), ". "))
    }
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}




#' @rdname wrapper_logistic_regression_core_interaction
#' 
#' @inheritParams wrapper_logistic_regression_core_interaction
#' @param strat1_var Name of the first stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_logistic_regression_core_interaction_strat <- function(data, response_var, interaction1_var, interaction2_var, covariate_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
  
  ### Covariates cannot include strata
  stopifnot(length(intersect(c(strat1_var, strat2_var), c(interaction1_var, interaction2_var, covariate_vars))) == 0)
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(response_var, interaction1_var, interaction2_var, covariate_vars, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  

  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  wrapper_res <- lapply(1:length(strata2_levels), function(j){
    # j = 3
    
    data_strata2 <- data[data[, strat2_var] == strata2_levels[j] & !is.na(data[, strat2_var]), ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    wrapper_res <- lapply(1:length(strata1_levels), function(i){
      # i = 1
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      wrapper_res <- wrapper_logistic_regression_core_interaction(data = data_strata1, response_var = response_var, interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars, variable_names = variable_names, caption = caption, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      bresults(wrapper_res) <- res
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      boutput(wrapper_res) <- out
      
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
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(stats::p.adjust(res$pvalue, method = "BH"))
  }
  
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



#' Logistic regression estimating effect of interaction between biomarker and treatment
#' 
#' @inheritParams wrapper_logistic_regression_core_interaction_strat
#' @param treatment_var Name of the variable with treatment information.
#' @param biomarker_vars Vector of biomarker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
#' @export
wrapper_logistic_regression_interaction <- function(data, response_var, treatment_var, biomarker_vars, adjustment_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include treatment
  stopifnot(length(intersect(treatment_var, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(treatment_var, biomarker_vars, adjustment_vars), c(strat1_var, strat2_var))) == 0)
  
  
  vars_class <- sapply(data[, c(treatment_var)], class)
  stopifnot(all(vars_class == "factor"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 1
    
    interaction1_var <- biomarker_vars[i]
    interaction2_var <- treatment_var
    covariate_vars <- adjustment_vars
    
    wrapper_res <- wrapper_logistic_regression_core_interaction_strat(data = data, response_var = response_var,  interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(stats::p.adjust(res$pvalue, method = "BH"))
  }
  
  
  
  ### Rename 'Covariate' column name to 'Biomarker'
  
  colnames(res)[colnames(res) == "covariate1"] <- "biomarker"
  colnames(out)[colnames(out) == "Covariate1"] <- "Biomarker"
  
  colnames(res)[colnames(res) == "covariate2"] <- "treatment"
  out$Covariate2 <- NULL
  
  
  colnames(out)[colnames(out) == "Effect1"] <- "Biomarker Effect"
  colnames(out)[colnames(out) == "Effect2"] <- "Treatment Effect"
  
  
  ### Replace NAs with "" for columns that are missing for numerical biomarkers
  
  out$`Biomarker Effect`[is.na(out$`Biomarker Effect`)] <- ""
  
  
  
  ### Generate caption
  
  
  if(is.null(caption)){
    
    caption <- paste0("Effect of interaction between ", "biomarker", " and ", "treatment", " on ", variable_names[response_var], ".")
    
    
    if(is.null(adjustment_vars)){
      
      caption <- paste0(caption, "Unadjusted analysis. Logistic regression model includes only the biomarker and treatment.")
      
    }else{
      
      caption <- paste0(caption, "Adjusted analysis. Logistic regression model includes the biomarker, treatment and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}









