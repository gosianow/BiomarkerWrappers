



#' Cox regression with simple additive model
#' 
#' This function can be used if one is interested in effects of multiple covariates in a multivariate model.
#' 
#' @param data Data frame.
#' @param tte_var Name of the time-to-event variable. This variable must be numeric.
#' @param censor_var Name of the censor variable. It has to be numeric and encode 1 for event and 0 for censor.
#' @param covariate_vars Vector with names of covariates that are included in the formula of the simple additive model.
#' @param strata_vars Vector with names of covariates that are used as strata.
#' @param return_vars Vector with names of covariates for which the statistics should be returned. If NULL, statistics are returned for all covariates.
#' @param variable_names Named vector with variable names. If not supplied, variable names are created by replacing in column names underscores with spaces.
#' @param caption Caption for the table with results.
#' @param force_empty_cols Logical. Whether to display output columns which are all empty.
#' @param print_mst Logical. Whether to print median survival time (MST).
#' @param print_total Logical. Whether to print total number of samples and total number of events.
#' @param print_pvalues Logical. Whether to print p-values.
#' @param print_adjpvalues Logical. Whether to print adjusted p-values.
#' @details 
#' If for a factor covariate that should be returned the reference level has zero count, results are set to NAs because this levels is not used as a reference which means that it is not possible to fit the model that we want.
#' @examples
#' 
#' data(bdata)
#' 
#' data <- bdata
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' covariate_vars <- c("Treatment_Arm", "GeneA", "IPI", "Cell_Of_Origin")
#' 
#' x <- wrapper_core_cox_regression_simple(data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars)
#' 
#' boutput(x)
#' 
#' bforest(x)
#' 
#' 
#' ### Fit a stratified model 
#' 
#' covariate_vars <- c("Treatment_Arm", "GeneA")
#' strata_vars <- c("IPI", "Cell_Of_Origin")
#' 
#' x <- wrapper_core_cox_regression_simple(data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, strata_vars = strata_vars)
#' 
#' boutput(x)
#' 
#' 
#' @export
wrapper_core_cox_regression_simple <- function(data, tte_var, censor_var, covariate_vars, strata_vars = NULL, return_vars = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_mst = TRUE, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  ## Time to event variable must be numeric
  stopifnot(length(tte_var) == 1)
  stopifnot(is.numeric(data[, tte_var]))
  
  ## Censor variable must be numeric and encode 1 for event and 0 for censor
  stopifnot(length(censor_var) == 1)
  stopifnot(is.numeric(data[, censor_var]) && all(data[, censor_var] %in% c(0, 1)))
  
  
  covariate_class <- sapply(data[, covariate_vars], class)
  stopifnot(all(covariate_class %in% c("factor", "numeric", "integer")))
  
  
  if(!is.null(strata_vars)){
    strata_class <- sapply(data[, strata_vars], class)
    stopifnot(all(strata_class %in% c("factor")))
  }
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_vars, strata_vars)]), ]
  
  stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate data frame with coefficient names and levels and information about reference groups
  # --------------------------------------------------------------------------
  
  coef_info <- lapply(1:length(covariate_vars), function(i){
    # i = 1
    
    if(covariate_class[i] %in% c("numeric", "integer")){
      
      out <- data.frame(covariate = covariate_vars[i], covariate_class = covariate_class[i], subgroup = "", reference = "", reference_indx = 0, n = NA, nevent = NA, propevent = NA, MST = NA, MST_CI95_lower = NA, MST_CI95_upper = NA, HR_non_empty = TRUE,
        stringsAsFactors = FALSE)
      
    }else{
      
      # -----------------------------------
      # Number of events
      # -----------------------------------
      
      ## We use this way of calculating number of events because survfit does not return results for levels with zero counts.
      
      ## Check if the first level has non zero counts so it can be used as a reference group in the regression
      ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
      tbl <- table(data[, covariate_vars[i]])
      
      if(tbl[1] > 0){
        reference_indx <- 1
        names(reference_indx) <- names(tbl[1])
      }else{
        reference_indx <- rev(which(tbl > 0))[1]
      }
      
      ## Calculate nevent 
      tbl_event <- table(data[data[, censor_var] == 1, covariate_vars[i]])
      prop_event <- tbl_event / tbl * 100
      ## Replace NaN with NA
      prop_event[is.na(prop_event)] <- NA


      out <- data.frame(covariate = covariate_vars[i], covariate_class = covariate_class[i], subgroup = levels(data[, covariate_vars[i]]), reference = names(reference_indx), reference_indx = as.numeric(reference_indx), n = as.numeric(tbl), nevent = as.numeric(tbl_event), propevent = as.numeric(prop_event),
        stringsAsFactors = FALSE)
      
      out$HR_non_empty <- c(FALSE, rep(TRUE, nrow(out) - 1))
      
      out$coefficient <- paste0(out$covariate, "=", out$subgroup)
      
      
      # -----------------------------------
      # MST
      # -----------------------------------
      
      ## Use the survfit function to calculate n.start, events and median when the biomarker variable is categorical
      ## ?print.survfit The median and its confidence interval are defined by drawing a horizontal line at 0.5 on the plot of the survival curve and its confidence bands. The intersection of the line with the lower CI band defines the lower limit for the median's interval, and similarly for the upper band. If any of the intersections is not a point the we use the center of the intersection interval, e.g., if the survival curve were exactly equal to 0.5 over an interval. When data is uncensored this agrees with the usual definition of a median.
      ## Using conf.type = "plain" to obtain the same results as computed by Biostats.
      
      f <- stats::as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", covariate_vars[i]))
      
      survfit_fit <- survival::survfit(f, data, conf.type = "plain")
      survfit_summ <- summary(survfit_fit)
      
      
      if(is.matrix(survfit_summ$table)){
        survfit_out <- as.data.frame.matrix(survfit_summ$table[, c("median", "0.95LCL", "0.95UCL")])
      }else{
        survfit_out <- data.frame(MST = survfit_summ$table["median"],
          MST_CI95_lower = survfit_summ$table["0.95LCL"],
          MST_CI95_upper = survfit_summ$table["0.95UCL"], row.names = out$coefficient[out$n > 0])
      }
      
      
      colnames(survfit_out) <- c("MST", "MST_CI95_lower", "MST_CI95_upper")
      survfit_out$coefficient <- rownames(survfit_out)
      
      
      out <- dplyr::left_join(out, survfit_out, by = "coefficient")
      out$coefficient <- NULL
      
    }
    
    
    out
    
    
  })
  
  coef_info <- plyr::rbind.fill(coef_info)
  
  
  coef_info$n_total <- nrow(data)
  coef_info$nevent_total <- sum(data[, censor_var])
  
  coef_info$coefficient <- paste0(coef_info$covariate, coef_info$subgroup)
  
  rownames(coef_info) <- coef_info$coefficient
  
  
  # --------------------------------------------------------------------------
  # Cox regression
  # --------------------------------------------------------------------------
  
  ## Create the formula
  
  formula_covariates <- paste0(covariate_vars, collapse = " + ")
  
  if(!is.null(strata_vars)){
    formula_strata <- paste0("strata(", paste0(strata_vars, collapse = ", "), ")")
    formula_model <- paste0(formula_covariates, " + ", formula_strata)
  }else{
    formula_model <- formula_covariates
  }
  
  
  f <- stats::as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_model))
  
  
  ## Fit the Cox model
  regression_fit <- NULL
  
  try(regression_fit <- survival::coxph(f, data), silent = TRUE)
  
  if(is.null(regression_fit)){
    regression_summ <- NULL
  }else{
    regression_summ <- summary(regression_fit)
  }
  
  
  # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
  # h(mm)
  
  
  # --------------------------------------------------------------------------
  # Parce the regression summary
  # --------------------------------------------------------------------------
  
  
  if(is.null(regression_summ)){
    
    conf_int <- data.frame(coefficient = coef_info$coefficient, HR = NA, HR_CI95_lower = NA, HR_CI95_upper = NA, stringsAsFactors = FALSE)
    coefficients <- data.frame(coefficient = coef_info$coefficient, pvalue = NA, stringsAsFactors = FALSE)
    
  }else{
    
    conf_int <- data.frame(coefficient = rownames(regression_summ$conf.int), regression_summ$conf.int[, c("exp(coef)", "lower .95", "upper .95"), drop = FALSE], stringsAsFactors = FALSE)
    colnames(conf_int) <- c("coefficient", "HR", "HR_CI95_lower", "HR_CI95_upper")
    
    coefficients <- data.frame(coefficient = rownames(regression_summ$coefficients), regression_summ$coefficients[, c("Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE)
    colnames(coefficients) <- c("coefficient", "pvalue")
    
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
    
    ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA becasue this level is not used as a reference in the fitted model.
    if(all(res$reference_indx > 1)){
      
      res[, colnames(res) %in% c("HR", "HR_CI95_lower", "HR_CI95_upper", "pvalue", "adj_pvalue")] <- NA
      
    }
    
    
  }
  
  
  # --------------------------------------------------------------------------
  # Prepare the output data frame that will be dispalyed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = as.character(res$subgroup),
    `Total N` = as.character(res$n_total),
    `Total Events` = as.character(res$nevent_total),
    `N` = format_difference(res$n, digits = 0),
    `Events` = format_counts_and_props_core(counts = res$nevent, props = res$propevent, digits = 1),
    `MST` = format_or(res$MST, digits = 1, non_empty = res$covariate_class == "factor"),
    `MST 95% CI` = format_CIs(res$MST_CI95_lower, res$MST_CI95_upper, digits = 1, non_empty = res$covariate_class == "factor"),
    `HR` = format_or(res$HR, non_empty = res$HR_non_empty),
    `HR 95% CI` = format_CIs(res$HR_CI95_lower, res$HR_CI95_upper, non_empty = res$HR_non_empty),
    `P-value` = format_pvalues(res$pvalue, non_empty = res$HR_non_empty),
    `Adj. P-value` = format_pvalues(res$adj_pvalue, non_empty = res$HR_non_empty),
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
  
  if(!print_mst || !force_empty_cols){
    out$`MST` <- NULL
    out$`MST 95% CI` <- NULL
  }
  
  if(!print_total){
    col_total <- grep("^Total", colnames(out), value = TRUE)
    for(i in seq_along(col_total)){
      out[, col_total[i]] <- NULL
    }
  }
  
  if(!print_pvalues){
    out$`P-value` <- NULL
    out$`Adj. P-value` <- NULL
  }
  
  if(!print_adjpvalues){
    out$`Adj. P-value` <- NULL
  }
  

  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Covariate effect on ", variable_names[tte_var], ". ", 
      ifelse(is.null(strata_vars), "Unstratified ", "Stratified "),
      "Cox regression model includes: ", paste0(variable_names[covariate_vars], collapse = ", "), ".")
    
    if(!is.null(strata_vars)){
      caption <- paste0(caption, 
        " Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ".")
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}





#' @rdname wrapper_core_cox_regression_simple 
#' 
#' @param strat1_var Name of the firts stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' 
#' @examples 
#' 
#' data(bdata)
#' 
#' data <- bdata
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' covariate_vars <- c("IPI", "GeneA")
#' 
#' strat1_var = "Cell_Of_Origin"
#' strat2_var = "Treatment_Arm"
#' 
#' 
#' x <- wrapper_core_cox_regression_simple_strat(data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, strat1_var = strat1_var, strat2_var = strat2_var)
#' 
#' boutput(x)
#' 
#' @export
wrapper_core_cox_regression_simple_strat <- function(data, tte_var, censor_var, covariate_vars, strata_vars = NULL, return_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_mst = TRUE, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_vars, strat1_var, strat2_var, strata_vars)]), ]
  stopifnot(nrow(data) > 0)
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  wrapper_res <- lapply(1:length(strata2_levels), function(j){
    # j = 1
    
    data_strata2 <- data[data[, strat2_var] %in% strata2_levels[j], ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    wrapper_res <- lapply(1:length(strata1_levels), function(i){
      # i = 3
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] %in% strata1_levels[i], ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      wrapper_res <- wrapper_core_cox_regression_simple(data = data_strata1, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, strata_vars = strata_vars, return_vars = return_vars, variable_names = variable_names, caption = caption, force_empty_cols = force_empty_cols, print_mst = print_mst, print_total = print_total, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      
      
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











#' Cox regression estimating biomarker effect 
#' 
#' @inheritParams wrapper_core_cox_regression_simple_strat
#' @param biomarker_vars Vector of biomaker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
#' @export
wrapper_cox_regression_biomarker <- function(data, tte_var, censor_var, biomarker_vars, adjustment_vars = NULL, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_mst = TRUE, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include biomarker variables
  stopifnot(length(intersect(biomarker_vars, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(biomarker_vars, adjustment_vars), c(strata_vars, strat1_var, strat2_var))) == 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 2
    
    covariate_vars <- c(biomarker_vars[i], adjustment_vars)
    return_vars <- biomarker_vars[i]
    
    
    wrapper_res <- wrapper_core_cox_regression_simple_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, strata_vars = strata_vars, return_vars = return_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_mst = print_mst, print_total = print_total, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
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
    
    caption <- paste0("Biomarker effect on ", variable_names[tte_var], ". ")
    
    if(is.null(adjustment_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only the biomarker.")
      
    }else if(!is.null(adjustment_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes the biomarker and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }else if(is.null(adjustment_vars) && !is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, stratified analysis. Cox regression model includes the biomarker and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, stratified analysis. Cox regression model includes the biomarker and ", paste0(variable_names[adjustment_vars], collapse = ", "), " and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}






#' Cox regression estimating treatment effect within biomaker subgroups
#' 
#' @inheritParams wrapper_core_cox_regression_simple_strat
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector with names of categorical biomarkers. When NULL, overall treatment effect is estimated. 
#' @param adjustment_vars Vector of covariate names used for adjustment.
#' @export
wrapper_cox_regression_treatment <- function(data, tte_var, censor_var, treatment_var, adjustment_vars = NULL, strata_vars = NULL, biomarker_vars = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_mst = TRUE, print_total = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
    
    
    wrapper_res <- wrapper_core_cox_regression_simple_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, strata_vars = strata_vars, return_vars = return_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_mst = print_mst, print_total = print_total, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
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
    
    caption <- paste0("Treatment effect on ", variable_names[tte_var], ". ")
    
    if(is.null(adjustment_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only the treatment.")
      
    }else if(!is.null(adjustment_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes the treatment and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }else if(is.null(adjustment_vars) && !is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, stratified analysis. Cox regression model includes the treatment and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, stratified analysis. Cox regression model includes the treatment and ", paste0(variable_names[adjustment_vars], collapse = ", "), " and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}



###############################################################################
# Model with interaction
###############################################################################




#' Cox regression with additive model with interaction
#' 
#' @inheritParams wrapper_core_cox_regression_simple
#' @param interaction1_var Name of the first interaction variable.
#' @param interaction2_var Name of the second interaction variable.
#' 
#' @examples 
#' 
#' data(bdata)
#' 
#' data <- bdata
#' 
#' data$GeneA_cat2 <- cut_core_2groups(data$GeneA)
#' 
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' 
#' interaction1_var <- "GeneA_cat2"
#' interaction2_var <- "Treatment_Arm"
#' covariate_vars <- c("IPI", "Cell_Of_Origin")
#' 
#' x <- wrapper_core_cox_regression_interaction(data, tte_var = tte_var, censor_var = censor_var, interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars)
#' 
#' boutput(x)
#' 
#' @export
wrapper_core_cox_regression_interaction <- function(data, tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars = NULL, strata_vars = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  ## Time to event variable must be numeric
  stopifnot(length(tte_var) == 1)
  stopifnot(is.numeric(data[, tte_var]))
  
  ## Censor variable must be numeric and encode 1 for event and 0 for censor
  stopifnot(length(censor_var) == 1)
  stopifnot(is.numeric(data[, censor_var]) && all(data[, censor_var] %in% c(0, 1)))
  
  
  covariate_class <- sapply(data[, c(interaction1_var, interaction2_var, covariate_vars)], class)
  stopifnot(all(covariate_class %in% c("factor", "numeric", "integer")))
  
  
  if(!is.null(strata_vars)){
    strata_class <- sapply(data[, strata_vars], class)
    stopifnot(all(strata_class %in% c("factor")))
  }
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars, strata_vars)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  # --------------------------------------------------------------------------
  # Cox regression
  # --------------------------------------------------------------------------
  
  
  ## Create the formula
  
  formula_covariates <- paste0(paste0(covariate_vars, collapse = " + "), " + ", interaction1_var, " * ", interaction2_var)
  
  if(!is.null(strata_vars)){
    formula_strata <- paste0("strata(", paste0(strata_vars, collapse = ", "), ")")
    formula_model <- paste0(formula_covariates, " + ", formula_strata)
  }else{
    formula_model <- formula_covariates
  }
  
  
  f <- stats::as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_model))
  
  
  ## Fit the Cox model
  regression_fit <- survival::coxph(f, data)
  regression_summ <- summary(regression_fit)
  
  
  # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
  # h(mm)
  
  
  # --------------------------------------------------------------------------
  ### Parce the regression summary for the interaction terms
  # --------------------------------------------------------------------------
  
  ## Generate data frame with coefficient names and levels and information about reference groups
  
  if(class(data[, interaction1_var]) %in% c("numeric", "integer") && class(data[, interaction2_var]) %in% c("numeric", "integer")){
    
    out <- data.frame(covariate1 = interaction1_var, subgroup1 = "", reference1 = "", reference1_indx = 0, covariate2 = interaction2_var, subgroup2 = "", reference2 = "", reference2_indx = 0, 
      stringsAsFactors = FALSE)
    
  }else if(class(data[, interaction1_var]) %in% c("numeric", "integer") && class(data[, interaction2_var]) %in% c("factor")){
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction2_var])
    
    if(tbl[1] > 0){
      reference_indx <- 1
      names(reference_indx) <- names(tbl[1])
    }else{
      reference_indx <- rev(which(tbl > 0))[1]
    }
    
    
    out <- data.frame(covariate1 = interaction1_var, subgroup1 = "", reference1 = "", reference1_indx = 0, covariate2 = interaction2_var, subgroup2 = levels(data[, interaction2_var]), reference2 = names(reference_indx), reference2_indx = as.numeric(reference_indx),
      stringsAsFactors = FALSE)
    
  }else if(class(data[, interaction1_var]) %in% c("factor") && class(data[, interaction2_var]) %in% c("numeric", "integer")){
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction1_var])
    
    if(tbl[1] > 0){
      reference_indx <- 1
      names(reference_indx) <- names(tbl[1])
    }else{
      reference_indx <- rev(which(tbl > 0))[1]
    }
    
    
    out <- data.frame(covariate1 = interaction1_var, subgroup1 = levels(data[, interaction1_var]), reference1 = names(reference_indx), reference1_indx = as.numeric(reference_indx), covariate2 = interaction2_var, subgroup2 = "", reference2 = "", reference2_indx = 0,
      stringsAsFactors = FALSE)
    
  }else{
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction1_var])
    
    if(tbl[1] > 0){
      reference1_indx <- 1
      names(reference1_indx) <- names(tbl[1])
    }else{
      reference1_indx <- rev(which(tbl > 0))[1]
    }
    
    tbl <- table(data[, interaction2_var])
    
    if(tbl[1] > 0){
      reference2_indx <- 1
      names(reference2_indx) <- names(tbl[1])
    }else{
      reference2_indx <- rev(which(tbl > 0))[1]
    }
    
    
    out <- data.frame(covariate1 = interaction1_var, subgroup1 = rep(levels(data[, interaction1_var]), times = nlevels(data[, interaction2_var])), reference1 = names(reference1_indx), reference1_indx = as.numeric(reference1_indx), covariate2 = interaction2_var, subgroup2 = rep(levels(data[, interaction2_var]), each = nlevels(data[, interaction1_var])), reference2 = names(reference2_indx), reference2_indx = as.numeric(reference2_indx),
      stringsAsFactors = FALSE)
    
    
  }
  
  
  coef_info <- out
  
  coef_info$coefficient <- paste0(coef_info$covariate1, coef_info$subgroup1, ":", coef_info$covariate2, coef_info$subgroup2)
  
  rownames(coef_info) <- coef_info$coefficient
  
  # --------------------------------------------------------------------------
  # Append results from regression
  # --------------------------------------------------------------------------
  
  conf_int <- data.frame(coefficient = rownames(regression_summ$conf.int), regression_summ$conf.int[, c("exp(coef)", "lower .95", "upper .95"), drop = FALSE], stringsAsFactors = FALSE)
  colnames(conf_int) <- c("coefficient", "HR", "HR_CI95_lower", "HR_CI95_upper")
  
  coefficients <- data.frame(coefficient = rownames(regression_summ$coefficients), regression_summ$coefficients[, c("Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE)
  colnames(coefficients) <- c("coefficient", "pvalue")
  
  coef_info$n_total <- regression_summ$n
  
  coef_info <- coef_info %>% 
    dplyr::left_join(conf_int, by = "coefficient") %>% 
    dplyr::left_join(coefficients, by = "coefficient")
  
  
  ## Calculate adjusted p-values using the Benjamini & Hochberg method
  coef_info$adj_pvalue <- stats::p.adjust(coef_info$pvalue, method = "BH")
  
  
  # --------------------------------------------------------------------------
  # Return results 
  # --------------------------------------------------------------------------
  
  
  res <- coef_info[coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
  
  ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA becasue eventually this level is not used as a reference in the fitted model.
  if(any(c(res$reference1_indx > 1, res$reference2_indx > 1))){
    
    res[, colnames(res) %in% c("HR", "HR_CI95_lower", "HR_CI95_upper", "pvalue", "adj_pvalue")] <- NA
    
  }
  
  
  
  # --------------------------------------------------------------------------
  ### Prepare the output data frame that will be dispalyed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate1 = variable_names[res$covariate1], 
    Effect1 = format_vs(res$subgroup1, res$reference1),
    Covariate2 = variable_names[res$covariate2], 
    Effect2 = format_vs(res$subgroup2, res$reference2),
    `Total N` = as.character(res$n_total),
    `HR` = format_or(res$HR),
    `HR 95% CI` = format_CIs(res$HR_CI95_lower, res$HR_CI95_upper),
    `P-value` = format_pvalues(res$pvalue),
    `Adj. P-value` = format_pvalues(res$adj_pvalue),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
    out$`Adj. P-value` <- NULL
  }
  
  if(!print_adjpvalues){
    out$`Adj. P-value` <- NULL
  }
  
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  if(is.null(caption)){
    
    caption <- paste0("Effect of interaction between ", variable_names[interaction1_var], " and ", variable_names[interaction2_var], " on ", variable_names[tte_var], ". ")
    
    if(is.null(covariate_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only: ", paste0(variable_names[c(interaction1_var, interaction2_var)], collapse = ", "), ".")
      
    }else if(!is.null(covariate_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes: ", paste0(variable_names[c(interaction1_var, interaction2_var, covariate_vars)], collapse = ", "), ". ")
      
    }else if(is.null(covariate_vars) && !is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, stratified analysis. Cox regression model includes: ", paste0(variable_names[c(interaction1_var, interaction2_var)], collapse = ", "), " and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, stratified analysis. Cox regression model includes: ", paste0(variable_names[c(interaction1_var, interaction2_var, covariate_vars)], collapse = ", "), " and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}




#' @rdname wrapper_core_cox_regression_interaction
#' 
#' @inheritParams wrapper_core_cox_regression_interaction
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_core_cox_regression_interaction_strat <- function(data, tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars = NULL, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
  stopifnot(length(intersect(c(strat1_var, strat2_var), c(interaction1_var, interaction2_var, covariate_vars, strata_vars))) == 0)
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars, strata_vars, strat1_var, strat2_var)]), ]
  
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
      
      
      wrapper_res <- wrapper_core_cox_regression_interaction(data = data_strata1, tte_var = tte_var, censor_var = censor_var, interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars, strata_vars = strata_vars, variable_names = variable_names, caption = caption, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      
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
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
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



#' Cox regression estimating interaction effect between biomaker and treatment
#' 
#' @inheritParams wrapper_core_cox_regression_interaction_strat
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector of biomaker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
#' @export
wrapper_cox_regression_interaction <- function(data, tte_var, censor_var, treatment_var, biomarker_vars, adjustment_vars = NULL, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL,  print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include treatment
  stopifnot(length(intersect(treatment_var, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(treatment_var, biomarker_vars, adjustment_vars, strata_vars), c(strat1_var, strat2_var))) == 0)
  
  
  vars_class <- sapply(data[, c(treatment_var)], class)
  stopifnot(all(vars_class == "factor"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  covariate_vars <- adjustment_vars
  
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 1
    
    interaction1_var <- biomarker_vars[i]
    interaction2_var <- treatment_var
    
    
    wrapper_res <- wrapper_core_cox_regression_interaction_strat(data = data, tte_var = tte_var, censor_var = censor_var, interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
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
  
  
  
  ### Rename 'Covariate' column name to 'Biomarker'
  
  colnames(res)[colnames(res) == "covariate1"] <- "biomarker"
  colnames(out)[colnames(out) == "Covariate1"] <- "Biomarker"
  
  colnames(res)[colnames(res) == "covariate2"] <- "treatment"
  out$Covariate2 <- NULL
  
  
  colnames(out)[colnames(out) == "Effect1"] <- "Biomarker Effect"
  colnames(out)[colnames(out) == "Effect2"] <- "Treatment Effect"
  
  
  ### Replace NAs with "" for columns that are missing for numerical biomarkers
  if("Biomarker Effect" %in% colnames(out)){
    out$`Biomarker Effect`[is.na(out$`Biomarker Effect`)] <- ""
  }
  
  
  
  
  ### Generate caption
  
  if(is.null(caption)){
    
    caption <- paste0("Effect of interaction between ", "biomarker", " and ", "treatment", " on ", variable_names[tte_var], ". ")
    
    if(is.null(covariate_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only: ", "biomarker and treatment", ".")
      
    }else if(!is.null(covariate_vars) && is.null(strata_vars)){
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes: biomarker,  treatment and ", paste0(variable_names[c(covariate_vars)], collapse = ", "), ". ")
      
    }else if(is.null(covariate_vars) && !is.null(strata_vars)){
      
      caption <- paste0(caption, "Unadjusted, stratified analysis. Cox regression model includes: ", "biomarker,  treatment", " and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, stratified analysis. Cox regression model includes: biomarker,  treatment and ", paste0(variable_names[c(covariate_vars)], collapse = ", "), " and stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}









