




#' Cox regression with simple additive model for Nested Case-Control NCC study 
#' 
#' @param data Data frame preprocessed for the NCC analysis with 'multipleNCC::wpl()'.
#' @param tte_var Name of the time-to-event variable. This variable must be numeric.
#' @param censor_var Name of the censor variable. It has to be numeric and encode 1 for event and 0 for censor.
#' @param covariate_vars Vector with names of covariates that are included in the formula of the simple additive model.
#' @param ncc_vars Vector of names of covariates that were measured in the NCC.
#' @param samplestat_var Name of variable indicating samplestat values. See 'multipleNCC::wpl()'.
#' @param m See 'multipleNCC::wpl()'.
#' @param match.var See 'multipleNCC::wpl()'. It has to be a matrix of continuous values. 
#' @param match.int See 'multipleNCC::wpl()'.
#' @param return_vars Vector with names of covariates for which the statistics should be returned. If NULL, statistics are returned for all covariates.
#' @param variable_names Named vector with variable names. If not supplied, variable names are created by replacing in column names underscores with spaces.
#' @param caption Caption for the table with results.
#' @param print_pvalues Logical. Whether to print p-values.
#' @details 
#' If for a factor covariate that should be returned the reference level has zero count, results are set to NAs because this levels is not used as a reference which means that it is not possible to fit the model that we want.
#' @export
wrapper_ncc_cox_regression_core_prognostic <- function(data, tte_var, censor_var, covariate_vars, ncc_vars, samplestat_var, m, match.var, match.int, return_vars = NULL, variable_names = NULL, caption = NULL, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  covariate_vars <- c(covariate_vars, ncc_vars)
  
  stopifnot(isValidAndUnreservedName(c(tte_var, censor_var, covariate_vars, samplestat_var, return_vars)))
  
  stopifnot(is.data.frame(data))
  
  stopifnot(type(match.var) == "double")
  
  ## Time to event variable must be numeric
  stopifnot(length(tte_var) == 1)
  stopifnot(is.numeric(data[, tte_var]))
  
  ## Censor variable must be numeric and encode 1 for event and 0 for censor
  stopifnot(length(censor_var) == 1)
  stopifnot(is.numeric(data[, censor_var]) && all(data[, censor_var] %in% c(0, 1, NA)))
  
  
  covariate_class <- sapply(data[, covariate_vars], class)
  stopifnot(all(covariate_class %in% c("factor", "numeric", "integer")))
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_vars)]), , drop = FALSE]
  
  # stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate data frame with coefficient names and levels and information about reference groups
  # --------------------------------------------------------------------------
  
  
  coef_info <- lapply(1:length(covariate_vars), function(i){
    # i = 1
    
    if(covariate_class[i] %in% c("numeric", "integer")){
      
      out <- data.frame(covariate = covariate_vars[i], covariate_class = covariate_class[i], subgroup = "", reference = "", reference_indx = 0, n = NA, HR_non_empty = TRUE,
        stringsAsFactors = FALSE)
      
    }else{
      
      # -----------------------------------
      # Number of events
      # -----------------------------------
      
      ## We use this way of calculating number of events because survfit does not return results for levels with zero counts.
      
      ## Check if the first level has non zero counts so it can be used as a reference group in the regression
      ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
      
      
      if(covariate_vars[i] %in% ncc_vars){
        tbl <- table(data[as.logical(data[, samplestat_var]), covariate_vars[i]])
      }else{
        tbl <- table(data[, covariate_vars[i]])
      }
      
      
      if(tbl[1] > 0){
        reference_indx <- 1
        names(reference_indx) <- names(tbl[1])
      }else{
        reference_indx <- rev(which(tbl > 0))[1]
      }
      
      
      out <- data.frame(covariate = covariate_vars[i], covariate_class = covariate_class[i], subgroup = levels(data[, covariate_vars[i]]), reference = names(reference_indx), reference_indx = as.numeric(reference_indx), n = as.numeric(tbl), 
        stringsAsFactors = FALSE)
      
      out$HR_non_empty <- c(FALSE, rep(TRUE, nrow(out) - 1))
      
      out$coefficient <- paste0(out$covariate, "=", out$subgroup)
      
      
    }
    
    
    out
    
    
  })
  
  coef_info <- plyr::rbind.fill(coef_info)
  
  
  coef_info$n_total <- nrow(data)
  
  coef_info$n_total[coef_info$covariate %in% ncc_vars] <- nrow(data[as.logical(data[, samplestat_var]), ])
  
  coef_info$coefficient <- paste0(coef_info$covariate, coef_info$subgroup)
  
  rownames(coef_info) <- coef_info$coefficient
  
  
  # --------------------------------------------------------------------------
  # Cox regression with Weighted partial likelihood for nested case-control data
  # --------------------------------------------------------------------------
  
  if(nrow(data) > 0){
    
    ## Create the formula
    
    formula_covariates <- paste0(covariate_vars, collapse = " + ")
    
    formula_model <- formula_covariates
    
    f <- stats::as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_model))
    
    
    ## Fit the Cox model with WPL
 
    ipw_fit <- NULL
    
    try(ipw_fit <- multipleNCC::wpl(f,
      data = data, 
      samplestat = data[, samplestat_var], 
      m = m, 
      weight.method = "KM",
      match.var = match.var,
      match.int = match.int), silent = TRUE)
    
    
    if(is.null(ipw_fit)){
      ipw_summ <- NULL
    }else{
      class(ipw_fit) <- "coxph"
      ipw_fit$score <- Inf
      ipw_summ <- summary(ipw_fit)
    }
    
    
    # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
    # h(mm)
    
    
  }else{
    
    ipw_summ <- NULL
    
  }
  
  
  # --------------------------------------------------------------------------
  # Parse the regression summary
  # --------------------------------------------------------------------------
  
  
  if(is.null(ipw_summ)){
    
    conf_int <- data.frame(coefficient = coef_info$coefficient, HR = NA, HR_CI95_lower = NA, HR_CI95_upper = NA, stringsAsFactors = FALSE)
    coefficients <- data.frame(coefficient = coef_info$coefficient, pvalue = NA, stringsAsFactors = FALSE)
    
  }else{
    
    conf_int <- data.frame(coefficient = rownames(ipw_summ$conf.int), ipw_summ$conf.int[, c("exp(coef)", "lower .95", "upper .95"), drop = FALSE], stringsAsFactors = FALSE, row.names = NULL)
    colnames(conf_int) <- c("coefficient", "HR", "HR_CI95_lower", "HR_CI95_upper")
    
    coefficients <- data.frame(coefficient = rownames(ipw_summ$coefficients), ipw_summ$coefficients[, c("Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE, row.names = NULL)
    colnames(coefficients) <- c("coefficient", "pvalue")
    
    ### Somehow there is 'x' prefix added 
    conf_int$coefficient <- gsub("^x", "", conf_int$coefficient)
    coefficients$coefficient <- gsub("^x", "", coefficients$coefficient)
    
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
    
    ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA because this level is not used as a reference in the fitted model.
    if(any(res$reference_indx > 1 | is.na(res$reference_indx))){
      
      res[, colnames(res) %in% c("HR", "HR_CI95_lower", "HR_CI95_upper", "pvalue", "adj_pvalue")] <- NA
      
    }
    
    
  }
  
  
  # --------------------------------------------------------------------------
  # Prepare the output data frame that will be displayed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  
  out1 <- data.frame(
    Covariate = variable_names[res$covariate], 
    Subgroup = as.character(res$subgroup),
    `Total N` = as.character(res$n_total),
    `N` = format_difference(res$n, digits = 0),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  out3 <- data.frame(
    `HR` = format_or(res$HR, non_empty = res$HR_non_empty),
    `HR 95% CI` = format_CIs(res$HR_CI95_lower, res$HR_CI95_upper, non_empty = res$HR_non_empty),
    `P-value` = format_pvalues(res$pvalue, non_empty = res$HR_non_empty),
    `Adj. P-value` = format_pvalues(res$adj_pvalue, non_empty = res$HR_non_empty),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  out <- cbind(out1, out3) 
  
  
  stopifnot(all(sapply(out, class) == "character"))
  

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
      "NCC Cox regression model includes: ", paste0(variable_names[covariate_vars], collapse = ", "), ".")
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}





#' Cox regression with additive model with interactions for Nested Case-Control NCC study
#' 
#' Returns HRs with treatment effects in the biomarker subgroups and p-value of the interaction term
#' 
#' @param data Data frame preprocessed for the NCC analysis with 'multipleNCC::wpl()'.
#' @param tte_var Name of the time-to-event variable. This variable must be numeric.
#' @param censor_var Name of the censor variable. It has to be numeric and encode 1 for event and 0 for censor.
#' @param interaction1_var Names of the first interaction variable. It would correspond to the biomarker. This variable must be a factor with two levels.
#' @param interaction2_var Name of the second interaction variable. It would correspond to the treatment arm. This variable must be a factor with two levels.
#' @param covariate_vars Vector with names of covariates that are included in the formula of the additive model.
#' @param ncc_vars Vector of names of covariates that were measured in the NCC.
#' @param samplestat_var Name of variable indicating samplestat values. See 'multipleNCC::wpl()'.
#' @param m See 'multipleNCC::wpl()'.
#' @param match.var See 'multipleNCC::wpl()'. It has to be a matrix of continuous values. 
#' @param match.int See 'multipleNCC::wpl()'.
#' @param variable_names Named vector with variable names. If not supplied, variable names are created by replacing in column names underscores with spaces.
#' @param caption Caption for the table with results.
#' @param print_pvalues Logical. Whether to print p-values.
#' @details 
#' If for a factor covariate that should be returned the reference level has zero count, results are set to NAs because this levels is not used as a reference which means that it is not possible to fit the model that we want.
#' @export
wrapper_ncc_cox_regression_core_predictive <- function(data, tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars = NULL, ncc_vars = NULL, samplestat_var, m, match.var, match.int, variable_names = NULL, caption = NULL, print_pvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  covariate_vars <- c(covariate_vars, ncc_vars)
  
  stopifnot(isValidAndUnreservedName(c(tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars, samplestat_var)))
  
  stopifnot(is.data.frame(data))
  
  stopifnot(type(match.var) == "double")
  
  ## Time to event variable must be numeric
  stopifnot(length(tte_var) == 1)
  stopifnot(is.numeric(data[, tte_var]))
  
  ## Censor variable must be numeric and encode 1 for event and 0 for censor
  stopifnot(length(censor_var) == 1)
  stopifnot(is.numeric(data[, censor_var]) && all(data[, censor_var] %in% c(0, 1, NA)))
  
  
  covariate_class <- sapply(data[, covariate_vars], class)
  stopifnot(all(covariate_class %in% c("factor", "numeric", "integer")))
  
  
  ### ARM and BIOMARKER must have two levels
  
  covariate_class <- sapply(data[, c(interaction1_var, interaction2_var)], class)
  stopifnot(all(covariate_class %in% c("factor")))
  
  stopifnot(all(sapply(data[, c(interaction1_var, interaction2_var)], nlevels) == 2))
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars)]), , drop = FALSE]
  
  # stopifnot(nrow(data) > 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate data frame with coefficient names and levels and information about reference groups
  # In the function we need that only for the BIOMARKER because we return treatment effects within biomarker subgroups
  # --------------------------------------------------------------------------
  
  
  # -----------------------------------
  # Number of samples 
  # -----------------------------------
  
  ## We use this way of calculating number of events because survfit does not return results for levels with zero counts.
  
  ## Check if the first level has non zero counts so it can be used as a reference group in the regression
  ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
  
  tbl <- table(data[as.logical(data[, samplestat_var]), interaction1_var])
  
  
  if(tbl[1] > 0){
    reference_indx <- 1
    names(reference_indx) <- names(tbl[1])
  }else{
    reference_indx <- rev(which(tbl > 0))[1]
  }
  
  
  out <- data.frame(covariate = interaction1_var, covariate_class = "factor", subgroup = levels(data[, interaction1_var]), reference = names(reference_indx), reference_indx = as.numeric(reference_indx), n = as.numeric(tbl), 
    stringsAsFactors = FALSE)
  
  out$HR_non_empty <- rep(TRUE, nrow(out))
  
  out$coefficient <- paste0(out$covariate, "=", out$subgroup)
  
  
  coef_info <- out
  
  
  coef_info$n_total <- nrow(data[as.logical(data[, samplestat_var]), ])
  
  coef_info$coefficient <- paste0(coef_info$covariate, coef_info$subgroup)
  
  rownames(coef_info) <- coef_info$coefficient
  
  
  # --------------------------------------------------------------------------
  # Cox regression with Weighted partial likelihood for nested case-control data
  # --------------------------------------------------------------------------
  
  if(nrow(data) > 0){
    
    ## Create the formula following the code from Joe 
    
    formula_covariates <- paste0(paste0(c(interaction2_var, covariate_vars), collapse = " + "), " + ", interaction2_var, "*", interaction1_var)
    
    formula_model <- formula_covariates
    
    f <- stats::as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_model))
    
    
    ## Fit the Cox model with WPL
    
    ipw_fit <- NULL
    
    try(ipw_fit <- multipleNCC::wpl(f,
      data = data, 
      samplestat = data[, samplestat_var], 
      m = m, 
      weight.method = "KM",
      match.var = match.var,
      match.int = match.int), silent = TRUE)
    
    
    # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
    # h(mm)
    
    
  }else{
    
    ipw_fit <- NULL
    
  }
  
  
  
  # --------------------------------------------------------------------------
  # Parse the regression summary
  # --------------------------------------------------------------------------
  
  
  if(is.null(ipw_fit)){
    
    conf_int <- data.frame(coefficient = coef_info$coefficient, HR = NA, HR_CI95_lower = NA, HR_CI95_upper = NA, stringsAsFactors = FALSE)
    coefficients <- data.frame(coefficient = coef_info$coefficient, pvalue = NA, stringsAsFactors = FALSE)
    
  }else{
    
    ### Code from Joe to compute the HR and CIs of the treatment effect in the biomarker subgroups 
    
    class(ipw_fit) <- "coxph"
    ipw_fit$score <- Inf
    
    ipw_summ <- summary(ipw_fit)
    
    CI <- 0.95
    z <- qnorm((1 + CI)/2, 0, 1)
    
    coef_indx <- length(ipw_fit$coefficients)
    
    vc <- vcov(ipw_fit)[1, coef_indx]
    old <- diag((vcov(ipw_fit)))[c(1, coef_indx)]
    old_var <- sum((old))
    se_new <- sqrt(old_var + 2*vc)
    
    est <- exp(coefficients(ipw_fit))
    low <- exp(coefficients(ipw_fit) - z*sqrt(diag(vcov(ipw_fit))))
    high <- exp(coefficients(ipw_fit) + z*sqrt(diag(vcov(ipw_fit))))
    
    cf <- coefficients(ipw_fit)[c(1, coef_indx)] %>% sum
    bm_hr <- exp(cf)
    bm_hr_low <- exp(cf - z*se_new)
    bm_hr_high <- exp(cf + z*se_new)
    
    # vals = rbind('WT'=c(est[c(1)],low[c(1)],high[c(1)]),'BM+'=c(bm_hr,bm_hr_low,bm_hr_high))
    # colnames(vals) = c("HR","Lower","Upper")
    # vals
    
    
    conf_int <- data.frame(coefficient = paste0(interaction1_var, levels(data[, interaction1_var])), 
      "exp(coef)" = c(est[1], bm_hr),
      "lower .95" = c(low[1], bm_hr_low),
      "upper .95" = c(high[1], bm_hr_high),
      stringsAsFactors = FALSE, row.names = NULL)
    
    colnames(conf_int) <- c("coefficient", "HR", "HR_CI95_lower", "HR_CI95_upper")
    
    
    coefficients <- data.frame(coefficient = paste0(interaction1_var, levels(data[, interaction1_var])),
      "Pr(>|z|)" = c(ipw_summ$coefficients[coef_indx, "Pr(>|z|)"], NA), 
      stringsAsFactors = FALSE, row.names = NULL)
    
    colnames(coefficients) <- c("coefficient", "pvalue")

    
  }
  
  
  # --------------------------------------------------------------------------
  # Append results from regression
  # --------------------------------------------------------------------------
  
  
  coef_info <- coef_info %>% 
    dplyr::left_join(conf_int, by = "coefficient") %>% 
    dplyr::left_join(coefficients, by = "coefficient")
  
  
  
  # --------------------------------------------------------------------------
  # Return results for all the covariates
  # --------------------------------------------------------------------------
  
  res <- coef_info
  
  # --------------------------------------------------------------------------
  # Prepare the output data frame that will be displayed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  
  out1 <- data.frame(
    Covariate = variable_names[res$covariate], 
    Subgroup = as.character(res$subgroup),
    `Total N` = as.character(res$n_total),
    `N` = format_difference(res$n, digits = 0),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  out3 <- data.frame(
    `HR` = format_or(res$HR, non_empty = res$HR_non_empty),
    `HR 95% CI` = format_CIs(res$HR_CI95_lower, res$HR_CI95_upper, non_empty = res$HR_non_empty),
    `P-value` = format_pvalues(res$pvalue, non_empty = 1),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  out <- cbind(out1, out3) 
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
    out$`Adj. P-value` <- NULL
  }
  
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Treatment effect on ", variable_names[tte_var], " within ", variable_names[interaction1_var], " subgroups. ", 
      "NCC Cox regression model includes: ", paste0(variable_names[c(covariate_vars, interaction1_var, interaction2_var)], collapse = ", "), " and interaction term between ", variable_names[interaction1_var], " and ", variable_names[interaction2_var], ".")
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}








#############################################################################
### Code from Joe
#############################################################################





### Prognostic from Joe 



# df_ncc = full_df[-which(full_df$fmiAssayed==0 & full_df$IDFS_EVNT==1), ]
# 
# df_ncc$assayed = as.numeric(df_ncc$sequencing & df_ncc$fmiAssayed)
# 
# 
# if(i==0) df_ncc = subset(df_ncc,ARM==0)
# if(i==1) df_ncc = subset(df_ncc,ARM==1)
# 
# fits <- lapply(which(colnames(df_ncc) %in% gg),function(j){
#   
#   print(colnames(df_ncc)[j])
#   df_ncc$val <- as.factor(df_ncc[,j])
#   
#   
#   ipw_fit <- multipleNCC::wpl(Surv(IDFS, event=IDFS_EVNT) ~ ARM +
#       hr_status + 
#       node_status + 
#       age + 
#       chemo_tx + val, 
#     
#     data = df_ncc, 
#     m = 3, 
#     
#     samplestat = df_ncc$assayed + df_ncc$IDFS_EVNT, 
#     
#     weight.method = "KM",
#     match.var = cbind(df_ncc$ARM, df_ncc$hr_status, df_ncc$node_status),
#     
#     match.int = rep(0,6))
#   
#   
#   class(ipw_fit) = "coxph"
#   ipw_fit$score = Inf
#   tidyres = tidy(ipw_fit,conf.int = TRUE,exponentiate = TRUE)[6:7,c(2,7,8)]
#   t(as.matrix(cbind(tidyres,p.value = tidy(ipw_fit)[6:7,6])))
# })
# fits <- Reduce("cbind",fits)
# colnames(fits) <- paste(rep(colnames(df_ncc)[which(colnames(df_ncc) %in% gg)],each=2),1:2,sep=":")
# #rownames(fits) <- c("HR","Lower CI","Upper CI","p.value")
# fits
# 
# 
# 
# 
# ### Predictive from Joe 
# 
# 
# # PIK3CA 
# 
# df_ncc = full_df[-which(full_df$fmiAssayed==0 & full_df$IDFS_EVNT==1),]
# df_ncc$assayed = as.numeric(df_ncc$sequencing & df_ncc$fmiAssayed)
# 
# fits <- sapply(which(colnames(df_ncc) %in% "pik3ca_mut"),function(j){
#   df_ncc$val <- df_ncc[,j]
#   ipw_fit <- wpl(Surv(IDFS, event=IDFS_EVNT) ~ARM +
#       hr_status + 
#       node_status + 
#       age + 
#       chemo_tx + ARM*val, data = df_ncc, 
#     m = 3, samplestat=df_ncc$assayed+df_ncc$IDFS_EVNT, weight.method = "KM",
#     match.var = cbind(df_ncc$ARM,df_ncc$hr_status,df_ncc$node_status),match.int=rep(0,6))#,
#   class(ipw_fit) = "coxph"
#   ipw_fit$score = Inf
#   
#   
#   coef = 7
#   
#   CI = 0.95
#   z <- qnorm((1 + CI)/2, 0, 1)
#   vc = vcov(ipw_fit)[1,coef]
#   old = diag((vcov(ipw_fit)))[c(1,coef)]
#   old_var = sum((old))
#   se_new = sqrt(old_var + 2*vc)
#   
#   est = exp(coefficients(ipw_fit))
#   low = exp(coefficients(ipw_fit) - z*sqrt(diag(vcov(ipw_fit))) )
#   high = exp(coefficients(ipw_fit) + z*sqrt(diag(vcov(ipw_fit))) )
#   
#   cf = coefficients(ipw_fit)[c(1,coef)] %>% sum
#   bm_hr = exp(cf) %>% round(digits=4)
#   bm_hr_low =  exp(cf - z*se_new ) %>% round(digits=4)
#   bm_hr_high = exp(cf + z*se_new ) %>% round(digits=4)
#   vals = rbind('WT'=c(est[c(1)],low[c(1)],high[c(1)]),'BM+'=c(bm_hr,bm_hr_low,bm_hr_high))
#   colnames(vals) = c("HR","Lower","Upper")
#   as.matrix(cbind(vals,p.value = tidy(ipw_fit)[7,6])) %>% t
# })
# fits <- t(fits) %>% matrix(nrow=2,ncol=4,byrow=TRUE)
# resPik3caInt <- fits
# rownames(resPik3caInt) <- c("Wt PIK3CA","Mut PIK3CA")
# colnames(resPik3caInt) <- c("HR","Lower CI","Upper CI","p.value")

























