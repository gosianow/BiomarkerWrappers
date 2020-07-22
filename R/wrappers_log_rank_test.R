#' @include wrappers_cox_regression.R
NULL



#' Log-rank test
#' 
#' @param data Data frame.
#' @param tte_var Name of the time-to-event variable. This variable must be numeric.
#' @param censor_var Name of the censor variable. It has to be numeric and encode 1 for event and 0 for censor.
#' @param covariate_var Name of covariate that defines subgroups where the survival is estimated.
#' @param strata_vars Vector with names of covariates that are used as strata.
#' @param variable_names Named vector with variable names. If not supplied, variable names are created by replacing in column names underscores with spaces.
#' @param caption Caption for the table with results.
#' @param print_nevent Logical. Whether to print numbers of events.
#' @param print_mst Logical. Whether to print median survival time (MST).
#' @param print_hr Logical. Whether to print hazard rations estimated with Cox regression.
#' @param print_pvalues Logical. Whether to print p-values.
#' 
#' @examples
#' 
#' data(bdata)
#' 
#' data <- bdata
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' covariate_var <- "Treatment_Arm"
#' 
#' x <- wrapper_core_log_rank_test_simple(data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var)
#' 
#' boutput(x)
#' 
#' 
#' ### Fit a stratified model 
#' 
#' strata_vars <- c("IPI", "Cell_Of_Origin")
#' 
#' x <- wrapper_core_log_rank_test_simple(data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, strata_vars = strata_vars)
#' 
#' boutput(x)
#' 
#' 
#' @export
wrapper_core_log_rank_test_simple <- function(data, tte_var, censor_var, covariate_var, strata_vars = NULL, variable_names = NULL, caption = NULL, print_nevent = TRUE, print_mst = TRUE, print_hr = TRUE, print_pvalues = TRUE){
  
  
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
  
  
  covariate_class <- sapply(data[, covariate_var], class)
  stopifnot(all(covariate_class %in% c("factor")))
  
  
  if(!is.null(strata_vars)){
    strata_class <- sapply(data[, strata_vars], class)
    stopifnot(all(strata_class %in% c("factor")))
  }
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_var, strata_vars)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Prepare data frame with results
  # --------------------------------------------------------------------------
  
  
  
  # -----------------------------------
  # Fit Cox regression to obtain Number of Events, MST and HR
  # -----------------------------------
  
  covariate_vars <- covariate_var
  return_vars <- covariate_var
  
  wrapper_res <- wrapper_core_cox_regression_simple(data = data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, strata_vars = strata_vars, return_vars = return_vars, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_mst = TRUE, print_pvalues = FALSE, print_adjpvalues = FALSE)
  
  
  res <- bresults(wrapper_res)
  out <- boutput(wrapper_res)
  
  
  
  # -----------------------------------
  # Log-rank test
  # -----------------------------------
  
  if(sum(res$n > 0) >= 2){
    ## Create the formula
    
    formula_covariates <- paste0(covariate_var, collapse = " + ")
    
    if(!is.null(strata_vars)){
      formula_strata <- paste0("strata(", paste0(strata_vars, collapse = ", "), ")")
      formula_model <- paste0(formula_covariates, " + ", formula_strata)
    }else{
      formula_model <- formula_covariates
    }
    
    
    f <- stats::as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_model))
    
    
    ## Fit the Cox model
    # regression_fit <- survival::coxph(f, data, ties = "exact")
    # regression_summ <- summary(regression_fit)
    # 
    # pvalue <- regression_summ$sctest["pvalue"]
    
    # mm <- model.matrix(stats::as.formula(paste0(" ~ ", formula_covariates)), data)
    # h(mm)
    
    
    ## Test Survival Curve Differences with Log-rank test
    survdiff_fit <- survival::survdiff(f, data)
    
    pvalue <- pchisq(survdiff_fit$chisq, length(survdiff_fit$n) - 1, lower.tail = FALSE)
    
  }else{
    pvalue <- NA
  }
  
  
  ### Add pvalue
  res$pvalue <- c(pvalue, rep(NA, nrow(res) - 1))
  ### Remove adj_pvalue
  res$adj_pvalue <- NULL
  
  
  # --------------------------------------------------------------------------
  # Prepare the output data frame that will be dispalyed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  out$`P-value` <- format_pvalues(res$pvalue)
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_nevent){
    out$`Total Events` <- NULL
    out$`Events` <- NULL
  }
  
  if(!print_mst){
    out$`MST` <- NULL
    out$`MST 95% CI` <- NULL
  }
  
  if(!print_hr){
    out$`HR` <- NULL
    out$`HR 95% CI` <- NULL
  }
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  

  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0(variable_names[covariate_var], " effect on ", variable_names[tte_var], ". ", 
      ifelse(is.null(strata_vars), "Unstratified ", "Stratified "),
      "log-rank test.")
    
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





#' @rdname wrapper_core_log_rank_test_simple 
#' 
#' @param strat1_var Name of the firts stratification variable.
#' @param strat2_var Name of the second stratification variable.
#' @param print_adjpvalues Logical. Whether to print adjusted p-values.
#' 
#' @examples 
#' 
#' data(bdata)
#' 
#' data <- bdata
#' tte_var <- "PFS"
#' censor_var <- "PFS_Event"
#' covariate_var <- "Treatment_Arm"
#' 
#' strat1_var = "Cell_Of_Origin"
#' 
#' 
#' x <- wrapper_core_log_rank_test_simple_strat(data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, strat1_var = strat1_var)
#' 
#' boutput(x)
#' 
#' @export
wrapper_core_log_rank_test_simple_strat <- function(data, tte_var, censor_var, covariate_var, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = TRUE, print_mst = TRUE, print_hr = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
  stopifnot(length(intersect(c(strat1_var, strat2_var), covariate_var)) == 0)
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(tte_var, censor_var, covariate_var, strat1_var, strat2_var, strata_vars)]), ]
  
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
      
      
      wrapper_res <- wrapper_core_log_rank_test_simple(data = data_strata1, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, strata_vars = strata_vars, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_mst = print_mst, print_hr = print_hr, print_pvalues = print_pvalues)
      
      
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











#' Log-rank test testing biomarker effect 
#' 
#' @inheritParams wrapper_core_log_rank_test_simple_strat
#' @param biomarker_vars Vector of biomaker names.
#' @export
wrapper_log_rank_test_biomarker <- function(data, tte_var, censor_var, biomarker_vars, strata_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = TRUE, print_mst = TRUE, print_hr = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
    
    
    wrapper_res <- wrapper_core_log_rank_test_simple_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_mst = print_mst, print_hr = print_hr, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
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
  
  colnames(res)[colnames(res) == "covariate"] <- "biomarker"
  colnames(out)[colnames(out) == "Covariate"] <- "Biomarker"
  
  
  ### Generate caption
  
  
  if(is.null(caption)){
    
    caption <- paste0("Biomarker effect on ", variable_names[tte_var], ". ")
    
    if(is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified log-rank test.")
      
    }else{
      
      caption <- paste0(caption, "Stratified log-rank test. Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}






#' Log-rank test testinging treatment effect within biomaker subgroups
#' 
#' @inheritParams wrapper_core_log_rank_test_simple_strat
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector with names of categorical biomarkers. When NULL, overall treatment effect is estimated. 
#' @export
wrapper_log_rank_test_treatment <- function(data, tte_var, censor_var, treatment_var, strata_vars = NULL, biomarker_vars = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = TRUE, print_mst = TRUE, print_hr = TRUE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
    
    
    wrapper_res <- wrapper_core_log_rank_test_simple_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_var = covariate_var, strata_vars = strata_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_mst = print_mst, print_hr = print_hr, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
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
    
    if(is.null(strata_vars)){
      
      caption <- paste0(caption, "Unstratified log-rank test.")
      
    }else{
      
      caption <- paste0(caption, "Stratified log-rank test. Stratification factors: ", paste0(variable_names[strata_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}









