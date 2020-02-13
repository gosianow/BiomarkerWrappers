



# data <- data_goya
# tte_var <- "PFS"
# censor_var <- "PFS_Censor"
# covariate_vars <- c("Treatment_Arm", "IPI_Caterories", "Cell_Of_Origin", "FCGR3A")
# return_vars = NULL
# variable_names = NULL
# caption = NULL
# print_nevent = FALSE
# print_pvalues = TRUE
# print_adjpvalues = TRUE


#' Cox regression with simple additive model
#' 
#' Cox regression with simple additive model.
#' 
#' @param data Data frame.
#' @param tte_var Name of the time-to-event variable. This variable must be numeric.
#' @param censor_var Name of the censor variable. It has to be numeric and encode 1 for event and 0 for censor.
#' @param covariate_vars Vector with names of covariate that should be included in the formula.
#' @param return_vars Vector with names of covariate that for which the statistics should be returned. If NULL, sattistics for all covariates are returned.
#' @details 
#' If for a factor covariate that should be returned the reference level has zero count, results are set to NA becasue this levels is not used as a reference which means that it is not possible to fit a model that we want.
wrapper_core_cox_regression_simple <- function(data, tte_var, censor_var, covariate_vars, return_vars = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
  
  
  ### Keep non-missing data
  
  data <- data[complete.cases(data[, c(tte_var, censor_var, covariate_vars)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  if(print_adjpvalues){
    print_pvalues <- TRUE
  }
  
  # --------------------------------------------------------------------------
  # Cox regression
  # --------------------------------------------------------------------------
  
  ## Create the formula
  formula_covariates <- paste0(covariate_vars, collapse = " + ")
  f <- as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_covariates))
  
  
  ## Fit the Cox model
  regression_fit <- survival::coxph(f, data)
  regression_summ <- summary(regression_fit)
  
  
  # mm <- model.matrix(as.formula(paste0(" ~ ", formula_covariates)), data)
  # h(mm)
  
  # --------------------------------------------------------------------------
  ### Parce the regression summary
  # --------------------------------------------------------------------------
  
  ## Generate data frame with coefficient names and levels and information about reference groups
  
  coef_info <- lapply(1:length(covariate_vars), function(i){
    # i = 1
    
    if(covariate_class[i] %in% c("numeric", "integer")){
      
      out <- data.frame(covariate = covariate_vars[i], levels = "", reference = "", reference_indx = 0, n_levels = 0, n_reference = 0, nevent_levels = 0, nevent_reference = 0,
        stringsAsFactors = FALSE)
      
      return(out)
      
    }else{
      
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
      
      out <- data.frame(covariate = covariate_vars[i], levels = levels(data[, covariate_vars[i]]), reference = names(reference_indx), reference_indx = as.numeric(reference_indx), n_levels = as.numeric(tbl), n_reference = as.numeric(tbl[reference_indx]), nevent_levels = as.numeric(tbl_event), nevent_reference = as.numeric(tbl_event[reference_indx]),
        stringsAsFactors = FALSE)
      
      return(out) 
      
    }
    
  })
  
  coef_info <- plyr::rbind.fill(coef_info)
  
  coef_info$coefficient <- paste0(coef_info$covariate, coef_info$levels)
  
  rownames(coef_info) <- coef_info$coefficient
  
  # --------------------------------------------------------------------------
  ## Append results from regression
  # --------------------------------------------------------------------------
  
  conf.int <- data.frame(coefficient = rownames(regression_summ$conf.int), regression_summ$conf.int[, c("exp(coef)", "lower .95", "upper .95")], stringsAsFactors = FALSE)
  colnames(conf.int) <- c("coefficient", "HR", "CI95_lower", "CI95_upper")
  
  coefficients <- data.frame(coefficient = rownames(regression_summ$coefficients), regression_summ$coefficients[, c("Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE)
  colnames(coefficients) <- c("coefficient", "pvalue")
  
  coef_info$n <- regression_summ$n
  coef_info$nevent <- regression_summ$nevent
  
  coef_info <- coef_info %>% 
    left_join(conf.int, by = "coefficient") %>% 
    left_join(coefficients, by = "coefficient")
  
  
  ## Calculate adjusted p-values using the Benjamini & Hochberg method
  coef_info$adj_pvalue <- p.adjust(coef_info$pvalue, method = "BH")
  
  
  # --------------------------------------------------------------------------
  ### Return results for the covariates of interest defined by `return_vars`
  # --------------------------------------------------------------------------
  
  if(is.null(return_vars)){
    
    ## Return results for all covariates
    res <- coef_info[coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
    
  }else{
    
    ## Return results for covariates defined by `return_vars`
    stopifnot(any(return_vars %in% coef_info$covariate))
    
    res <- coef_info[coef_info$covariate %in% return_vars & coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
    
    ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA becasue eventually this level is not used as a reference in the fitted model.
    if(any(res$reference_indx > 1)){
      
      res[, -which(colnames(res) %in% c("covariate", "levels", "reference", "reference_indx", "coefficient"))] <- NA
      
    }
    
    
  }
  
  # --------------------------------------------------------------------------
  ### Prepare the output data frame that will be dispalyed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Effect = format_vs(res$levels, res$reference),
    `Total n` = as.character(res$n),
    `Subgroup n` = format_vs(res$n_levels, res$n_reference),
    `Total events` = as.character(res$nevent),
    `Subgroup events` = format_vs(res$nevent_levels, res$nevent_reference),
    `HR` = as.character(round(res$HR, 2)),
    `95% CI` = format_CIs(res$CI95_lower, res$CI95_upper),
    `P-value` = format_pvalues(res$pvalue),
    `Adj. P-value` = format_pvalues(res$adj_pvalue),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  stopifnot(all(sapply(out, class) == "character"))
  
  if(!print_nevent){
    out$`Total events` <- NULL
    out$`Subgroup events` <- NULL
  }else{
    if(all(out$`Subgroup events` == "")){
      out$`Subgroup events` <- NULL
    }
  }
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  if(!print_adjpvalues){
    out$`Adj. P-value` <- NULL
  }
  
  ### For all numerical covariate Effect and Subgroup n are empty and we do not display them.
  if(all(out$Effect == "")){
    out$Effect <- NULL
  }
  
  if(all(out$`Subgroup n` == "")){
    out$`Subgroup n` <- NULL
  }
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Covariate effect on ", variable_names[tte_var], ". ", 
      "Cox regression model includes ", paste0(variable_names[covariate_vars], collapse = ", "), ".")
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  
  
  bout <- BclassRegression(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}









# data <- data_goya
# tte_var <- "PFS"
# censor_var <- "PFS_Censor"
# covariate_vars <- c("IPI_Caterories", "FCGR3A")
# return_vars = NULL
# 
# strat1_var = "Treatment_Arm"
# strat2_var = NULL
# 
# variable_names = NULL
# caption = NULL
# print_nevent = FALSE
# print_pvalues = TRUE
# print_adjpvalues = TRUE



#' @inheritParams wrapper_core_cox_regression_simple
#' 
#' Cox regression with simple additive model within subgroups defined by strata.
#' 
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
wrapper_core_cox_regression_simple_strat <- function(data, tte_var, censor_var, covariate_vars, return_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
  
  data <- data[complete.cases(data[, c(tte_var, censor_var, covariate_vars, strat1_var, strat2_var)]), ]
  
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
      
      
      wrapper_res <- wrapper_core_cox_regression_simple(data = data_strata1, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, return_vars = return_vars, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      res <- Bresults(wrapper_res)
      out <- Boutput(wrapper_res)
      
      
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      Bresults(wrapper_res) <- res
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      Boutput(wrapper_res) <- out
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
    
    wrapper_res <- BclassRegression(results = res, output = out, caption = Bcaption(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
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
  
  
  wrapper_res <- BclassRegression(results = res, output = out, caption = Bcaption(wrapper_res[[1]]))
  
  return(wrapper_res)
  
  
}











#' @inheritParams wrapper_core_cox_regression_simple_strat
#' 
#' Cox regression estimating biomarker effect 
#' 
#' @param biomarker_vars Vector of biomaker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
wrapper_cox_regression_biomarker <- function(data, tte_var, censor_var, biomarker_vars, adjustment_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include biomarker variables
  stopifnot(length(intersect(biomarker_vars, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(biomarker_vars, adjustment_vars), c(strat1_var, strat2_var))) == 0)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 2
    
    covariate_vars <- c(biomarker_vars[i], adjustment_vars)
    return_vars <- biomarker_vars[i]
    
    
    wrapper_res <- wrapper_core_cox_regression_simple_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, return_vars = return_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
    return(wrapper_res)
    
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
  }
  
  
  ### Replace NAs with "" for columns that are missing for numerical biomarkers
  
  out$Effect[is.na(out$Effect)] <- ""
  out$`Subgroup n`[is.na(out$`Subgroup n`)] <- ""
  if("Subgroup events" %in% colnames(out)){
    out$`Subgroup events`[is.na(out$`Subgroup events`)] <- ""
  }
  
  
  
  ### Rename 'Covariate' column name to 'Biomarker'
  
  colnames(res)[colnames(res) == "covariate"] <- "biomarker"
  colnames(out)[colnames(out) == "Covariate"] <- "Biomarker"
  
  
  ### Generate caption
  
  
  if(is.null(caption)){
    
    caption <- paste0("Biomarker effect on ", variable_names[tte_var], ". ")
    
    
    if(is.null(adjustment_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only the biomarker.")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes the biomarker and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  wrapper_res <- BclassRegression(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}








# data <- data_goya
# tte_var <- "PFS"
# censor_var <- "PFS_Censor"
# 
# treatment_var = "Treatment_Arm"
# biomarker_vars <- c("FCGR2B_cat2", "FCGR3A_cat2")
# adjustment_vars <- "IPI_Caterories"
# 
# strat2_var = "Cell_Of_Origin"
# 
# variable_names = NULL
# caption = NULL
# print_nevent = FALSE
# print_pvalues = TRUE
# print_adjpvalues = TRUE




#' @inheritParams wrapper_core_cox_regression_simple_strat
#' 
#' Cox regression estimating treatment effect within biomaker subgroups
#' 
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector of biomaker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
wrapper_cox_regression_treatment <- function(data, tte_var, censor_var, treatment_var, biomarker_vars, adjustment_vars = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  ### TODO Allow biomarker_vars = NULL.
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  ### Adjustment variables cannot include treatment
  stopifnot(length(intersect(treatment_var, adjustment_vars)) == 0)
  
  ### Model variables cannot include strata variables
  stopifnot(length(intersect(c(treatment_var, adjustment_vars), c(biomarker_vars, strat2_var))) == 0)
  
  
  vars_class <- sapply(data[, c(treatment_var, biomarker_vars)], class)
  stopifnot(all(vars_class == "factor"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  
  wrapper_res <- lapply(1:length(biomarker_vars), function(i){
    # i = 1
    
    covariate_vars <- c(treatment_var, adjustment_vars)
    return_vars <- treatment_var
    strat1_var <- biomarker_vars[i]
    
    wrapper_res <- wrapper_core_cox_regression_simple_strat(data = data, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, return_vars = return_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    
    
    ### Rename strata column name to 'Biomarker Subgroup'
    
    res <- Bresults(wrapper_res)
    out <- Boutput(wrapper_res)
    
    colnames(res)[colnames(res) == biomarker_vars[i]] <- "biomarker_subgroup"
    colnames(out)[colnames(out) == variable_names[biomarker_vars[i]]] <- "Biomarker Subgroup"
    
    ### Treatment is the same for all the biomarkers and biomarker info is missing. Thus, we replace covariate with biomarker name. Swap Biomarker Subgroup with Biomarker.
    
    res[, colnames(res) == "covariate"] <- biomarker_vars[i]
    out[, colnames(out) == "Covariate"] <- variable_names[biomarker_vars[i]]
    
    colnames(res)[colnames(res) == "covariate"] <- "biomarker"
    colnames(out)[colnames(out) == "Covariate"] <- "Biomarker"
    
    res <- dplyr::select(res, c(strat2_var, "biomarker", "biomarker_subgroup"), everything())
    out <- dplyr::select(out, c(variable_names[strat2_var], "Biomarker", "Biomarker Subgroup"), everything())
    
    
    Bresults(wrapper_res) <- res
    Boutput(wrapper_res) <- out
    
    
    return(wrapper_res)
    
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
  }
  
  
  ### Change the name to Treatment Effect
  
  colnames(out)[colnames(out) == "Effect"] <- "Treatment Effect"
  
  
  ### Generate caption
  
  
  if(is.null(caption)){
    
    caption <- paste0("Treatment effect on ", variable_names[tte_var], ". ")
    
    
    if(is.null(adjustment_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only the treatment.")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes the treatment and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  wrapper_res <- BclassRegression(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}



###############################################################################
# Model with interaction
###############################################################################



# data <- data_goya
# tte_var <- "PFS"
# censor_var <- "PFS_Censor"
# 
# interaction1_var <- "FCGR3A_cat2"
# interaction2_var <- "Treatment_Arm"
# covariate_vars <- c("IPI_Caterories", "Cell_Of_Origin")
# 
# 
# variable_names = NULL
# caption = NULL
# print_nevent = FALSE
# print_pvalues = TRUE
# print_adjpvalues = TRUE



#' @inheritParams wrapper_core_cox_regression_simple
#' 
#' Cox regression with additive model with interaction
#' 
#' @param interaction1_var Data frame.
#' @param interaction2_var Name of the time-to-event variable. This variable must be numeric.
wrapper_core_cox_regression_interaction <- function(data, tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
  
  
  ### Keep non-missing data
  
  data <- data[complete.cases(data[, c(tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  if(print_adjpvalues){
    print_pvalues <- TRUE
  }
  
  
  # --------------------------------------------------------------------------
  # Cox regression
  # --------------------------------------------------------------------------
  
  
  ## Create the formula
  formula_covariates <- paste0(paste0(covariate_vars, collapse = " + "), " + ", interaction1_var, " * ", interaction2_var)
  f <- as.formula(paste0("Surv(", tte_var, ", ", censor_var, ") ~ ", formula_covariates))
  
  
  ## Fit the Cox model
  regression_fit <- survival::coxph(f, data)
  regression_summ <- summary(regression_fit)
  
  
  # mm <- model.matrix(as.formula(paste0(" ~ ", formula_covariates)), data)
  # h(mm)
  
  
  # --------------------------------------------------------------------------
  ### Parce the regression summary for the interaction terms
  # --------------------------------------------------------------------------
  
  ## Generate data frame with coefficient names and levels and information about reference groups
  
  if(class(interaction1_var) %in% c("numeric", "integer") && class(interaction2_var) %in% c("numeric", "integer")){
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = "", reference1 = "", reference1_indx = 0, covariate2 = interaction2_var, levels2 = "", reference2 = "", reference2_indx = 0, 
      stringsAsFactors = FALSE)
    
  }else if(class(interaction1_var) %in% c("numeric", "integer") && class(interaction2_var) %in% c("factor")){
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction2_var])
    
    if(tbl[1] > 0){
      reference_indx <- 1
      names(reference_indx) <- names(tbl[1])
    }else{
      reference_indx <- rev(which(tbl > 0))[1]
    }
    
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = "", reference1 = "", reference1_indx = 0, covariate2 = interaction2_var, levels2 = levels(data[, interaction2_var]), reference2 = names(reference_indx), reference2_indx = as.numeric(reference_indx),
      stringsAsFactors = FALSE)
    
  }else if(class(interaction1_var) %in% c("factor") && class(interaction2_var) %in% c("numeric", "integer")){
    
    ## Check if the first level has non zero counts so it can be used as a reference group in the regression
    ## Actually, when the first level has zero counts, then the last level with non-zero counts is used as a reference 
    tbl <- table(data[, interaction1_var])
    
    if(tbl[1] > 0){
      reference_indx <- 1
      names(reference_indx) <- names(tbl[1])
    }else{
      reference_indx <- rev(which(tbl > 0))[1]
    }
    
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = levels(data[, interaction1_var]), reference1 = names(reference_indx), reference1_indx = as.numeric(reference_indx), covariate2 = interaction2_var, levels2 = "", reference2 = "", reference2_indx = 0,
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
    
    
    out <- data.frame(covariate1 = interaction1_var, levels1 = rep(levels(data[, interaction1_var]), times = nlevels(data[, interaction2_var])), reference1 = names(reference1_indx), reference1_indx = as.numeric(reference1_indx), covariate2 = interaction2_var, levels2 = rep(levels(data[, interaction2_var]), each = nlevels(data[, interaction1_var])), reference2 = names(reference2_indx), reference2_indx = as.numeric(reference2_indx),
      stringsAsFactors = FALSE)
    
    
  }
  
  
  coef_info <- out
  
  coef_info$coefficient <- paste0(coef_info$covariate1, coef_info$levels1, ":", coef_info$covariate2, coef_info$levels2)
  
  rownames(coef_info) <- coef_info$coefficient
  
  # --------------------------------------------------------------------------
  ## Append results from regression
  # --------------------------------------------------------------------------
  
  conf.int <- data.frame(coefficient = rownames(regression_summ$conf.int), regression_summ$conf.int[, c("exp(coef)", "lower .95", "upper .95")], stringsAsFactors = FALSE)
  colnames(conf.int) <- c("coefficient", "HR", "CI95_lower", "CI95_upper")
  
  coefficients <- data.frame(coefficient = rownames(regression_summ$coefficients), regression_summ$coefficients[, c("Pr(>|z|)"), drop = FALSE], stringsAsFactors = FALSE)
  colnames(coefficients) <- c("coefficient", "pvalue")
  
  coef_info$n <- regression_summ$n
  coef_info$nevent <- regression_summ$nevent
  
  coef_info <- coef_info %>% 
    left_join(conf.int, by = "coefficient") %>% 
    left_join(coefficients, by = "coefficient")
  
  
  ## Calculate adjusted p-values using the Benjamini & Hochberg method
  coef_info$adj_pvalue <- p.adjust(coef_info$pvalue, method = "BH")
  
  
  # --------------------------------------------------------------------------
  ### Return results 
  # --------------------------------------------------------------------------
  
  
  res <- coef_info[coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
  
  ## If for a factor covariate that should be returned the (first) reference level has zero count, results are set to NA becasue eventually this level is not used as a reference in the fitted model.
  if(any(res$reference_indx > 1)){
    
    res[, -which(colnames(res) %in% c("covariate1", "levels1", "reference1", "reference1_indx", "covariate2", "levels2", "reference2", "reference2_indx", "coefficient"))] <- NA
    
  }
  
  
  
  # --------------------------------------------------------------------------
  ### Prepare the output data frame that will be dispalyed. All columns in `out` are characters.
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate1 = variable_names[res$covariate1], 
    Effect1 = format_vs(res$levels1, res$reference1),
    Covariate2 = variable_names[res$covariate2], 
    Effect2 = format_vs(res$levels2, res$reference2),
    `Total n` = as.character(res$n),
    `Total events` = as.character(res$nevent),
    `HR` = as.character(round(res$HR, 2)),
    `95% CI` = format_CIs(res$CI95_lower, res$CI95_upper),
    `P-value` = format_pvalues(res$pvalue),
    `Adj. P-value` = format_pvalues(res$adj_pvalue),
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_nevent){
    out$`Total events` <- NULL
  }
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  if(!print_adjpvalues){
    out$`Adj. P-value` <- NULL
  }
  
  ### For all numerical covariate Effect and Subgroup n are empty and we do not display them.
  if(all(out$Effect == "")){
    out$Effect <- NULL
  }
  
  if(all(out$`Subgroup n` == "")){
    out$`Subgroup n` <- NULL
  }
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Effect of interaction between ", variable_names[interaction1_var], " and ", variable_names[interaction2_var], " on ", variable_names[tte_var], ". ")
    
    if(!is.null(covariate_vars)){
      caption <- paste0(caption, paste0("Cox regression model additionally includes ", paste0(variable_names[covariate_vars], collapse = ", "), ". "))
    }
    
  }
  
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  
  
  bout <- BclassRegression(results = res, output = out, caption = caption)
  
  
  return(bout)
  
}




#' @inheritParams wrapper_core_cox_regression_interaction
#' 
#' Cox regression with additive model with interaction within subgroups defined by strata
#' 
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
wrapper_core_cox_regression_interaction_strat <- function(data, tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
  
  data <- data[complete.cases(data[, c(tte_var, censor_var, interaction1_var, interaction2_var, covariate_vars, strat1_var, strat2_var)]), ]
  
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
      
      
      wrapper_res <- wrapper_core_cox_regression_interaction(data = data_strata1, tte_var = tte_var, censor_var = censor_var, interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      
      res <- Bresults(wrapper_res)
      out <- Boutput(wrapper_res)
      
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      Bresults(wrapper_res) <- res
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      Boutput(wrapper_res) <- out
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
    
    wrapper_res <- BclassRegression(results = res, output = out, caption = Bcaption(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
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
  
  
  wrapper_res <- BclassRegression(results = res, output = out, caption = Bcaption(wrapper_res[[1]]))
  
  return(wrapper_res)
  
  
}



#' @inheritParams wrapper_core_cox_regression_interaction_strat
#' 
#' Cox regression estimating interaction effect between biomaker and treatment
#' 
#' @param treatment_var Name of column with treatment information.
#' @param biomarker_vars Vector of biomaker names.
#' @param adjustment_vars Vector of covariate names used for adjustment.
wrapper_cox_regression_interaction <- function(data, tte_var, censor_var, treatment_var, biomarker_vars, adjustment_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
    
    wrapper_res <- wrapper_core_cox_regression_interaction_strat(data = data, tte_var = tte_var, censor_var = censor_var, interaction1_var = interaction1_var, interaction2_var = interaction2_var, covariate_vars = covariate_vars, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, Bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, Boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
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
    
    caption <- paste0("Effect of interaction between ", "biomarker", " and ", "treatment", " on ", variable_names[tte_var], ". ")
    
    
    if(is.null(adjustment_vars)){
      
      caption <- paste0(caption, "Unadjusted, unstratified analysis. Cox regression model includes only the biomarker and treatment.")
      
    }else{
      
      caption <- paste0(caption, "Adjusted, unstratified analysis. Cox regression model includes the biomarker, treatment and ", paste0(variable_names[adjustment_vars], collapse = ", "), ". ")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  wrapper_res <- BclassRegression(results = res, output = out, caption = caption)
  
  return(wrapper_res)
  
  
}









