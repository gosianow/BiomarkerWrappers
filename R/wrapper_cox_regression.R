



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


#' Cox regression with additive model
#' 
#' Cox regression with additive model.
#' 
#' @param data Data frame.
#' @param tte_var Name of the time-to-event variable. This variable must be numeric.
#' @param censor_var Name of the censor variable. It has to be numeric and encode 1 for event and 0 for censor.
#' @param covariate_vars Vector with names of covariate that should be included in the formula.
#' @param return_vars Vector with names of covariate that for which the statistics should be returned. If NULL, sattistics for all covariates are returned.
#' @details 
#' If for a factor covariate that should be returned the reference level has zero count, results are set to NA becasue this levels is not used as a reference which means that it is not possible to fit a model that we want.
wrapper_core_cox_regression <- function(data, tte_var, censor_var, covariate_vars, return_vars = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
  
  
  if(!print_nevent){
    out$`Total events` <- NULL
    out$`Subgroup events` <- NULL
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
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
  }
  
  
  bout <- BiomarkerCoreCoxRegression(results = res, output = out, caption = caption)
  
  
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



wrapper_core_cox_regression_strat <- function(data, tte_var, censor_var, covariate_vars, return_vars = NULL, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, print_nevent = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
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
      
      
      wrapper_res <- wrapper_core_cox_regression(data = data_strata1, tte_var = tte_var, censor_var = censor_var, covariate_vars = covariate_vars, return_vars = return_vars, variable_names = variable_names, caption = caption, print_nevent = print_nevent, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      res <- BiomarkerResults(wrapper_res)
      out <- BiomarkerOutput(wrapper_res)
      
      
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      BiomarkerResults(wrapper_res) <- res
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      BiomarkerOutput(wrapper_res) <- out
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, BiomarkerResults))
    out <- plyr::rbind.fill(lapply(wrapper_res, BiomarkerOutput))
    
    wrapper_res <- BiomarkerCoreCoxRegression(results = res, output = out, caption = BiomarkerCaption(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, BiomarkerResults))
  out <- plyr::rbind.fill(lapply(wrapper_res, BiomarkerOutput))
  
  
  ### Remove dummy columns
  
  if(strat2_var == "strat2_dummy"){
    res$strat2_dummy <- NULL
    out$`strat2 dummy` <- NULL
  }
  if(strat1_var == "strat1_dummy"){
    res$strat1_dummy <- NULL
    out$`strat1 dummy` <- NULL
  }
  
  
  wrapper_res <- BiomarkerCoreCoxRegression(results = res, output = out, caption = BiomarkerCaption(wrapper_res[[1]]))
  
  return(wrapper_res)
  
  
}


























