




data <- data_goya

tte_var <- "PFS"
censor_var <- "PFS_Censor"


covariate_vars <- c("Treatment_Arm", "IPI_Caterories", "Cell_Of_Origin", "FCGR3A")

return_vars = NULL
return_vars = "Treatment_Arm"







data <- data_goya

data$IPI_Caterories[data$IPI_Caterories %in% c("High")] <- NA

table(data$IPI_Caterories)




# data$IPI_Caterories <- factor(data$IPI_Caterories, levels = c("High", "Low-Intermediate", "High-Intermediate"))



covariate_vars <- c("IPI_Caterories")

return_vars = NULL





variable_names = NULL
caption = NULL
font_size = 11
print_pvalues = TRUE
print_adjpvalues = TRUE




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
wrapper_core_cox_regression <- function(data, tte_var, censor_var, covariate_vars, return_vars = NULL, variable_names = NULL, caption = NULL, font_size = 11, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
  
  
  ### Parce the regression summary
  
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
  
  
  ## Append results from regression
  
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
  
  
  
  ### Return results for the covariates of interest defined by `return_vars`
  
  if(is.null(return_vars)){
    
    ## Return results for all covariates
    res <- coef_info[coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
    
  }else{
    
    ## Return results for covariates defined by `return_vars`
    stopifnot(any(return_vars %in% coef_info$covariate))
    
    res <- coef_info[coef_info$covariate %in% return_vars & coef_info$coefficient %in% rownames(regression_summ$coefficients), , drop = FALSE]
    
    ## If for a factor covariate that should be returned the reference level has zero count, results are set to NA becasue this levels is not used as a reference which means that it is not possible to fit a model that we want.
    if(any(res$reference_indx > 1)){
      
      res[, -which(colnames(res) %in% c("covariate", "levels", "reference", "reference_indx", "coefficient"))] <- NA
      
    }
    
    
  }
  
  
  ### Prepare the output data frame that will be dispalyed. All columns in `out` are characters.
  
  
  
  
  
  
}



















