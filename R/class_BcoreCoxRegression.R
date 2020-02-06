#' @include class_Bclass.R
NULL

###############################################################################
### BcoreCoxRegression class
###############################################################################


#' @rdname Bclass-class
BcoreCoxRegression <- setClass("BcoreCoxRegression", 
  contains = "Bclass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BcoreCoxRegression", function(object){
  # It has to return TRUE when valid object!
  
  ### There can be more checks about the exact covariates in the data frames.
  
  if(nrow(object@results) == nrow(object@output)){
    out <- TRUE
  }else{
    return(paste0("Different number of rows for 'results' and 'output'!"))
  }
  
  return(out)
  
})



################################################################################
### Show method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("Bkable", "BcoreCoxRegression", function(x, caption = NULL, font_size = 11, block_vars = NULL){
  
  
  res <- Bresults(x)
  out <- Boutput(x)
  
  if(is.null(caption)){
    caption <- Bcaption(x)
  }
  
  
  kable <- knitr::kable(out, caption = caption, booktabs = TRUE, linesep = "", row.names = FALSE) %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "bordered"), latex_options = c("HOLD_position"), full_width = FALSE, font_size = font_size) %>% 
    kableExtra::column_spec(which(colnames(out) %in% c("HR", "P-value", "Adj. P-value")), bold = TRUE) 
  
  
  
  # --------------------------------------------------------------------------
  ## Define rows that should be colored
  # --------------------------------------------------------------------------
  
  ### By default color per covariate/biomarker block
  if(is.null(block_vars)){
    
    block_vars <- colnames(res)[which(colnames(res) %in% c("covariate", "biomarker"))]
    
  }
  
  
  ## Identify every second changing character group
  which_row_spec <- indicate_blocks(d = res, block_vars = block_vars, return = "seq")
  
  
  if(!is.null(which_row_spec)){
    kable <- kable %>% 
      kableExtra::row_spec(which_row_spec, background = "#F2F4F4")
  }
  
  
  return(kable)
  
})



setMethod("show", "BcoreCoxRegression", function(object){
  
  print(Bkable(object))
  
})

























