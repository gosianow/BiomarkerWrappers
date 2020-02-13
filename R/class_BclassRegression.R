#' @include class_Bclass.R
NULL

###############################################################################
### BclassRegression class
###############################################################################


#' @rdname Bclass-class
BclassRegression <- setClass("BclassRegression", 
  contains = "Bclass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BclassRegression", function(object){
  # It has to return TRUE when valid object!
  
  ### TODO There can be more checks about the exact covariates in the data frames.
  
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
setMethod("Bkable", "BclassRegression", function(x, caption = NULL, font_size = 11, block_vars = NULL){
  
  
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
    
    block_vars <- colnames(res)[which(colnames(res) %in% c("covariate", "biomarker", "covariate1"))]
    
  }
  
  
  ## Identify every second changing character group
  which_row_spec <- indicate_blocks(d = res, block_vars = block_vars, return = "block")
  
  
  if(!is.null(which_row_spec)){
    kable <- kable %>% 
      kableExtra::row_spec(which_row_spec, background = "#F2F4F4")
  }
  
  
  return(kable)
  
})



setMethod("show", "BclassRegression", function(object){
  
  print(Bkable(object))
  
})

























