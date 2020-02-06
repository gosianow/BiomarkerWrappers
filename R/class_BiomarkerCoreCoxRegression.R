#' @include class_BiomarkerClass.R
NULL

###############################################################################
### BiomarkerCoreCoxRegression class
###############################################################################


#' @rdname BiomarkerClass-class
BiomarkerCoreCoxRegression <- setClass("BiomarkerCoreCoxRegression", 
  contains = "BiomarkerClass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BiomarkerCoreCoxRegression", function(object){
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


#' @rdname BiomarkerClass-class
#' @export
setMethod("BiomarkerKable", "BiomarkerCoreCoxRegression", function(x, caption = NULL, font_size = 11, block_vars = NULL){
  
  
  res <- BiomarkerResults(x)
  out <- BiomarkerOutput(x)
  
  if(is.null(caption)){
    caption <- BiomarkerCaption(x)
  }
  
  
  kable <- knitr::kable(out, caption = caption, booktabs = TRUE, linesep = "", row.names = FALSE) %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "bordered"), latex_options = c("HOLD_position"), full_width = FALSE, font_size = font_size) %>% 
    kableExtra::column_spec(which(colnames(out) %in% c("HR", "P-value", "Adj. P-value")), bold = TRUE) 
  
  
  
  # --------------------------------------------------------------------------
  ## Define rows that should be colored
  # --------------------------------------------------------------------------
  
  
  ## Strata variables and the covariate description are appended at the beginning of the 'res' data.frame, before the 'levels' column.
  
  which_columns <- colnames(res)[1:(which(colnames(res) == "levels") - 1)]
  
  if(is.null(block_vars)){
    ### By default color per covariate block
    block_vars <- which_columns[length(which_columns)]
    
  }
  
  stopifnot(all(block_vars %in% which_columns))
  
  ## Identify every second changing character group
  
  
  which_row_spec <- indicate_blocks(d = res, block_vars = block_vars, return = "seq")
  
  if(!is.null(which_row_spec)){
    kable <- kable %>% 
      kableExtra::row_spec(which_row_spec, background = "#F2F4F4")
  }
  
  
  return(kable)
  
})



setMethod("show", "BiomarkerCoreCoxRegression", function(object){
  
  print(BiomarkerKable(object))
  
})

























