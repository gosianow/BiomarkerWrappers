#' @include class_Bclass.R
NULL

###############################################################################
### BclassTesting class
###############################################################################


#' @rdname Bclass-class
BclassTesting <- setClass("BclassTesting", 
  contains = "Bclass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BclassTesting", function(object){
  # It has to return TRUE when valid object!
  
  ### TODO There can be more checks about the exact covariates in the data frames.
  
  out <- TRUE
  
  return(out)
  
})

################################################################################
### Replacing methods
################################################################################

### results

#' @rdname Bclass-class
#' @export
setMethod("Bresults<-", "BclassTesting", function(x, value){
  
  BclassTesting(results = value, output = Boutput(x), caption = Bcaption(x), header = Bheader(x))
  
})

### output


#' @rdname Bclass-class
#' @export
setMethod("Boutput<-", "BclassTesting", function(x, value){
  
  BclassTesting(results = Bresults(x), output = value, caption = Bcaption(x), header = Bheader(x))
  
})


### caption

#' @rdname Bclass-class
#' @export
setMethod("Bcaption<-", "BclassTesting", function(x, value){
  
  BclassTesting(results = Bresults(x), output = Boutput(x), caption = value, header = Bheader(x))
  
})


### header

#' @rdname Bclass-class
#' @export
setMethod("Bheader<-", "BclassTesting", function(x, value){
  
  BclassTesting(results = Bresults(x), output = Boutput(x), caption = Bcaption(x), header = value)
  
})



################################################################################
### Show method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("Bkable", "BclassTesting", function(x, caption = NULL, header = NULL, font_size = 11, block_vars = NULL){
  
  
  res <- Bresults(x)
  out <- Boutput(x)
  
  
  if(is.null(caption)){
    caption <- Bcaption(x)
  }
  
  if(is.null(header)){
    header <- Bheader(x)
  }
  
  
  kable <- knitr::kable(out, caption = caption, booktabs = TRUE, linesep = "", row.names = FALSE) %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "bordered"), latex_options = c("HOLD_position"), full_width = FALSE, font_size = font_size) %>% 
    kableExtra::column_spec(which(colnames(out) %in% c("HR", "OR", "FC", "P-value", "Adj. P-value")), bold = TRUE) %>% 
    add_header_above(header)
  
  
  
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



setMethod("show", "BclassTesting", function(object){
  
  print(Bkable(object))
  
})

























