#' @include class_Bclass.R 
NULL


###############################################################################
### BclassCharacteristics class
###############################################################################


#' @rdname Bclass-class
BclassCharacteristics <- setClass("BclassCharacteristics", 
  contains = "Bclass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BclassCharacteristics", function(object){
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
setMethod("Bresults<-", "BclassCharacteristics", function(x, value){
  
  BclassCharacteristics(results = value, output = Boutput(x), caption = Bcaption(x), header = Bheader(x))
  
})




### output


#' @rdname Bclass-class
#' @export
setMethod("Boutput<-", "BclassCharacteristics", function(x, value){
  
  BclassCharacteristics(results = Bresults(x), output = value, caption = Bcaption(x), header = Bheader(x))
  
})





### caption



#' @rdname Bclass-class
#' @export
setMethod("Bcaption<-", "BclassCharacteristics", function(x, value){
  
  BclassCharacteristics(results = Bresults(x), output = Boutput(x), caption = value, header = Bheader(x))
  
})


### header



#' @rdname Bclass-class
#' @export
setMethod("Bheader<-", "BclassCharacteristics", function(x, value){
  
  BclassCharacteristics(results = Bresults(x), output = Boutput(x), caption = Bcaption(x), header = value)
  
})






################################################################################
### Show method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("Bkable", "BclassCharacteristics", function(x, caption = NULL, header = NULL, font_size = 11){
  
  
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
    add_header_above(header)
  
  
  
  # --------------------------------------------------------------------------
  ## Define rows that should be colored
  # --------------------------------------------------------------------------
  
  ## Rows with names of covariates. They contain "".
  
  which_row_spec <- which(rowSums(out == "") >= 1)
  
  which_row_spec_isna <- which(out[, 1] %in% c("Total (non-NA)", "NAs"))
  
  
  if(!is.null(which_row_spec)){
    kable <- kable %>% 
      kableExtra::row_spec(which_row_spec, bold = TRUE, background = "#F0F8FF") %>% 
      kableExtra::row_spec(which_row_spec_isna, background = "#fafcff", color = "#666666")
  }
  
  
  return(kable)
  
})





setMethod("show", "BclassCharacteristics", function(object){
  
  print(Bkable(object))
  
})

























