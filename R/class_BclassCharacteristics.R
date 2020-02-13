#' @include class_Bclass.R 
NULL

###############################################################################
### Set class unions
###############################################################################

### To allow header = NULL
setClassUnion("integerNULL", c("integer", "NULL"))


###############################################################################
### BclassCharacteristics class
###############################################################################


#' @rdname BclassCharacteristics-class
BclassCharacteristics <- setClass("BclassCharacteristics", 
  slots = c(header = "integerNULL"),
  contains = "Bclass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BclassCharacteristics", function(object){
  # It has to return TRUE when valid object!
  
  ### TODO There can be more checks about the exact covariates in the data frames.
  
  if(nrow(object@results) == nrow(object@output)){
    out <- TRUE
  }else{
    return(paste0("Different number of rows for 'results' and 'output'!"))
  }
  
  if(!is.null(object@header)){
    if(sum(object@header) == ncol(object@output)){
      out <- TRUE
    }else{
      return(paste0("Sum of values in 'header' must be equal to the number of columns in 'output'!"))
    }
  }
  
  
  return(out)
  
})


################################################################################
### Accessor methods
################################################################################


#' @rdname BclassCharacteristics-class
#' @export
setGeneric("Bheader", function(x, ...) standardGeneric("Bheader"))


#' @rdname BclassCharacteristics-class
#' @export
setMethod("Bheader", "BclassCharacteristics", function(x) x@header )


#' @rdname BclassCharacteristics-class
#' @export
setMethod("Bheader", "NULL", function(x) NULL )


#' @rdname BclassCharacteristics-class
#' @export
setGeneric("Bheader<-", function(x, value) standardGeneric("Bheader<-"))


#' @rdname BclassCharacteristics-class
#' @export
setMethod("Bheader<-", "BclassCharacteristics", function(x, value){
  
  BclassRegression(results = Bresults(x), output = Boutput(x), caption = Bcaption(x), header = value)
  
})






################################################################################
### Show method
################################################################################


#' @rdname BclassCharacteristics-class
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
      kableExtra::row_spec(which_row_spec_isna, background = "#fafcff")
  }
  
  
  return(kable)
  
})





setMethod("show", "BclassCharacteristics", function(object){
  
  print(Bkable(object))
  
})

























