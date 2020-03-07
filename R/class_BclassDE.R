#' @include class_Bclass.R
NULL

###############################################################################
### BclassDE class
###############################################################################


#' @rdname Bclass-class
BclassDE <- setClass("BclassDE", 
  contains = "Bclass")


# --------------------------------------------------------------------------
# Object validation
# --------------------------------------------------------------------------


setValidity("BclassDE", function(object){
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
setMethod("bresults<-", "BclassDE", function(x, value){
  
  BclassDE(results = value, output = boutput(x), caption = bcaption(x), header = bheader(x))
  
})

### output


#' @rdname Bclass-class
#' @export
setMethod("boutput<-", "BclassDE", function(x, value){
  
  BclassDE(results = bresults(x), output = value, caption = bcaption(x), header = bheader(x))
  
})


### caption

#' @rdname Bclass-class
#' @export
setMethod("bcaption<-", "BclassDE", function(x, value){
  
  BclassDE(results = bresults(x), output = boutput(x), caption = value, header = bheader(x))
  
})


### header

#' @rdname Bclass-class
#' @export
setMethod("bheader<-", "BclassDE", function(x, value){
  
  BclassDE(results = bresults(x), output = boutput(x), caption = bcaption(x), header = value)
  
})



################################################################################
### Show method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("bkable", "BclassDE", function(x, caption = NULL, header = NULL, font_size = NULL, full_width = NULL){
  
  
  res <- bresults(x)
  out <- boutput(x)
  
  if(nrow(res) == 0){
    message <- paste0("\n", bcaption(x), "\n\n")
    return(message)
  }
  
  
  if(is.null(caption)){
    caption <- bcaption(x)
  }
  
  if(is.null(header)){
    header <- bheader(x)
  }
  
  if(is.null(font_size)){
    font_size <- getOption("bkable_font_size", default = 11)
  }
  
  if(is.null(full_width)){
    full_width <- getOption("bkable_full_width", default = TRUE)
  }
  
  
  kable <- knitr::kable(out, caption = caption, booktabs = TRUE, linesep = "", row.names = FALSE) %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "bordered", "striped"), latex_options = c("HOLD_position", "striped"), full_width = full_width, font_size = font_size)
  
  
  if(!is.null(header)){
    kable <- kable %>% 
      add_header_above(header)
  }
  
  
  return(kable)
  
})



setMethod("show", "BclassDE", function(object){
  
  x <- bkable(object)
  
  if(is.character(x)){
    cat(x)
  }else{
    print(x)
  }
  
})






















