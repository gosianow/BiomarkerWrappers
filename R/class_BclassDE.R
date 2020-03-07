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
setMethod("Bresults<-", "BclassDE", function(x, value){
  
  BclassDE(results = value, output = Boutput(x), caption = Bcaption(x), header = Bheader(x))
  
})

### output


#' @rdname Bclass-class
#' @export
setMethod("Boutput<-", "BclassDE", function(x, value){
  
  BclassDE(results = Bresults(x), output = value, caption = Bcaption(x), header = Bheader(x))
  
})


### caption

#' @rdname Bclass-class
#' @export
setMethod("Bcaption<-", "BclassDE", function(x, value){
  
  BclassDE(results = Bresults(x), output = Boutput(x), caption = value, header = Bheader(x))
  
})


### header

#' @rdname Bclass-class
#' @export
setMethod("Bheader<-", "BclassDE", function(x, value){
  
  BclassDE(results = Bresults(x), output = Boutput(x), caption = Bcaption(x), header = value)
  
})



################################################################################
### Show method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("Bkable", "BclassDE", function(x, caption = NULL, header = NULL, font_size = NULL, full_width = NULL){
  
  
  res <- Bresults(x)
  out <- Boutput(x)
  
  if(nrow(res) == 0){
    message <- paste0("\n", Bcaption(x), "\n\n")
    return(message)
  }
  
  
  if(is.null(caption)){
    caption <- Bcaption(x)
  }
  
  if(is.null(header)){
    header <- Bheader(x)
  }
  
  if(is.null(font_size)){
    font_size <- getOption("Bkable_font_size", default = 11)
  }
  
  if(is.null(full_width)){
    full_width <- getOption("Bkable_full_width", default = TRUE)
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
  
  x <- Bkable(object)
  
  if(is.character(x)){
    cat(x)
  }else{
    print(x)
  }
  
})






















