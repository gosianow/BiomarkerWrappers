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
### Show method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("bkable", "BclassCharacteristics", function(x, caption = NULL, header = NULL, font_size = NULL, full_width = NULL){
  
  
  res <- bresults(x)
  out <- boutput(x)
  
  
  if(is.null(caption)){
    caption <- bcaption(x)
  }
  
  if(is.null(header)){
    header <- bheader(x)
  }
  
  ### Use the globally defined font_size
  if(is.null(font_size)){
    font_size <- getOption("bkable_font_size", default = 11)
  }
  
  if(is.null(full_width)){
    full_width <- getOption("bkable_full_width", default = TRUE)
  }
  
  ### We add format as an option because kableExtra::column_spec does not work if format is not specified in the kable call
  format <- getOption("knitr.table.format", default = "html")
  
  ### Make unique labels to have the right numbering of the tables. See https://github.com/rstudio/bookdown/issues/1097
  
  kable <- knitr::kable(out, format = format, caption = caption, booktabs = TRUE, linesep = "", row.names = FALSE, label = gsub("\\.", "", make.names(Sys.time()))) %>%
    kableExtra::kable_styling(bootstrap_options = c("condensed", "bordered"), latex_options = c("HOLD_position"), full_width = full_width, font_size = font_size) 
  
  
  if(!is.null(header)){
    kable <- kable %>% 
      add_header_above(header)
  }
  
  
  
  # --------------------------------------------------------------------------
  ## Define rows that should be colored
  # --------------------------------------------------------------------------
  
  ## Rows with names of covariates. They contain "".
  
  which_row_spec <- which(rowSums(out[, -1, drop = FALSE] == "") >= 1)
  
  which_row_spec_isna <- which(out[, 1] %in% c("N", "NAs"))
  
  
  kable <- kable %>% 
    kableExtra::row_spec(which_row_spec, bold = TRUE, background = "#d7ecff") %>% 
    kableExtra::row_spec(which_row_spec_isna, background = "#f0f8ff", color = "#666666")
  
  
  
  return(kable)
  
})





setMethod("show", "BclassCharacteristics", function(object){
  
  if(interactive()){
    print(bkable(object))
  }else{
    cat(bkable(object))
  }
  
})

























