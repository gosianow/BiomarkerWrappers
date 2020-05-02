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
### Show method and bkable
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("bkable", "BclassTesting", function(x, caption = NULL, header = NULL, block_vars = NULL, font_size = NULL, full_width = NULL){
  
  
  res <- bresults(x)
  out <- boutput(x)
  
  
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
    kableExtra::kable_styling(bootstrap_options = c("condensed", "bordered"), latex_options = c("HOLD_position"), full_width = full_width, font_size = font_size) %>% 
    kableExtra::column_spec(which(colnames(out) %in% c("HR", "OR", "Difference", "P-value", "Adj. P-value")), bold = TRUE) 
  
  
  if(!is.null(header)){
    kable <- kable %>% 
      add_header_above(header)
  }
  
  
  
  
  
  # --------------------------------------------------------------------------
  ## Define rows that should be colored
  # --------------------------------------------------------------------------
  
  ### By default color per covariate/biomarker block
  if(is.null(block_vars)){
    
    block_vars <- colnames(res)[which(colnames(res) %in% c("covariate", "biomarker", "covariate1"))]
    
  }
  
  
  
  ## If more than one covariate/biomarker, identify every second changing character group
  which_row_spec <- indicate_blocks(d = res, block_vars = block_vars, return = "block")
  
  if(is.null(which_row_spec)){
    ## If one covariate/biomarker, make all rows colored
    which_row_spec <- as.numeric(1:nrow(res))
  }
  
  
  kable <- kable %>% 
    kableExtra::row_spec(which_row_spec, background = "#F2F4F4")
  
  
  
  return(kable)
  
})



setMethod("show", "BclassTesting", function(object){
  
  print(bkable(object))
  
})



################################################################################
### bforest method
################################################################################


#' @rdname Bclass-class
#' @export
setMethod("bforest", "BclassTesting", function(x, mean_var = NULL, lower_var = NULL, upper_var = NULL, block_vars = NULL, xlab = NULL, clip = c(0, 5), lineheight = "auto"){
  
  # lineheight = unit(1, "cm")
  
  
  out <- boutput(x)
  res <- bresults(x)
  caption <- bcaption(x)
  
  ### ----------------------------------------------------------------------
  ### Some checks
  ### ----------------------------------------------------------------------
  
  if(is.null(mean_var)){
    mean_var <- colnames(res)[which(colnames(res) %in% c("difference", "HR", "OR"))]
  }
  stopifnot(length(mean_var) == 1)
  stopifnot(mean_var %in% colnames(res))
  
  if(is.null(lower_var)){
    lower_var <- colnames(res)[which(colnames(res) %in% c("CI95_lower"))]
  }
  stopifnot(length(lower_var) == 1)
  stopifnot(lower_var %in% colnames(res))
  
  if(is.null(upper_var)){
    upper_var <- colnames(res)[which(colnames(res) %in% c("CI95_upper"))]
  }
  stopifnot(length(upper_var) == 1)
  stopifnot(upper_var %in% colnames(res))
  
  
  if(is.null(xlab)){
    xlab <- mean_var
  }
  
  
  ### ----------------------------------------------------------------------
  
  res[res[, upper_var] %in% Inf, upper_var] <- NA
  
  labeltext <- rbind(colnames(out), out)
  
  
  ### ----------------------------------------------------------------------
  ### To separate Biomarkers with a horizontal line
  
  
  ## By default color per covariate/biomarker block
  if(is.null(block_vars)){
    block_vars <- colnames(res)[which(colnames(res) %in% c("covariate", "biomarker", "covariate1"))]
  }
  stopifnot(all(block_vars %in% colnames(res)))
  
  
  
  ## Find numbers of rows where the biomarker blocks end. This value has to be shifted by 2
  line_row <- indicate_blocks(res, block_vars = block_vars, return = "line") + 2
  
  line_row <- line_row[-length(line_row)]
  
  
  hrzl_lines <- list()
  
  if(length(line_row) >= 1){
    hrzl_lines <- lapply(line_row, function(x){
      gpar(col = "#b4b4b4", lwd = 0.5)
    })
    names(hrzl_lines) <- line_row
  }
  
  hrzl_lines[["2"]] <- gpar(col = "#444444", lwd = 1)
  
  
  ### ----------------------------------------------------------------------
  
  
  
  forestplot::forestplot(labeltext, mean = c(NA, res[, mean_var]), lower = c(NA, res[, lower_var]), upper = c(NA, res[, upper_var]), 
    is.summary = c(TRUE, rep(FALSE, nrow(res))), xlog = FALSE, xlab = xlab, zero = 1,
    title = caption,
    col = fpColors(box = "darkblue", line = "darkblue"), 
    boxsize = 0.3,
    hrzl_lines = hrzl_lines, 
    graphwidth = unit(10, "cm"), colgap = unit(6, "mm"),
    lineheight = lineheight,
    lwd.ci = 2, lwd.xaxis = 2, lwd.zero = 2, 
    txt_gp = forestplot::fpTxtGp(xlab = gpar(fontsize = 20), ticks = gpar(fontsize = 18)), 
    mar = unit(c(20, rep(5, times = 3)), "mm"), 
    clip = clip, 
    ci.vertices = TRUE)
  
  
})






















