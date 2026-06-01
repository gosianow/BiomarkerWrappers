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
  
  ### We add format as an option because kableExtra::column_spec does not work if format is not specified in the kable call
  format <- getOption("knitr.table.format", default = "html")
  
  ### Make unique labels to have the right numbering of the tables. See https://github.com/rstudio/bookdown/issues/1097
  
  kable <- knitr::kable(out, format = format, caption = caption, booktabs = TRUE, linesep = "", row.names = FALSE, label = gsub("\\.", "", make.names(Sys.time()))) %>%
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
  
  if(interactive()){
    print(bkable(object))
  }else{
    cat(bkable(object))
  }
  
})



################################################################################
### bforest method
################################################################################



#' @rdname Bclass-class
#' @export
setMethod("bforest", "BclassTesting", function(x, mean_var = NULL, lower_var = NULL, upper_var = NULL, block_vars = NULL, xlab = NULL, xlog = FALSE, clip = NULL, xticks = NULL, lineheight = "auto", label_fontsize = 14, caption_width = 150, shapes_gp = forestplot::fpShapesGp()){
  
  # ------------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------------
  
  pick_one <- function(res_colnames, candidates, label){
    hit <- candidates[candidates %in% res_colnames]
    if(length(hit) < 1){
      stop(paste0("Could not infer ", label, ". None of these columns found: ", paste(candidates, collapse = ", ")), call. = FALSE)
    }
    if(length(hit) > 1){
      stop(paste0("Multiple matches for ", label, ": ", paste(candidates, collapse = ", ")), call. = FALSE)
    }
    hit
  }
  
  log_ticks_125 <- function(lims){
    if(any(!is.finite(lims)) || any(lims <= 0)){
      stop("For xlog=TRUE, limits must be finite and > 0.", call. = FALSE)
    }
    e <- seq(floor(log10(lims[1])), ceiling(log10(lims[2])))
    sort(unique(as.vector(outer(c(1, 2, 5), 10^e, `*`))))
  }
  
  cap_range <- function(res_range, clip){
    c(max(res_range[1], clip[1]), min(res_range[2], clip[2]))
  }
  
  # ------------------------------------------------------------------
  # Pull data from object
  # ------------------------------------------------------------------
  
  out <- boutput(x)
  res <- bresults(x)
  caption <- bcaption(x)
  
  if(!is.null(caption_width)){
    caption <- stringr::str_wrap(caption, width = caption_width)
  }
  
  res_colnames <- colnames(res)
  
  # ------------------------------------------------------------------
  # Infer columns
  # ------------------------------------------------------------------
  
  if(is.null(mean_var)){
    mean_var <- pick_one(res_colnames, c("HR", "OR", "difference"), "mean_var")
  }
  if(is.null(lower_var)){
    lower_var <- pick_one(res_colnames, c("HR_CI95_lower", "OR_CI95_lower", "difference_CI95_lower"), "lower_var")
  }
  if(is.null(upper_var)){
    upper_var <- pick_one(res_colnames, c("HR_CI95_upper", "OR_CI95_upper", "difference_CI95_upper"), "upper_var")
  }
  
  
  # Validate they exist and are scalar
  stopifnot(length(mean_var) == 1, length(lower_var) == 1, length(upper_var) == 1)
  stopifnot(all(c(mean_var, lower_var, upper_var) %in% res_colnames))
  
  if(is.null(xlab)){
    xlab <- mean_var
  } 
  
  is_ratio <- mean_var %in% c("HR", "OR")
  
  if(is_ratio){
    zero <- 1
  }else{
    zero <- 0  
  }
  
  # ------------------------------------------------------------------
  # Clip defaults
  # ------------------------------------------------------------------
  
  # If xticks supplied, let them define clip unless user explicitly set clip
  
  if(!is.null(xticks) && is.null(clip)){
    clip <- range(xticks, finite = TRUE)
  }
  
  if(is.null(clip)){
    if(is_ratio){
      clip <- if(xlog) c(0.1, 10) else c(0, 4)
    } else {
      clip <- c(-40, 40)
    }
  }
  
  if(xlog && any(clip <= 0)){
    stop("For xlog=TRUE, clip must be strictly > 0.", call. = FALSE)
  }
  
  if(!is.null(xticks) && xlog && any(xticks <= 0)){
    stop("For xlog=TRUE, xticks must be strictly > 0.", call. = FALSE)
  }
  
  # ------------------------------------------------------------------
  # Build xticks if missing
  # ------------------------------------------------------------------
  
  if(is.null(xticks)){
    
    res_range <- range(c(res[[lower_var]], res[[upper_var]]), na.rm = TRUE, finite = TRUE)
    
    res_range <- cap_range(res_range, clip)
    
    center <- if(is_ratio) 1 else 0
    
    if(!xlog){
      # Use pretty breaks; then ensure center is included
      xticks <- sort(unique(c(center, pretty(res_range, n = 7))))
    } else {
      xticks <- log_ticks_125(res_range)
      xticks <- xticks[xticks >= clip[1] & xticks <= clip[2]]
      xticks <- sort(unique(c(1, xticks)))
    }
  }
  
  
  # ------------------------------------------------------------------
  # Prepare plotting data (don’t mutate res)
  # ------------------------------------------------------------------
  
  res_plot <- res
  
  # Clip CIs and mean to plotting range
  
  if(is_ratio){

    offset_upper <- clip[2] / 100
    offset_lower <- clip[1] / 100
    
    # res_plot[[upper_var]] <- ifelse(!is.na(res_plot[[upper_var]]) & res_plot[[upper_var]] > clip[2], clip[2], res_plot[[upper_var]])
    res_plot[[upper_var]] <- ifelse(is.na(res_plot[[upper_var]]) | res_plot[[upper_var]] > clip[2] + offset_upper, clip[2] + offset_upper, res_plot[[upper_var]])
    # res_plot[[lower_var]] <- ifelse(!is.na(res_plot[[lower_var]]) & res_plot[[lower_var]] < clip[1], clip[1], res_plot[[lower_var]])
    res_plot[[lower_var]] <- ifelse(is.na(res_plot[[lower_var]]) | res_plot[[lower_var]] < clip[1] - offset_lower, clip[1] - offset_lower, res_plot[[lower_var]])
    res_plot[[mean_var]] <- ifelse(!is.na(res_plot[[mean_var]]) & res_plot[[mean_var]] > clip[2], clip[2], res_plot[[mean_var]])
    res_plot[[mean_var]] <- ifelse(!is.na(res_plot[[mean_var]]) & res_plot[[mean_var]] < clip[1], clip[1], res_plot[[mean_var]])

  }else{
    # set +/-Inf to NA
    # res_plot[[upper_var]][is.infinite(res_plot[[upper_var]])] <- NA
    # res_plot[[lower_var]][is.infinite(res_plot[[lower_var]])] <- NA
    # res_plot[[mean_var]][is.infinite(res_plot[[mean_var]])] <- NA
    
    offset <- (clip[2] - clip[1]) / 100
    
    # res_plot[[upper_var]] <- ifelse(!is.na(res_plot[[upper_var]]) & res_plot[[upper_var]] > clip[2], clip[2], res_plot[[upper_var]])
    res_plot[[upper_var]] <- ifelse(is.na(res_plot[[upper_var]]) | res_plot[[upper_var]] > clip[2] + offset, clip[2] + offset, res_plot[[upper_var]])
    # res_plot[[lower_var]] <- ifelse(!is.na(res_plot[[lower_var]]) & res_plot[[lower_var]] < clip[1], clip[1], res_plot[[lower_var]])
    res_plot[[lower_var]] <- ifelse(is.na(res_plot[[lower_var]]) | res_plot[[lower_var]] < clip[1] - offset, clip[1] - offset, res_plot[[lower_var]])
    res_plot[[mean_var]] <- ifelse(!is.na(res_plot[[mean_var]]) & res_plot[[mean_var]] > clip[2], clip[2], res_plot[[mean_var]])
    res_plot[[mean_var]] <- ifelse(!is.na(res_plot[[mean_var]]) & res_plot[[mean_var]] < clip[1], clip[1], res_plot[[mean_var]])
    
  }
  
  # ------------------------------------------------------------------
  # Label text
  # ------------------------------------------------------------------
  
  # Ensure matrix-like structure for forestplot labeltext
  labeltext <- rbind(colnames(out), as.matrix(out))
  
  # ------------------------------------------------------------------
  # Horizontal lines for blocks
  # ------------------------------------------------------------------
  
  
  if(is.null(block_vars)){
    # preserve your original default candidates, but deterministic priority
    # (take all that exist)
    candidates <- c("covariate", "biomarker", "covariate1")
    block_vars <- candidates[candidates %in% colnames(res_plot)]
  }
  stopifnot(length(block_vars) >= 1)
  stopifnot(all(block_vars %in% colnames(res_plot)))
  
  line_row <- indicate_blocks(res_plot, block_vars = block_vars, return = "line") + 2
  if(length(line_row) > 0) line_row <- line_row[-length(line_row)]  # drop last
  
  hrzl_lines <- list()
  if(length(line_row) >= 1){
    hrzl_lines <- setNames(
      rep(list(grid::gpar(col = "#b4b4b4", lwd = 0.5)), length(line_row)),
      as.character(line_row)
    )
  }
  # thick line under header row
  hrzl_lines[["2"]] <- grid::gpar(col = "#444444", lwd = 1)
  
  
  # ------------------------------------------------------------------
  # Plot
  # ------------------------------------------------------------------
  
  
  p <- forestplot::forestplot(labeltext,
    mean  = c(NA, res_plot[[mean_var]]),
    lower = c(NA, res_plot[[lower_var]]),
    upper = c(NA, res_plot[[upper_var]]),
    is.summary = c(TRUE, rep(FALSE, nrow(res_plot))),
    xlab = xlab,
    zero = zero,
    title = caption,
    col = forestplot::fpColors(box = "darkblue", line = "darkblue"),
    boxsize = 0.4,
    hrzl_lines = hrzl_lines,
    graphwidth = grid::unit(10, "cm"),
    colgap = grid::unit(6, "mm"),
    lineheight = lineheight,
    lwd.ci = 2,
    lwd.xaxis = 2,
    lwd.zero = 2.5,
    txt_gp = forestplot::fpTxtGp(
      label = grid::gpar(fontsize = label_fontsize),
      xlab  = grid::gpar(fontsize = 24),
      ticks = grid::gpar(fontsize = 22),
      title = grid::gpar(fontsize = label_fontsize, fontface = "bold")
    ),
    mar = grid::unit(c(5, 5, 5, 5), "mm"),
    clip = clip,
    ci.vertices = TRUE,
    align = "l",
    xticks = xticks,
    xticks.digits = 2,
    xlog = xlog,
    shapes_gp = shapes_gp
  )
  
  print(p)
  
  invisible(p)
  
})




