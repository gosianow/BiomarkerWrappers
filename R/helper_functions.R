

#' Use round or signif
#' 
#' @keywords internal
round_signif <- function(x, digits = 2){
  ifelse(abs(x) >= 1, round(x, digits = digits), signif(x, digits = digits))
}




#' Save content of an eSet into a CSV file 
#' 
#' @export
eSet2csv <- function(es, file, digits = 2){
  
  fdata <- Biobase::fData(es)
  pdata <- Biobase::pData(es)
  expr <- round_signif(Biobase::exprs(es), digits = digits)
  
  
  if(!all(colnames(expr) == pdata[, 1])){
    stop("Column names of expression must be the same as the fisrt column in pData.")
  }
  
  
  ## Prepare the top table
  if(ncol(pdata) == 0){
    top_table <- data.frame()
  }else{
    top_table <- cbind(matrix("", nrow = ncol(pdata), ncol = ncol(fdata)), colnames(pdata), t(pdata))
  }
  
  utils::write.table(top_table, file = file, sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE, append = FALSE)
  
  ## Add empty row
  if(nrow(top_table) > 0){
    utils::write.table("", file = file, sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE, append = TRUE)
  }
  
  ## Prepare bottom table
  
  bottom_header <- t(c(colnames(fdata), "", colnames(expr)))
  
  utils::write.table(bottom_header, file = file, sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  bottom_table <- data.frame(fdata, "", as.data.frame.matrix(expr), stringsAsFactors = FALSE, check.names = FALSE)
  
  utils::write.table(bottom_table, file = file, sep = ",", quote = TRUE, row.names = FALSE, col.names = FALSE, append = TRUE)
  
  
}






### Based on this AWESOME post http://michaeljw.com/blog/post/subchunkify/

# subchunkify <- function(g, fig_height=7, fig_width=5) {
#   g_deparsed <- paste0(deparse(
#     function() {g}
#   ), collapse = '')
#   
#   sub_chunk <- paste0("
#   `","``{r, fig.height=", fig_height, ", fig.width=", fig_width, ", echo=FALSE}",
#     "\n(", 
#     g_deparsed
#     , ")()",
#     "\n`","``
#   ")
#   
#   cat(knitr::knit(text = knitr::knit_expand(text = sub_chunk), quiet = TRUE))
# }





#' My version of subchunkify 
#' 
#' @export
subchunkify <- function(g, fig.height = 5, fig.width = 6, chunk_name = NULL, envir = parent.frame(), display_subchunk = FALSE){
  ### Necessary to set for knit_child
  ## knitr::opts_knit$set(output.dir = getwd())
  
  # -----------------------------------------------------------------------
  # Some checks
  # -----------------------------------------------------------------------
  
  stopifnot(length(fig.height) == 1)
  stopifnot(fig.height > 0)
  
  stopifnot(length(fig.width) == 1)
  stopifnot(fig.width > 0)
  
  # -----------------------------------------------------------------------
  # Subchunkify
  # -----------------------------------------------------------------------
  
  g_deparsed <- paste0(deparse(substitute(g)), collapse = "")
  
  ### Make unique labels. In the end when using knit_child no need for unique labels. However, chunk_name must be unique.
  # time <- paste0(gsub("[[:punct:] ]", "_", format(Sys.time())), "_", floor(runif(1) * 10000))
  time <- NULL
  
  sub_chunk <- paste0("\n```{r ", chunk_name, time, ", fig.height=", fig.height, ", fig.width=", fig.width, ", echo=FALSE}\n",
    g_deparsed,
    "\n```\n")
  
  if(display_subchunk){
    cat(sub_chunk)
  }
  
  ### Using knit_child is essential. Otherwise, with knit() error about unique chunk names.
  
  cat(knitr::knit_child(text = sub_chunk, quiet = TRUE, envir = envir))
  
}









#' Calculate logFC for expression data that will be plotted in a heatmap 
#' 
#' @param x eSet object
#' @export
wrapper_calculate_sample_logFC <- function(x, comparison_var, subgroup_var = NULL){
  
  
  if(is.null(subgroup_var)){
    Biobase::pData(x)$dummy_subgroup_var <- factor("dummy_subgroup_var")
    subgroup_var <- "dummy_subgroup_var"
  }
  
  subgroup_levels <- levels(Biobase::pData(x)[, subgroup_var])
  
  comparison_levels <- levels(Biobase::pData(x)[, comparison_var])
  
  reference_level <- comparison_levels[1]
  
  tbl <- table(Biobase::pData(x)[, comparison_var])
  
  stopifnot(tbl[reference_level] > 0)
  
  
  expr_heatmap <- lapply(1:length(subgroup_levels), function(i){
    # i = 1
    
    x_sub <- x[, Biobase::pData(x)[, subgroup_var] == subgroup_levels[i]]
    
    expr_sub <- Biobase::exprs(x_sub)
    
    expr_lfc <- expr_sub - rowMeans(expr_sub[, Biobase::pData(x_sub)[, comparison_var] == reference_level], na.rm = TRUE)
    
    expr_lfc
    
    
  }) 
  
  
  expr_heatmap <- do.call("cbind", expr_heatmap)
  
  
  ### Bring the same sample order as in x
  
  expr_heatmap <- expr_heatmap[, colnames(Biobase::exprs(x)), drop = FALSE]
  
  ### Return an eSet
  
  out <- Biobase::ExpressionSet(assayData = expr_heatmap, 
    phenoData = Biobase::phenoData(x), 
    featureData = Biobase::featureData(x))
  
  
  out
  
}


#' My version of duplicated
#' 
#' @export
duplicated2 <- function(x, value = FALSE, na.rm = FALSE){
  
  if(is.factor(x)){
    x <- as.character(x)
  }
  
  ### Make the NAs to count as unique
  if(na.rm){
    x[is.na(x)] <- paste0("X...NA...", 1:sum(is.na(x)))
  }
  
  ifduplicated <- duplicated(x, fromLast = FALSE) | duplicated(x, fromLast = TRUE)
  
  if(value){
    out <- x[ifduplicated]
  }else{
    out <- ifduplicated
  }
  
  return(out)
  
}



#' My version of order for data frame
#' 
#' @param x Data frame
#' @export
order2 <- function(x, na.last = TRUE, decreasing = FALSE){
  
  stopifnot(is.data.frame(x))
  
  n <- ncol(x)
  
  if(length(na.last) == 1){
    na.last <- rep(na.last, n)
  }else{
    stopifnot(length(na.last) == n)  
  }
  
  if(length(decreasing) == 1){
    decreasing <- rep(decreasing, n)
  }else{
    stopifnot(length(decreasing) == n)  
  }
  
  
  new_order <- seq_len(nrow(x))
  
  for(i in rev(seq_len(n))){
    # i = 1
    
    oo <- order(x[new_order, i], na.last = na.last[i], decreasing = decreasing[i])
    
    new_order <- new_order[oo]
    
  }
  
  return(new_order)
  
}







#' Calculate number of non-NA elements in a vector 
#' 
#' @keywords internal
length_nonNA <- function(x){
  sum(!is.na(x))
}




#' Read gmt file and return a list of genes
#' 
#' @export
wrapper_read_gmt <- function(filename){
  x <- GSEABase::getGmt(filename)
  x <- GSEABase::geneIds(x)
  return(x)
}


#' Write a list of genes into a gmt file
#' 
#' @export
wrapper_write_gmt <- function(x, filename){
  
  gsc <- GSEABase::GeneSetCollection(lapply(1:length(x), function(i){
    # i = 1
    GSEABase::GeneSet(setName = names(x[i]), geneIds = x[[i]])
    
  }))
  
  GSEABase::toGmt(gsc, filename)
  
  invisible(NULL)
  
}









#' Spit a list of plots into chunks and plot them in a grid layout
#' 
#' @export
wrapper_print_plot_grid <- function(plotlist, nsplit = NULL, ncol = 2, nrow = NULL){
  
  if(is.null(nsplit)){
    nsplit <- length(plotlist)
  }
  
  indx <- seq_along(plotlist)
  indx_split <- split(indx, ceiling(seq_along(indx) / nsplit))
  
  for(i in seq_along(indx_split)){
    print(plot_grid(plotlist = plotlist[indx_split[[i]]], ncol = ncol, nrow = nrow))
  }
  
  # invisible()
  
}





#' Stratify data into quartiles
#' 
#' @param x Vector of continuous values to cut.
#' 
#' @export
wrapper_cut_quartiles <- function(x, labels = c("[0%, 25%]", "(25%, 50%]", "(50%, 75%]", "(75%, 100%]")){
  
  stopifnot(length(labels) == 4)
  
  out <- NULL
  
  try(out <- ggplot2::cut_number(x, n = 4), silent = TRUE)
  
  if(!is.null(out)){
    out <- factor(out, labels = labels)
  }else{
    out <- rep(NA, length(x))
  }
  
  return(out)
  
}


#' Dichotomize data by median 
#' 
#' @param x Vector of continuous values to cut.
#' 
#' @export
wrapper_cut_median <- function(x, labels = c("<=MED", ">MED")){
  
  stopifnot(length(labels) == 2)
  
  out <- ggplot2::cut_number(x, n = 2)
  
  out <- factor(out, labels = labels)
  
  return(out)
  
}

#' Stratify data into two groups
#' 
#' @param x Vector of continuous values to cut.
#' 
#' @export
wrapper_cut_2groups <- function(x, probs = 0.5, cutoff = NULL, labels = c("low", "high")){
  
  stopifnot(length(probs) == 1)
  stopifnot(length(labels) == 2)
  
  if(is.null(cutoff)){
    cutoff <- stats::quantile(x, probs = probs, na.rm = TRUE)
  }
  
  stopifnot(length(cutoff) == 1)
  
  out <- factor(ifelse(x <= cutoff, labels[1], labels[2]), levels = labels)
  
  return(out)
  
}


#' @rdname wrapper_cut_quartiles
#' @export
wrapper_cut_quartiles_strat <- function(x, strata, labels = c("[0%, 25%]", "(25%, 50%]", "(50%, 75%]", "(75%, 100%]")){
  
  stopifnot(is.factor(strata))
  
  strata_levels <- levels(strata)
  
  out <- factor(rep(NA, length(x)), levels = labels)
  
  for(i in 1:length(strata_levels)){
    # i = 1
    
    indx_sub <- which(strata %in% strata_levels[i])
    
    out[indx_sub] <- wrapper_cut_quartiles(x[indx_sub], labels = labels)
    
  }
  
  return(out)
  
}



#' @rdname wrapper_cut_median
#' @export
wrapper_cut_median_strat <- function(x, strata, labels = c("<=MED", ">MED")){
  
  stopifnot(is.factor(strata))
  
  strata_levels <- levels(strata)
  
  out <- factor(rep(NA, length(x)), levels = labels)
  
  for(i in 1:length(strata_levels)){
    # i = 1
    
    indx_sub <- which(strata %in% strata_levels[i])
    
    out[indx_sub] <- wrapper_cut_median(x[indx_sub], labels = labels)
    
  }
  
  return(out)
  
}



#' @rdname wrapper_cut_2groups
#' @export
wrapper_cut_2groups_strat <- function(x, strata, probs = rep(0.5, nlevels(strata)), cutoff = NULL, labels = c("low", "high")){
  
  stopifnot(is.factor(strata))
  
  stopifnot(length(probs) == nlevels(strata))
  
  if(!is.null(cutoff)){
    stopifnot(length(cutoff) == nlevels(strata))
  }
  
  strata_levels <- levels(strata)
  
  out <- factor(rep(NA, length(x)), levels = labels)
  
  for(i in seq_along(strata_levels)){
    # i = 1
    
    indx_sub <- which(strata %in% strata_levels[i])
    
    out[indx_sub] <- wrapper_cut_2groups(x[indx_sub], probs = probs[i], cutoff = cutoff[i], labels = labels) 
    
  }
  
  return(out)
  
}









#' Internal function used in bkable to highlight the rows 
#' 
#' @keywords internal
indicate_blocks <- function(d, block_vars, return = "block"){
  
  stopifnot(length(return) == 1)
  stopifnot(return %in% c("block", "line", "empty"))
  
  
  all_data <- apply(d[, block_vars, drop = FALSE], 1, paste, collapse = ".")
  
  all_rle <- c(0, cumsum(rle(all_data)$lengths))
  
  
  if(return == "block"){
    
    if(length(all_rle) > 2){
      
      max_rle_indx <- ifelse(length(all_rle) %% 2 == 0, length(all_rle), length(all_rle) - 1)
      
      out <- unlist(lapply(seq(1, max_rle_indx, by = 2), function(i){
        seq(all_rle[i] + 1, all_rle[i + 1], by = 1)
      }))
      
    }else{
      out <- NULL
    }
    
    
  }else if(return == "line"){
    
    out <- all_rle[-1]
    
    
  }else if(return == "empty"){
    
    ## Indicate first lines that are unique
    out <- c(1, all_rle[-1] + 1)
    out <- out[-length(out)]
    
    ## But make empty all besides them
    
    out <- setdiff(1:nrow(d), out)
    
    
  }
  
  
  return(out)
  
}






#' Format header for kable
#' 
#' @param all_colnames All column names 
#' @param header_colnames Column names where the header should be placed.
#' @param header_name Name of the header.
#' 
#' @examples 
#' 
#' \dontrun{
#' all_colnames <- c("Covariate", "Subgroup", "A", "B", "C", "OR", "P-value")
#' header_colnames <- c("A", "B", "C")
#' header_name <- "Gene A"
#' 
#' format_header(all_colnames, header_colnames, header_name)
#' 
#' 
#' all_colnames <- c("A", "B", "C", "OR", "P-value")
#' header_colnames <- c("A", "B", "C")
#' header_name <- "Gene A"
#' 
#' format_header(all_colnames, header_colnames, header_name)
#' }
#' @keywords internal
format_header <- function(all_colnames, header_colnames, header_name){
  
  
  header_indx <- which(all_colnames %in% header_colnames)
  num_mid_cols <- length(header_indx)
  stopifnot(num_mid_cols > 0)
  
  
  if(header_indx[1] > 1){
    
    num_start_cols <- header_indx[1] - 1
    
    if(header_indx[num_mid_cols] < length(all_colnames)){
      
      num_end_cols <- length(all_colnames) - header_indx[num_mid_cols]
      
      header <- c(num_start_cols, num_mid_cols, num_end_cols)
      header <- as.integer(header)
      names(header) <- c(" ", header_name, " ") 
      
    }else{
      
      header <- c(num_start_cols, num_mid_cols)
      header <- as.integer(header)
      names(header) <- c(" ", header_name)
      
    }
    
  }else{
    
    
    if(header_indx[num_mid_cols] < length(all_colnames)){
      
      num_end_cols <- length(all_colnames) - header_indx[num_mid_cols]
      
      header <- c(num_mid_cols, num_end_cols)
      header <- as.integer(header)
      names(header) <- c(header_name, " ") 
      
    }else{
      
      header <- c(num_mid_cols)
      header <- as.integer(header)
      names(header) <- c(header_name)
      
    }
    
    
  }
  
  stopifnot(sum(header) == length(all_colnames))
  
  header
  
  
}






#' Format p-values
#' 
#' @param x Vector of p-values to be formatted.
#' @param digits Number of digits after decimal to display.
#' @param asterisk Logical, whether to indicate significance levels with asterisks. Encoding: `***` p-value < 0.001, `**` p-value < 0.01, `*` p-value < 0.05, `.` p-value < 0.1.
#' @param non_empty Vector defining which values should be displayed despite being NAs.
#' @keywords internal 
format_pvalues <- function(x, digits = 4, asterisk = TRUE, non_empty = NULL){
  
  
  if(!is.null(non_empty)){
    if(is.logical(non_empty)){
      stopifnot(length(non_empty) == length(x))
    }else{
      non_empty_logical <- rep(FALSE, length(x))
      non_empty_logical[non_empty] <- TRUE
      non_empty <- non_empty_logical
    }
  }
  
  
  # digits = 4
  # asterisk <- TRUE
  # x <- c(0.2, 0.05, 0.034534, 1.366332e-05, 1.366332e-04, NA)
  
  
  if(sum(is.na(x)) == length(x) && is.null(non_empty)){
    output <- rep("", length(x))
    return(output)
  }else if(sum(is.na(x)) == length(x) && !is.null(non_empty)){
    output <- rep("", length(x))
    output[non_empty] <- "NA"
    return(output)
  }
  
  min_pval <- 1/10^digits
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = FALSE)
  output[x < min_pval] <- paste0("<", formatC(min_pval, format = "f", digits = digits))
  output[is.na(x)] <- "NA"
  
  if(is.null(non_empty)){
    output[output == "NA"] <- ""
  }else{
    output[output == "NA" & !non_empty] <- ""
  }
  
  
  if(asterisk){
    
    pval_asterisk <- ifelse(x < 0.001, " ***", ifelse(x < 0.01, " **", ifelse(x < 0.05, " *", ifelse(x < 0.1, " .", ""))))
    
    pval_asterisk[is.na(pval_asterisk)] <- ""
    
    output <- paste0(output, pval_asterisk)
    
  }
  
  return(output)
  
}



#' Format p-values using scientific format
#' 
#' @param x Vector of p-values to be formatted.
#' @param digits Number of digits after decimial to display.
#' @param asterisk Logical, whether to indicate significance levels with asterisks. Encoding: `***` p-value < 0.001, `**` p-value < 0.01, `*` p-value < 0.05, `.` p-value < 0.1.
#' @keywords internal 
format_pvalues2 <- function(x, digits = 4, asterisk = TRUE){
  
  
  
  # digits = 4
  # asterisk <- TRUE
  # x <- c(0.2, 0.05, 0.034534, 1.366332e-05, 1.366332e-04, NA, 7.174163e-16, 1.501826e-06, 6.642127e-10)
  
  
  if(sum(is.na(x)) == length(x)){
    return(rep("", length(x)))
  }
  
  
  min_pval <- 1/10^digits
  
  output_format_non_scientific <- formatC(x, format = "f", digits = digits, drop0trailing = FALSE)
  
  output_format_scientific <- paste0("<", formatC(10 ^ -round(-log10(x)), format = "e", digits = 0))
  
  
  output <- ifelse(x < min_pval, output_format_scientific, output_format_non_scientific)
  
  output[is.na(x)] <- ""
  
  
  
  if(asterisk){
    
    pval_asterisk <- ifelse(x < 0.001, " ***", ifelse(x < 0.01, " **", ifelse(x < 0.05, " *", ifelse(x < 0.1, " .", ""))))
    
    pval_asterisk[is.na(pval_asterisk)] <- ""
    
    output <- paste0(output, pval_asterisk)
    
  }
  
  
  return(output)
  
}







#' Format Odds Ratios
#' 
#' @param x Vector with odds ratios.
#' @param digits Number of decimal places.
#' @param non_empty Vector defining which values should be displayed despite being NAs.
#' @keywords internal 
format_or <- function(x, digits = 2, non_empty = NULL){
  
  if(!is.null(non_empty)){
    if(is.logical(non_empty)){
      stopifnot(length(non_empty) == length(x))
    }else{
      non_empty_logical <- rep(FALSE, length(x))
      non_empty_logical[non_empty] <- TRUE
      non_empty <- non_empty_logical
    }
  }
  
  
  if(sum(is.na(x)) == length(x) && is.null(non_empty)){
    output <- rep("", length(x))
    return(output)
  }else if(sum(is.na(x)) == length(x) && !is.null(non_empty)){
    output <- rep("", length(x))
    output[non_empty] <- "NA"
    return(output)
  }
  
  min_val <- 1/10^digits
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = FALSE)
  output[x < min_val] <- paste0("<", formatC(min_val, format = "f", digits = digits))
  output[x %in% 0] <- "0"
  output[is.na(x)] <- "NA"
  
  if(is.null(non_empty)){
    output[output == "NA"] <- ""
  }else{
    output[output == "NA" & !non_empty] <- ""
  }
  
  
  return(output)
  
}



#' Format difference 
#' 
#' @param x Vector with differences.
#' @param digits Number of decimal places.
#' @param non_empty Vector defining which values should be displayed despite being NAs.
#' @keywords internal
format_difference <- function(x, digits = 2, non_empty = NULL){
  
  
  if(!is.null(non_empty)){
    if(is.logical(non_empty)){
      stopifnot(length(non_empty) == length(x))
    }else{
      non_empty_logical <- rep(FALSE, length(x))
      non_empty_logical[non_empty] <- TRUE
      non_empty <- non_empty_logical
    }
  }
  
  if(sum(is.na(x)) == length(x) && is.null(non_empty)){
    output <- rep("", length(x))
    return(output)
  }else if(sum(is.na(x)) == length(x) && !is.null(non_empty)){
    output <- rep("", length(x))
    output[non_empty] <- "NA"
    return(output)
  }
  
  min_val <- 1/10^digits
  
  output <- formatC(x, format = "f", digits = digits, drop0trailing = FALSE)
  output[is.na(x)] <- "NA"
  
  if(is.null(non_empty)){
    output[output == "NA"] <- ""
  }else{
    output[output == "NA" & !non_empty] <- ""
  }
  
  return(output)
  
}



#' Format CIs (confidence intervals)
#' 
#' @param CI_lower Vector with lower CIs.
#' @param CI_upper Vector with upper CIs.
#' @param digits Number of decimal places.
#' @param non_empty Vector defining which values should be displayed despite being NAs.
#' @keywords internal
format_CIs <- function(CI_lower, CI_upper, digits = 2, non_empty = NULL){
  
  stopifnot(length(CI_lower) == length(CI_upper))
  
  if(!is.null(non_empty)){
    if(is.logical(non_empty)){
      stopifnot(length(non_empty) == length(CI_lower))
    }else{
      non_empty_logical <- rep(FALSE, length(CI_lower))
      non_empty_logical[non_empty] <- TRUE
      non_empty <- non_empty_logical
    }
  }
  
  if(sum(is.na(CI_lower)) == length(CI_lower) && sum(is.na(CI_upper)) == length(CI_upper) && is.null(non_empty)){
    output <- rep("", length(CI_lower))
    return(output)
  }else if(sum(is.na(CI_lower)) == length(CI_lower) && sum(is.na(CI_upper)) == length(CI_upper) && !is.null(non_empty)){
    output <- rep("", length(CI_lower))
    output[non_empty] <- "NA"
    return(output)
  }
  
  output <- paste0("(", formatC(CI_lower, format = "f", digits = digits, drop0trailing = FALSE), " - ", formatC(CI_upper, format = "f", digits = digits, drop0trailing = FALSE), ")")
  
  output[is.na(CI_lower) & is.na(CI_upper)] <- "NA"
  
  if(is.null(non_empty)){
    output[output == "NA"] <- ""
  }else{
    output[output == "NA" & !non_empty] <- ""
  }
  
  return(output)
  
}


#' Format CIs (confidence intervals)
#' 
#' @param x Data frame.
#' @param digits Number of decimal places.
#' @param colnames New colnames.
#' @param non_empty Vector defining which values should be displayed despite being NAs.
#' @keywords internal
format_CIs_df <- function(x, digits = 2, colnames = NULL, non_empty = NULL){
  
  output <- data.frame(format_CIs(x[, 1], x[, 2], digits = digits, non_empty = non_empty), stringsAsFactors = FALSE)
  colnames(output) <- colnames
  
  return(output)
  
}




#' Format versus
#' 
#' @param level Vector with levels.
#' @param reference Vector with references.
#' @keywords internal
format_vs <- function(level, reference){
  
  output <- paste0(level, " vs ", reference)
  
  output[level == "" & reference == ""] <- ""
  output[level == 0 & reference == 0] <- ""
  output[is.na(level) & is.na(reference)] <- ""
  
  
  return(output)
  
}



#' Format proportions
#' 
#' @param props Vector with proportions.
#' @param digits Number of decimal places when rounding proportions.
#' @keywords internal
format_props <- function(props, digits = 2){
  
  
  out <- ifelse(is.na(props), "", paste0(" (", formatC(as.numeric(props), format = "f", digits = digits, drop0trailing = FALSE), "%)"))
  
  out
  
  
}


#' Paste counts and proportions corresponding to one subgroup
#' 
#' @param counts Vector with counts.
#' @param props Vector with proportions.
#' @param digits Number of decimal places when rounding proportions.
#' @keywords internal
format_counts_and_props <- function(counts, props, digits = 2){
  
  
  out <- paste0(ifelse(is.na(counts), "", counts), ifelse(is.na(props), "", paste0(" (", formatC(as.numeric(props), format = "f", digits = digits, drop0trailing = FALSE), "%)")))
  
  # Remove white spaces from the beginning and the end of a string
  
  out <- stringr::str_trim(out, side = "both")
  
  out
  
  
  
}




#' Paste counts and proportions corresponding to one subgroup
#' 
#' @param counts Data frame with counts.
#' @param props Data frame with proportions.
#' @param digits Number of decimal places when rounding proportions.
#' @keywords internal
format_counts_and_props_df <- function(counts, props, digits = 2, prefix_counts = "counts_"){
  
  ### Match pattern at the beginning
  pattern <- paste0("^", prefix_counts)
  
  output_names <- gsub(pattern, "", colnames(counts))
  
  stopifnot(all(dim(counts) == dim(props)))
  
  
  output <- lapply(1:nrow(counts), function(i){
    # i = 1
    
    out <- format_counts_and_props(counts = counts[i, ], props = props[i, ], digits = digits)
    
    return(out)
    
  })
  
  output <- data.frame(do.call("rbind", output), stringsAsFactors = FALSE)
  
  colnames(output) <- output_names
  
  return(output)
  
}





#' Format summary statistics such as N, mean, median, min, max
#' 
#' @param summ Vector with values to format.
#' @param digits Vector with digits used in formatC.
#' 
#' @examples 
#' 
#' \dontrun{
#' summ <- c(10, 23.6442)
#' digits <- c(0, 2)
#' 
#' format_summ(summ, digits)
#' }
#' 
#' @keywords internal
format_summ <- function(summ, digits = 2){
  
  summ <- as.numeric(summ)
  
  if(length(digits) == 1){
    digits <- rep(digits, times = length(summ))
  }else{
    stopifnot(length(digits) == length(summ))
  }
  
  out <- sapply(seq_along(summ), function(i){
    
    out <- ifelse(is.na(summ[i]), "", formatC(summ[i], format = "f", digits = digits[i], drop0trailing = FALSE))
    
  })
  
  
  return(out)
  
}




#' Format summary statistics such as N, mean, median, min, max
#' 
#' @param summ Data frame.
#' @param per Whether to lapply per row or per column.
#' @param digits Vector with digits used in formatC.
#' @keywords internal
format_summ_df <- function(summ, per = "row", digits = 2){
  
  stopifnot(is.data.frame(summ) || is.matrix(summ))
  
  stopifnot(per %in% c("row", "col"))
  
  
  if(per == "row"){
    
    output <- lapply(1:nrow(summ), function(i){
      # i = 1
      
      out <- format_summ(summ[i, ], digits = digits)
      
    })
    
    output <- data.frame(do.call("rbind", output), stringsAsFactors = FALSE)
    
  }else{
    
    output <- lapply(1:ncol(summ), function(i){
      # i = 1
      
      out <- format_summ(summ[, i], digits = digits)
      
    })
    
    output <- data.frame(do.call("cbind", output), stringsAsFactors = FALSE)
    
  }
  
  colnames(output) <- colnames(summ)
  
  return(output)
  
}



#' Create variable names
#' 
#' Make sure that variable names exists for all the variables in the data frame. Creates from scratch or adds missing variable names for variables in data.
#' 
#' @param data Data frame.
#' @param variable_names Named vector of variable names corresponding to variables in data. This vector does not have to contain names for all the variables in data. If names for some variables are missing, they will be created. If NULL, variable names are created by subtracting underscore from the column names of data.
#' @return Named vector of (optionally unique) variable names for all variables from data.
#' @export
format_variable_names <- function(data, variable_names = NULL, unique = FALSE){
  
  
  new_variable_names <- gsub("_", " ", colnames(data))
  names(new_variable_names) <- colnames(data)
  
  
  if(!is.null(variable_names)){
    
    ## If variable_names is provided, it has to be a named vector
    stopifnot(!is.null(names(variable_names)))
    
    ## If there are multiple names corresponding to one variable, we keep the first one 
    if(sum(duplicated(names(variable_names))) > 0){
      warning("\nSome variables have multiple names assigned. The first one will be kept. Please double check if this is the name that you want to use.\n")
    }
    
    variable_names <- variable_names[!duplicated(names(variable_names))]
    
    mm <- match(names(new_variable_names), names(variable_names))
    new_variable_names[!is.na(mm)] <- variable_names[stats::na.omit(mm)]
    
  }
  
  
  if(unique){
    
    ## There cannot be duplicated names. If so, fix it.
    
    if(sum(duplicated(new_variable_names)) > 0){
      
      dupl_names <- unique(new_variable_names[duplicated(new_variable_names)])
      
      fixed_names <- lapply(seq_along(dupl_names), function(i){
        # i = 1
        
        dupl_indx <- which(new_variable_names == dupl_names[i])
        
        new_names <- paste(dupl_names[i], 1:length(dupl_indx))
        
        out <- data.frame(dupl_indx = dupl_indx, new_names = new_names, stringsAsFactors = FALSE)
        
        
      })
      
      fixed_names <- plyr::rbind.fill(fixed_names)
      
      new_variable_names[fixed_names[, "dupl_indx"]] <- fixed_names[, "new_names"]
      
    }
    
  }
  
  
  
  return(new_variable_names)
  
  
}



#' Format or create colors for a vector with categorical values
#' 
#' Make sure that unique named colors are created for each value.
#' 
#' @param x Vector of categorical values for which we want to specify colors.
#' @param colors Vector of colors longer or equal the number of unique levels of x. Can be named or non-named. If NULL, colors are generated.
#' @param palette Vector of at least two colors used to create a color palette with 'colorRampPalette' or name of a RColorBrewer palette e.g. "Oranges", "Spectral".
#' @return Named vector of unique colors for all unique values of x.
#' 
#' @examples 
#' 
#' x <- c("low", "high")
#' 
#' colors <- c("high" = "red", "low" = "grey")
#' palette <- NULL
#' allow_duplicated <- FALSE
#' 
#' 
#' format_colors_cat(x, colors = colors, palette = palette, allow_duplicated = allow_duplicated)
#' 
#' 
#' x <- c("<1", "<2", "<3", "<4")
#'
#' colors <- NULL
#' palette <- "RdYlBu"
#' allow_duplicated <- FALSE
#' 
#' 
#' out <- format_colors_cat(x, colors = colors, palette = palette, allow_duplicated = allow_duplicated)
#' 
#' barplot(rep(1, length(out)), col = out)
#' 
#' 
#' @export
format_colors_cat <- function(x, colors = NULL, palette = NULL, rev = FALSE, allow_duplicated = TRUE){
  
  
  x <- x[!is.na(x)]
  
  if(!is.factor(x)){
    x <- factor(x, levels = unique(x))
  }
  
  levels_x <- levels(x)
  
  
  
  if(is.null(colors)){
    
    if(is.null(palette)){
      
      ### d3 20 - light color first
      
      # colors_default <- c("#aec7e8", "#1f77b4",  "#ffbb78", "#ff7f0e", "#98df8a", "#2ca02c", "#ff9896", "#d62728", "#c5b0d5", "#9467bd", "#c49c94", "#8c564b", "#f7b6d2", "#e377c2", "#c7c7c7", "#7f7f7f", "#dbdb8d", "#bcbd22", "#9edae5", "#17becf")
      
      # barplot(rep(1, length(colors_default)), col = colors_default)
      
      
      
      ### paired
      
      # colors_default <- c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928")
      
      # barplot(rep(1, length(colors_default)), col = colors_default)
      
      
      
      ### Mix d3 20 - light color first with paired colors from brewer.pal
      
      colors_default <- c("#aec7e8", "#1F78B4", "#FDBF6F", "#FF7F00", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#c5b0d5", "#9467bd", "#c49c94", "#8c564b", "#f7b6d2", "#e377c2", "#c7c7c7", "#7f7f7f", "#dbdb8d", "#bcbd22", "#9edae5", "#17becf")
      
      # barplot(rep(1, length(colors_default)), col = colors_default)
      
      
      stopifnot(length(levels_x) <= 20)
      
      out <- colors_default[1:length(levels_x)]
      names(out) <- levels_x
      
      
    }else{
      
      if(length(palette) == 1){
        
        
        if(palette %in% c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")){
          
          n <- 9
          
          out <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, palette)[-c(1, n)])(length(levels_x))
          
        }else{
          
          n <- 11  
          
          out <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, palette)[-c(1, 5, 6, 7, n)])(length(levels_x))
          
        }
        
        
        
        
        if(rev){
          out <- rev(out)
        }
        
        names(out) <- levels_x
        
        
        # barplot(rep(1, length(out)), col = out)
        
        
      }else{
        out <- grDevices::colorRampPalette(palette)(length(levels_x))
        names(out) <- levels_x
      }
      
      
    }
    
    
    # barplot(rep(1, length(out)), col = out)
    
    
  }else{
    
    stopifnot(length(colors) >= length(levels_x))
    
    if(is.null(names(colors))){
      out <- colors[1:length(levels_x)]
      names(out) <- levels_x
    }else{
      stopifnot(all(levels_x %in% names(colors)))
      out <- colors[levels_x]
    }
    
    ### Colors have to be unique for ggsurvplot. Otherwise, it does not work.
    if(!allow_duplicated){
      stopifnot(sum(duplicated(out)) == 0)
    }
    
    
  }
  
  
  return(out)
  
  
}




#' @rdname format_colors_cat
#' @export
format_colors <- function(x, colors = NULL, palette = NULL, rev = FALSE, allow_duplicated = TRUE){
  
  format_colors_cat(x, colors = colors, palette = palette, rev = rev, allow_duplicated = allow_duplicated)
  
}






#' @rdname format_colors_cat
#' @param strata Vector of categorical values of stratification groups.
#' @param palette List with palettes for the different strata.
#' @return Named vector of unique colors for all unique values of x.
#' @export
format_colors_cat_strata <- function(x, strata = NULL, palette = NULL, rev = FALSE){
  
  
  x <- x[!is.na(x)]
  
  if(!is.factor(x)){
    x <- factor(x, levels = unique(x))
  }
  
  levels_x <- levels(x)
  
  
  
  if(is.null(strata)){
    strata <- "strata_dummy"
  }
  
  strata <- strata[!is.na(strata)]
  
  if(!is.factor(strata)){
    strata <- factor(strata, levels = unique(strata))
  }
  
  levels_strata <- levels(strata)
  
  
  if(is.null(palette)){
    
    palette <- lapply(seq_along(levels_strata), function(i){
      c("Oranges", "Blues", "Greens", "Purples", "Reds", "Greys")[i]
    })
    
  }else{
    
    stopifnot(length(palette) == length(levels_strata))
    
  }
  
  
  
  out <- lapply(seq_along(levels_strata), function(i){
    # i = 1 
    
    x <- paste0(levels_strata[i], ", ", levels_x)
    
    
    out <- format_colors_cat(x, palette = palette[[i]], rev = rev, allow_duplicated = FALSE)
    
    
  })
  
  
  out <- unlist(out, use.names = TRUE)
  
  names(out) <- gsub("^strata_dummy, ", "", names(out))
  
  
  out
  
  
  
}





compute_upper_whisker <- function(x, range = 1.5){
  
  
  first_quartile <- quantile(x, probs = 0.25, na.rm = TRUE)
  third_quartile <- quantile(x, probs = 0.75, na.rm = TRUE)
  
  IQR <- third_quartile - first_quartile
  
  upper_whisker <- as.numeric(third_quartile + range * IQR)
  
  
}



compute_lower_whisker <- function(x, range = 1.5){
  
  
  first_quartile <- quantile(x, probs = 0.25, na.rm = TRUE)
  third_quartile <- quantile(x, probs = 0.75, na.rm = TRUE)
  
  IQR <- third_quartile - first_quartile
  
  lower_whisker <- as.numeric(first_quartile - range * IQR)
  
  
}






#' Generate colors for ComplexHeatmap for numerical variables 
#' 
#' @param x Vector of numerical values for which we want to specify colors.
#' @param palette Vector of at least two colors used to create a color palette with 'colorRampPalette' or name of a RColorBrewer palette e.g. "Oranges", "Spectral".
#' @return Vector of unique colors for all unique values of x.
#' 
#' @examples 
#' 
#' x <- rnorm(20)
#' 
#' palette <- "Spectral"
#' 
#' format_colors_num(x, trim_values = 2.5)
#' 
#' 
#' @export
format_colors_num <- function(x, centered = TRUE, palette = NULL, rev = FALSE, trim_values = NULL, trim_prop = NULL, trim_range = NULL){
  
  
  x <- c(x)
  x <- x[!is.na(x)]
  
  if(is.null(palette)){
    if(centered){
      palette <- c("dodgerblue1", "dodgerblue3", "dodgerblue4", "black", "firebrick4", "firebrick3", "firebrick1")
      # barplot(rep(1, length(palette)), col = palette)
    }else{
      palette <- rev(grDevices::hcl.colors(10, palette = "viridis"))
    }
  }
  
  
  if(length(palette) == 1){
    
    if(palette %in% c("Blues", "BuGn", "BuPu", "GnBu", "Greens", "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd")){
      
      n <- 9
      
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, palette)[-c(1, n)])(n)
      
    }else{
      
      n <- 11  
      
      colors <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n, palette)[-c(1, n)])(n)
      
    }
    
    if(rev){
      colors <- rev(colors)
    }
    
    # barplot(rep(1, length(colors)), col = colors)
    
    
  }else{
    
    colors <- palette
    
  }
  
  
  if(centered){
    
    
    if(is.null(trim_values) && is.null(trim_prop) && is.null(trim_range)){
      max_abs_value <- max(abs(range(x, na.rm = TRUE)))
    }
    
    
    if(!is.null(trim_range)){
      ### Use whiskers
      
      lower_whisker <- compute_lower_whisker(x, range = trim_range)
      upper_whisker <- compute_upper_whisker(x, range = trim_range)
      
      max_abs_value <- max(abs(c(lower_whisker, upper_whisker)))
    }
    
    
    if(!is.null(trim_prop)){
      ### Use quantiles 
      max_abs_value <- max(abs(quantile(x, probs = c(trim_prop, 1 - trim_prop), na.rm = TRUE)))
    }
    
    
    if(!is.null(trim_values)){
      max_abs_value <- max(trim_values)
    }
    
    
    breaks <- seq(-max_abs_value, max_abs_value, length.out = length(colors))
    
    out <- circlize::colorRamp2(breaks, colors)
    
    
  }else{
    
    
    
    if(is.null(trim_values) && is.null(trim_prop) && is.null(trim_range)){
      range_value <- range(x, na.rm = TRUE)
    }
    
    
    if(!is.null(trim_range)){
      ### Use whiskers
      
      lower_whisker <- compute_lower_whisker(x, range = trim_range)
      upper_whisker <- compute_upper_whisker(x, range = trim_range)
      
      range_value <- c(lower_whisker, upper_whisker)
    }
    
    
    if(!is.null(trim_prop)){
      ### Use quantiles 
      range_value <- quantile(x, probs = c(trim_prop, 1 - trim_prop), na.rm = TRUE)
    }
    
    
    if(!is.null(trim_values)){
      range_value <- trim_values
    }
    
    
    breaks <- seq(range_value[1], range_value[2], length.out = length(colors))
    
    out <- circlize::colorRamp2(breaks, colors)
    
    
    
  }
  
  
  return(out)
  
  
}























#' Format or create colors for factor and numerical variables in a data frame
#' 
#' 
#' @param data Data frame.
#' @param color_list Optional color list with color definitions for some of the variables.
#' @param colors Vector of colors longer or equal the number of levels. Can be named or non-named. If NULL, colors are created.
#' @param palette_cat Vector of at least two colors used to create a color palette with 'colorRampPalette' or name of a RColorBrewer palette with 9 colors e.g. "Oranges" for categorical variables.
#' @param palette_num Vector of at least two colors used to create a color palette with 'colorRampPalette' or name of a RColorBrewer palette with 9 colors e.g. "Oranges" for numerical variables.
#' @return Named list with colors for factor and numerical variables in a data frame.
#' @export
format_color_list <- function(data, color_list = NULL, colors = NULL, palette_cat = NULL, palette_num = NULL, allow_duplicated = TRUE){
  
  
  
  
  
  
  
}















#' Calculate break time used in KM plots
#' 
#' @param x Vector with time-to-event data.
#' @param n_breaks Number of breaks.
#' @examples 
#' \dontrun{
#' data(bdata)
#' 
#' calculate_break_time(bdata$PFS)
#' }
#' @keywords internal
calculate_break_time <- function(x, n_breaks = 10){
  
  break_time_by <- max(x, na.rm = TRUE) %/% n_breaks
  
  break_time_by
  
}







#' Density Values for Smooth Density Plots
#' 
#' It produces a vector containing density values which encode the local densities at each point in a scatterplot.
#' 
#' @details It is based on function densCols from grDevices.
#' @return Vector with density level.
#' @export
densVals <- function(x, y = NULL, nbin = 128, bandwidth){
  
  xy <- grDevices::xy.coords(x, y, setLab = FALSE)
  
  select <- is.finite(xy$x) & is.finite(xy$y)
  
  x <- cbind(xy$x, xy$y)[select, ]
  
  map <- grDevices:::.smoothScatterCalcDensity(x, nbin, bandwidth)
  
  mkBreaks <- function(u) u - diff(range(u))/(length(u) - 1)/2
  
  xbin <- cut(x[, 1], mkBreaks(map$x1), labels = FALSE)
  
  ybin <- cut(x[, 2], mkBreaks(map$x2), labels = FALSE)
  
  dens <- map$fhat[cbind(xbin, ybin)]
  
  dens[is.na(dens)] <- 0
  
  colpal <- cut(dens, length(dens), labels = FALSE)
  
  vals <- rep(NA, length(select))
  
  vals[select] <- colpal
  
  vals
  
}










