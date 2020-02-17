



# data <- data_goya
# 
# col_var <- "Cell_Of_Origin2"
# row_var <- "FCGR2B"
# 
# variable_names = NULL
# caption = NULL
# 
# print_pvalues = TRUE
# 
# method = "kruskal"


#' Kruskal–Wallis H test or Wilcoxon Rank-Sum test
#' 
#' @param data Data frame.
wrapper_core_kruskal_test <- function(data, col_var, row_var, method = "kruskal", variable_names = NULL, caption = NULL, print_pvalues = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(method %in% c("kruskal", "wilcox"))
  
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(col_var, row_var), drop = FALSE]
  
  
  stopifnot(length(col_var) == 1)
  stopifnot(is.factor(data[, col_var]))
  stopifnot(nlevels(data[, col_var]) >= 2)
  
  if(method == "wilcox"){
    stopifnot(nlevels(data[, col_var]) == 2)
  }
  
  
  
  stopifnot(length(row_var) == 1)
  stopifnot(is.numeric(data[, row_var]) || is.integer(data[, row_var]))
  
  ### Keep non-missing data
  
  data <- data[complete.cases(data[, c(col_var, row_var)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics and do testing
  # --------------------------------------------------------------------------
  
  N = aggregate(data[, row_var], list(subgroup = data[, col_var]), FUN = length, drop = FALSE)[, 2]
  
  Median = aggregate(data[, row_var], list(subgroup = data[, col_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = aggregate(data[, row_var], list(subgroup = data[, col_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = aggregate(data[, row_var], list(subgroup = data[, col_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = aggregate(data[, row_var], list(subgroup = data[, col_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- t(data.frame(N, Median, Mean, Min, Max))
  colnames(summdf) <- levels(data[, col_var])
  
  
  
  tbl <- table(data[, col_var])
  
  if(sum(tbl > 1) >= 2){
    
    test_res <- NULL
    
    if(method == "wilcox"){
      ## Wilcoxon Rank-Sum test
      levels_col_var <- levels(data[, col_var])
      
      try(test_res <- wilcox.test(x = data[data[, col_var] == levels_col_var[1], row_var], 
        y = data[data[, col_var] == levels_col_var[2], row_var]), silent = TRUE)
      
    }else if(method == "kruskal"){
      ## Kruskal–Wallis H test
      
      try(test_res <- kruskal.test(x = data[, row_var], g = data[, col_var]), silent = TRUE)
      
    }
    
    
    if(is.null(test_res)){
      pvalue <- NA
      fc <- NA
    }else{
      pvalue <- test_res$p.value
      if(length(tbl) == 2 && sum(tbl > 1) >= 2){
        fc <- Median[2]/Median[1]
      }else{
        fc <- NA
      }
    }
    
    
  }else{
    pvalue <- NA
    fc <- NA
  }
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  res <- data.frame(covariate = row_var,
    statistic = rownames(summdf),
    summdf,
    fc = c(fc, rep(NA, nrow(summdf) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(summdf) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Statistic = res$statistic, 
    
    format_summ(summ = summdf),
    
    FC = format_or(res$fc, digits = 2),
    
    `P-value` = format_pvalues(res$pvalue), 
    
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  ### If all FC are empty, do not display that column.
  if(all(out$FC == "")){
    out$FC <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  num_start_cols <- 2
  num_end_cols <- sum(c("FC", "P-value") %in% colnames(out))
  
  
  header <- c(num_start_cols, nlevels(data[, col_var]), num_end_cols)
  header <- as.integer(header)
  names(header) <- c(" ", variable_names[col_var], " ")
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    if(method == "wilcox"){
      caption <- paste0("Wilcoxon Rank-Sum test.")
    }else if(method == "kruskal"){
      caption <- paste0("Kruskal–Wallis H test.")
    }
    
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}




















