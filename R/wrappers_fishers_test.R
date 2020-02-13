



# data <- data_goya
# 
# x_var <- "Cell_Of_Origin"
# y_var <- "Ann_Arbor_Stage"
# 
# # x_var <- "Cell_Of_Origin2"
# # y_var <- "FCGR2B_cat2"
# 
# variable_names = NULL
# caption = NULL
# 
# margin = 1
# 
# print_pvalues = TRUE
# print_adjpvalues = TRUE




#' Fisher's test
#' 
#' @param data Data frame.
wrapper_core_fishers_test <- function(data, x_var, y_var, variable_names = NULL, caption = NULL, margin = 1, print_pvalues = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(x_var) == 1)
  stopifnot(is.factor(data[, x_var]))
  stopifnot(nlevels(data[, x_var]) >= 2)
  
  stopifnot(length(y_var) == 1)
  stopifnot(is.factor(data[, y_var]))
  stopifnot(nlevels(data[, y_var]) >= 2)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions and do testing
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, y_var], data[, x_var])
  
  prop <- prop.table(tbl, margin = margin) * 100
  
  
  if(sum(margin.table(tbl, margin = 1) > 1) >= 2 && sum(margin.table(tbl, margin = 2) > 1) >= 2){
    ## Pearson’s Chi-squared test: get difference of proportions of succeses and CI for difference
    ## Success category is defined by the first column
    # test_res <- prop.test(tbl)
    
    ## Fisher's exact test: get odds rations and CI for OR
    test_res <- NULL
    try(test_res <- fisher.test(tbl), silent = TRUE)
    if(is.null(test_res)){
      test_res <- fisher.test(tbl, simulate.p.value = TRUE)
    }
    
    if(is.null(test_res)){
      pvalue <- NA
      or <- NA
    }else{
      pvalue <- test_res$p.value
      or <- test_res$estimate
      if(is.null(or)){
        or <- NA
      }
    }
    
  }else{
    pvalue <- NA
    or <- NA
  }
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  countdf <- as.data.frame.matrix(tbl)
  colnames(countdf) <- paste0("counts_", colnames(countdf))
  
  propdf <- as.data.frame.matrix(prop)
  colnames(propdf) <- paste0("proportions_", colnames(propdf))
  
  
  res <- data.frame(covariate = y_var,
    subgroup = rownames(tbl),
    countdf, propdf,
    or = c(or, rep(NA, nrow(tbl) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(tbl) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = res$subgroup, 
    
    format_counts_and_props(counts = countdf, props = propdf),
    
    OR = format_or(res$or, digits = 2),
    `P-value` = format_pvalues(res$pvalue), 
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  ### If all OR are empty, do not display that column.
  if(all(out$OR == "")){
    out$OR <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  num_start_cols <- 2
  num_end_cols <- sum(c("OR", "P-value") %in% colnames(out))
  
  
  header <- c(num_start_cols, nlevels(data[, x_var]), num_end_cols)
  header <- as.integer(header)
  names(header) <- c(" ", variable_names[x_var], " ")
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Fisher's exact test.")
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}




















