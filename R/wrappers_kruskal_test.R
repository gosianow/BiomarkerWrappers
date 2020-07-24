




#' Kruskal-Wallis H test or Wilcoxon Rank-Sum test
#' 
#' Returns a table where stratification subgroups are in columns and distribution summary statistics for the numerical variable are in rows. 
#' 
#' @param data Data frame.
#' @param num_var Name of a numerical variable.
#' @param cat_var Name of a categorical variable. That variable must be a factor with at least two levels.
#' @export
wrapper_core_kruskal_test_col_cat <- function(data, num_var, cat_var, method = "kruskal", variable_names = NULL, caption = NULL, display_statistics = c("N", "Median"), force_empty_cols = FALSE, print_pvalues = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(method %in% c("kruskal", "wilcox"))
  
  stopifnot(length(display_statistics) >= 1)
  stopifnot(display_statistics %in% c("N", "Median", "Mean", "Min", "Max"))
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  stopifnot(length(cat_var) == 1)
  stopifnot(is.factor(data[, cat_var]))
  stopifnot(nlevels(data[, cat_var]) >= 2)
  
  if(method == "wilcox"){
    stopifnot(nlevels(data[, cat_var]) == 2)
  }
  
  
  stopifnot(length(num_var) == 1)
  stopifnot(is.numeric(data[, num_var]) || is.integer(data[, num_var]))
  
  ### Keep non-missing data
  data <- data[stats::complete.cases(data[, c(cat_var, num_var)]), c(cat_var, num_var), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics and do testing
  # --------------------------------------------------------------------------
  
  N = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = length, drop = FALSE)[, 2]
  
  Median = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- data.frame(N, Median, Mean, Min, Max)
  summdf <- summdf[, display_statistics, drop = FALSE]
  summdf <- t(summdf)
  colnames(summdf) <- levels(data[, cat_var])
  
  
  
  tbl <- table(data[, cat_var])
  
  if(sum(tbl > 1) >= 2){
    
    test_res <- NULL
    
    if(method == "wilcox"){
      ## Wilcoxon Rank-Sum test
      levels_cat_var <- levels(data[, cat_var])
      
      try(test_res <- stats::wilcox.test(x = data[data[, cat_var] == levels_cat_var[1], num_var], 
        y = data[data[, cat_var] == levels_cat_var[2], num_var]), silent = TRUE)
      
    }else if(method == "kruskal"){
      ## Kruskal-Wallis H test
      
      try(test_res <- stats::kruskal.test(x = data[, num_var], g = data[, cat_var]), silent = TRUE)
      
    }
    
    
    if(is.null(test_res)){
      pvalue <- NA
      difference <- NA
    }else{
      pvalue <- test_res$p.value
      if(length(tbl) == 2 && sum(tbl > 1) >= 2){
        difference <- Median[2] - Median[1]
      }else{
        difference <- NA
      }
    }
    
    
  }else{
    pvalue <- NA
    difference <- NA
  }
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  res <- data.frame(covariate = num_var,
    statistic = rownames(summdf),
    summdf,
    difference = c(difference, rep(NA, nrow(summdf) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(summdf) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  digits <- rep(2, length(display_statistics))
  digits[display_statistics == "N"] <- 0
  
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Statistic = res$statistic, 
    format_summ(summ = summdf, digits = digits, per = "col"),
    Difference = format_difference(res$difference, digits = 2, non_empty = 1),
    `P-value` = format_pvalues(res$pvalue, non_empty = 1), 
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference %in% c("", "NA")) && !force_empty_cols){
    out$Difference <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  
  header <- format_header(all_colnames = colnames(out), header_colnames = colnames(summdf), header_name = variable_names[cat_var])
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    if(method == "wilcox"){
      caption <- paste0("Wilcoxon Rank-Sum test.")
    }else if(method == "kruskal"){
      caption <- paste0("Kruskal-Wallis H test.")
    }
    
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}



#' Kruskal-Wallis H test or Wilcoxon Rank-Sum test
#' 
#' Returns a table where stratification subgroups are in rows and distribution summary statistics for the numerical variable are in columns. 
#' 
#' @inheritParams wrapper_core_kruskal_test_col_cat
#' @export
wrapper_core_kruskal_test_col_num <- function(data, num_var, cat_var, method = "kruskal", variable_names = NULL, caption = NULL, display_statistics = c("N", "Median"), force_empty_cols = FALSE, print_pvalues = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(method %in% c("kruskal", "wilcox"))
  
  stopifnot(length(display_statistics) >= 1)
  stopifnot(display_statistics %in% c("N", "Median", "Mean", "Min", "Max"))
  
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(num_var, cat_var), drop = FALSE]
  
  
  stopifnot(length(cat_var) == 1)
  stopifnot(is.factor(data[, cat_var]))
  stopifnot(nlevels(data[, cat_var]) >= 2)
  
  if(method == "wilcox"){
    stopifnot(nlevels(data[, cat_var]) == 2)
  }
  
  
  stopifnot(length(num_var) == 1)
  stopifnot(is.numeric(data[, num_var]) || is.integer(data[, num_var]))
  
  
  ### Keep non-missing data
  
  data <- data[stats::complete.cases(data[, c(num_var, cat_var)]), ]
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics and do testing
  # --------------------------------------------------------------------------
  
  
  N = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = length, drop = FALSE)[, 2]
  
  Median = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- data.frame(N, Median, Mean, Min, Max)
  summdf <- summdf[, display_statistics, drop = FALSE]
  rownames(summdf) <- levels(data[, cat_var])
  
  
  
  tbl <- table(data[, cat_var])
  
  if(sum(tbl > 1) >= 2){
    
    test_res <- NULL
    
    if(method == "wilcox"){
      ## Wilcoxon Rank-Sum test
      levels_num_var <- levels(data[, cat_var])
      
      try(test_res <- stats::wilcox.test(x = data[data[, cat_var] == levels_num_var[1], num_var], 
        y = data[data[, cat_var] == levels_num_var[2], num_var]), silent = TRUE)
      
    }else if(method == "kruskal"){
      ## Kruskal-Wallis H test
      
      try(test_res <- stats::kruskal.test(x = data[, num_var], g = data[, cat_var]), silent = TRUE)
      
    }
    
    
    if(is.null(test_res)){
      pvalue <- NA
      difference <- NA
    }else{
      pvalue <- test_res$p.value
      if(length(tbl) == 2 && sum(tbl > 1) >= 2){
        difference <- Median[2] - Median[1]
      }else{
        difference <- NA
      }
    }
    
    
  }else{
    pvalue <- NA
    difference <- NA
  }
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  
  res <- data.frame(covariate = cat_var,
    subgroup = rownames(summdf),
    summdf,
    difference = c(difference, rep(NA, nrow(summdf) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(summdf) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  digits <- rep(2, length(display_statistics))
  digits[display_statistics == "N"] <- 0
  
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = res$subgroup, 
    format_summ(summ = summdf, per = "row", digits = digits),
    Difference = format_difference(res$difference, digits = 2, non_empty = 1),
    `P-value` = format_pvalues(res$pvalue, non_empty = 1), 
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference %in% c("", "NA")) && !force_empty_cols){
    out$Difference <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  header <- format_header(all_colnames = colnames(out), header_colnames = colnames(summdf), header_name = variable_names[num_var])
  
  

  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    if(method == "wilcox"){
      caption <- paste0("Wilcoxon Rank-Sum test.")
    }else if(method == "kruskal"){
      caption <- paste0("Kruskal-Wallis H test.")
    }
    
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}



#' @rdname wrapper_core_kruskal_test_col_cat
#' @inheritParams wrapper_core_kruskal_test_col_cat
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_core_kruskal_test_col_cat_strat <- function(data, num_var, cat_var, strat1_var = NULL, strat2_var = NULL, method = "kruskal", variable_names = NULL, caption = NULL, display_statistics = c("N", "Median"), force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check on strat vars
  # --------------------------------------------------------------------------
  
  if(!is.null(strat1_var)){
    stopifnot(length(strat1_var) == 1)
    stopifnot(is.factor(data[, strat1_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat1_dummy"] <- factor("strat1_dummy")
    strat1_var <- "strat1_dummy"
  }
  
  if(!is.null(strat2_var)){
    stopifnot(length(strat2_var) == 1)
    stopifnot(is.factor(data[, strat2_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat2_dummy"] <- factor("strat2_dummy")
    strat2_var <- "strat2_dummy"
  }
  
  ### Strata cannot include cat_vars
  stopifnot(length(intersect(c(strat1_var, strat2_var), cat_var)) == 0)
  
  if(!print_pvalues){
    print_adjpvalues <- FALSE
  }
  
  ### Keep non-missing data
  data <- data[stats::complete.cases(data[, c(num_var, cat_var, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  # --------------------------------------------------------------------------
  # Calculations within strata
  # --------------------------------------------------------------------------
  
  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  wrapper_res <- lapply(1:length(strata2_levels), function(j){
    # j = 1
    
    data_strata2 <- data[data[, strat2_var] == strata2_levels[j] & !is.na(data[, strat2_var]), ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    wrapper_res <- lapply(1:length(strata1_levels), function(i){
      # i = 1
      
      # print(paste(j, i))
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      wrapper_res <- wrapper_core_kruskal_test_col_cat(data = data_strata1, num_var = num_var, cat_var = cat_var, method = method, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues)
      
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      hdr <- bheader(wrapper_res)
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      
      ## Update header by adding 2 corresponding to the two strat variables to the first position
      hdr[1] <- hdr[1] + 2
      
      rownames(res) <- NULL
      rownames(out) <- NULL
      
      wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res), header = hdr)
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = bheader(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  hdr <- bheader(wrapper_res[[1]])
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if(print_adjpvalues){
    out$`Adj. P-value` <- format_pvalues(stats::p.adjust(res$pvalue, method = "BH"))
    ## Update header
    hdr[3] <- hdr[3] + 1
  }
  

  
  ### Set repeating Strata names to empty
  out[out$Covariate == "", variable_names[strat1_var]] <- ""
  out[out$Covariate == "", variable_names[strat2_var]] <- ""
  
  ### Remove dummy columns
  
  hdr_shift <- 0
  
  if(strat2_var == "strat2_dummy"){
    res$strat2_dummy <- NULL
    out$`strat2 dummy` <- NULL
    hdr_shift <- hdr_shift + 1
  }
  if(strat1_var == "strat1_dummy"){
    res$strat1_dummy <- NULL
    out$`strat1 dummy` <- NULL
    hdr_shift <- hdr_shift + 1
  }
  
  ## Update header 
  hdr[1] <- hdr[1] - hdr_shift
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
  
}



#' @rdname wrapper_core_kruskal_test_col_num
#' @inheritParams wrapper_core_kruskal_test_col_num
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_core_kruskal_test_col_num_strat <- function(data, num_var, cat_var, strat1_var = NULL, strat2_var = NULL, method = "kruskal", variable_names = NULL, caption = NULL, display_statistics = c("N", "Median"), force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Check on strat vars
  # --------------------------------------------------------------------------
  
  if(!is.null(strat1_var)){
    stopifnot(length(strat1_var) == 1)
    stopifnot(is.factor(data[, strat1_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat1_dummy"] <- factor("strat1_dummy")
    strat1_var <- "strat1_dummy"
  }
  
  if(!is.null(strat2_var)){
    stopifnot(length(strat2_var) == 1)
    stopifnot(is.factor(data[, strat2_var]))
  }else{
    ### Add dummy variable to data
    data[, "strat2_dummy"] <- factor("strat2_dummy")
    strat2_var <- "strat2_dummy"
  }
  
  ### Strata cannot include cat_vars
  stopifnot(length(intersect(c(strat1_var, strat2_var), cat_var)) == 0)
  
  if(!print_pvalues){
    print_adjpvalues <- FALSE
  }
  
  ### Keep non-missing data
  data <- data[stats::complete.cases(data[, c(num_var, cat_var, strat1_var, strat2_var)]), ]
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  # --------------------------------------------------------------------------
  # Calculations within strata
  # --------------------------------------------------------------------------
  
  strata1_levels <- levels(data[, strat1_var])
  strata2_levels <- levels(data[, strat2_var])
  
  
  wrapper_res <- lapply(1:length(strata2_levels), function(j){
    # j = 1
    
    data_strata2 <- data[data[, strat2_var] == strata2_levels[j] & !is.na(data[, strat2_var]), ]
    
    if(nrow(data_strata2) == 0){
      return(NULL)
    }
    
    
    wrapper_res <- lapply(1:length(strata1_levels), function(i){
      # i = 1
      
      # print(paste(j, i))
      
      data_strata1 <- data_strata2[data_strata2[, strat1_var] == strata1_levels[i] & !is.na(data_strata2[, strat1_var]), ]
      
      if(nrow(data_strata1) == 0){
        return(NULL)
      }
      
      
      wrapper_res <- wrapper_core_kruskal_test_col_num(data = data_strata1, num_var = num_var, cat_var = cat_var, method = method, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues)
      
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      hdr <- bheader(wrapper_res)
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      
      ## Update header by adding 2 corresponding to the two strat variables to the first position
      hdr[1] <- hdr[1] + 2
      
      rownames(res) <- NULL
      rownames(out) <- NULL
      
      wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res), header = hdr)
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = bheader(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  hdr <- bheader(wrapper_res[[1]])
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if(print_adjpvalues){
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
    ## Update header
    hdr[3] <- hdr[3] + 1
  }
  
  
  ### Set repeating Strata names to empty
  out[out$Covariate == "", variable_names[strat1_var]] <- ""
  out[out$Covariate == "", variable_names[strat2_var]] <- ""
  
  ### Remove dummy columns
  
  hdr_shift <- 0
  
  if(strat2_var == "strat2_dummy"){
    res$strat2_dummy <- NULL
    out$`strat2 dummy` <- NULL
    hdr_shift <- hdr_shift + 1
  }
  if(strat1_var == "strat1_dummy"){
    res$strat1_dummy <- NULL
    out$`strat1 dummy` <- NULL
    hdr_shift <- hdr_shift + 1
  }
  
  ## Update header
  hdr[1] <- hdr[1] - hdr_shift
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
  
}






#' Kruskal-Wallis H test or Wilcoxon Rank-Sum test
#' 
#' Run Kruskal-Wallis H test or Wilcoxon Rank-Sum test for multiple covariates
#' 
#' @inheritParams wrapper_core_kruskal_test_col_cat_strat
#' @param num_vars Vector with names of numerical variables. If it has length >= 1, then 'cat_var' must be of length 1, and stratification subgroups are displayed in columns and statistics in rows.
#' @param cat_vars Vector with names of categorical variables. If it has length >= 1, then 'num_var' must be of length 1, and stratification subgroups are displayed in rows and statistics in columns.
#' @export
wrapper_kruskal_test <- function(data, num_vars, cat_vars, strat1_var = NULL, strat2_var = NULL, method = "kruskal", variable_names = NULL, caption = NULL, display_statistics = c("N", "Median"), display_in_column = "cat", force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  stopifnot(length(num_vars) >= 1)
  stopifnot(length(cat_vars) >= 1)
  
  
  if(length(num_vars) > 1){
    stopifnot(length(cat_vars) == 1)
    display_in_column <- "cat"
  }
  
  if(length(cat_vars) > 1){
    stopifnot(length(num_vars) == 1)
    display_in_column <- "num"
  }
  
  
  stopifnot(length(display_in_column) == 1)
  stopifnot(display_in_column %in% c("cat", "num"))
  
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  
  if(display_in_column == "num"){
    
    stopifnot(length(num_vars) == 1)
    
    wrapper_res <- lapply(1:length(cat_vars), function(i){
      # i = 1
      
      # print(i)
      
      num_var <- num_vars
      cat_var <- cat_vars[i]
      
      wrapper_res <- wrapper_core_kruskal_test_col_num_strat(data = data, num_var = num_var, cat_var = cat_var, strat1_var = strat1_var, strat2_var = strat2_var, method = method, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      return(wrapper_res)
      
    })
    
    
  }else{
    
    
    wrapper_res <- lapply(1:length(num_vars), function(i){
      # i = 1
      
      num_var <- num_vars[i]
      cat_var <- cat_vars
      
      wrapper_res <- wrapper_core_kruskal_test_col_cat_strat(data = data, num_var = num_var, cat_var = cat_var, strat1_var = strat1_var, strat2_var = strat2_var, method = method, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
      
      return(wrapper_res)
      
    })
    
  }
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  hdr <- bheader(wrapper_res[[1]])
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(stats::p.adjust(res$pvalue, method = "BH"))
  }
  
 
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference %in% c("", "NA")) && !force_empty_cols){
    out$Difference <- NULL
    ### Update header
    if(hdr[3] == 1){
      hdr <- hdr[-3]
    }else{
      hdr[3] <- hdr[3] - 1
    }
  }
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
}






































