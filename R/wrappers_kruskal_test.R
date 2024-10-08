




#' Kruskal-Wallis H test or Wilcoxon Rank-Sum test (know also as Wilcoxon-Mann-Whitney test) or t-test
#' 
#' Returns a table where stratification subgroups are in columns and distribution summary statistics for the numerical variable are in rows. 
#' 
#' @param data Data frame.
#' @param num_var Name of a numerical variable.
#' @param cat_var Name of a categorical variable. That variable must be a factor with at least two levels.
#' @param method Test to be used. Possible values: "kruskal", "wilcox", "t".
#' @param display_statistics Vector of possible values: "N", "Median", "Mean", "Min", "Max", "First.Quartile", "Third.Quartile".
#' @param paired Logical. Paired test is possible only for wilcox and t. The data must be ordered by the pairing variable, usually subject ID. 
#' @export
wrapper_kruskal_test_core_col_cat <- function(data, num_var, cat_var, method = "kruskal", alternative = "two.sided", paired = FALSE, variable_names = NULL, caption = NULL, display_statistics = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, drop = FALSE){
  
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(method %in% c("kruskal", "wilcox", "t"))
  
  ### Paired test is possible only for wilcox and t
  if(paired){
    if(method == "kruskal"){
      paired <- FALSE
    }
  }
  
  if(is.null(display_statistics)){
    if(method == "t"){
      display_statistics <- c("N", "Mean")
    }else{
      display_statistics <- c("N", "Median")
    }
    
  }else{
    stopifnot(length(display_statistics) >= 1)
    stopifnot(display_statistics %in% c("N", "Median", "Mean", "Min", "Max", "First.Quartile", "Third.Quartile"))
  }
  

  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  stopifnot(length(cat_var) == 1)
  stopifnot(is.factor(data[, cat_var]))
  stopifnot(nlevels(data[, cat_var]) >= 2)
  
  if(method %in% c("wilcox", "t")){
    stopifnot(nlevels(data[, cat_var]) == 2)
  }
  
  
  stopifnot(length(num_var) == 1)
  stopifnot(is.numeric(data[, num_var]) || is.integer(data[, num_var]))
  
  ### Keep non-missing data
  if(!paired){
    data <- data[stats::complete.cases(data[, c(cat_var, num_var)]), , drop = FALSE]
  }
  
  ### Drop unused levels except the reference level
  if(drop){
    level_ref <- levels(data[, cat_var])[1]
    data[, cat_var] <- factor(data[, cat_var])
    data[, cat_var] <- factor(data[, cat_var], levels = union(level_ref, levels(data[, cat_var])))
  }
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics and do testing
  # --------------------------------------------------------------------------
  
  # N = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = length, drop = FALSE)[, 2]
  
  N <- as.numeric(table(data[stats::complete.cases(data[, c(cat_var, num_var)]), cat_var]))
  
  Median = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  First.Quartile = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = quantile, probs = 0.25, na.rm = TRUE, drop = FALSE)[, 2]
  Third.Quartile = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = quantile, probs = 0.75, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- data.frame(N, Median, Mean, Min, Max, First.Quartile, Third.Quartile)
  summdf <- summdf[, display_statistics, drop = FALSE]
  summdf <- t(summdf)
  colnames(summdf) <- levels(data[, cat_var])
  
  
  
  tbl <- table(data[, cat_var])
  
  if(sum(tbl > 1) >= 2){
    
    test_res <- NULL
    
    if(method == "wilcox"){
      ## Wilcoxon Rank-Sum test
      levels_cat_var <- levels(data[, cat_var])
      
      try(test_res <- stats::wilcox.test(x = data[data[, cat_var] == levels_cat_var[2], num_var], 
        y = data[data[, cat_var] == levels_cat_var[1], num_var], alternative = alternative, paired = paired), silent = TRUE)
      
    }else if(method == "kruskal"){
      ## Kruskal-Wallis H test
      
      try(test_res <- stats::kruskal.test(x = data[, num_var], g = data[, cat_var]), silent = TRUE)
      
    }else if(method == "t"){
      ## t-test
      levels_cat_var <- levels(data[, cat_var])
      
      try(test_res <- stats::t.test(x = data[data[, cat_var] == levels_cat_var[2], num_var], 
        y = data[data[, cat_var] == levels_cat_var[1], num_var], alternative = alternative, paired = paired), silent = TRUE)
      
    }
    
    
    if(is.null(test_res)){
      pvalue <- NA
      difference <- NA
    }else{
      pvalue <- test_res$p.value
      if(length(tbl) == 2 && sum(tbl > 1) >= 2){
        if(method == "t"){
          difference <- Mean[2] - Mean[1]
        }else{
          difference <- Median[2] - Median[1]
        }
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
    format_summ_df(summ = summdf, digits = digits, per = "col"),
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
    }else if(method == "t"){
      caption <- paste0("t-test.")
    }
    
    if(paired){
      caption <- paste0("Paired ", caption)
    }
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}



#' Kruskal-Wallis H test or Wilcoxon Rank-Sum test (know also as Wilcoxon-Mann-Whitney test) or t-test
#' 
#' Returns a table where stratification subgroups are in rows and distribution summary statistics for the numerical variable are in columns. 
#' 
#' @inheritParams wrapper_kruskal_test_core_col_cat
#' @export
wrapper_kruskal_test_core_col_num <- function(data, num_var, cat_var, method = "kruskal", alternative = "two.sided", paired = FALSE, pairwise = FALSE, variable_names = NULL, caption = NULL, display_statistics = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, drop = FALSE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(method %in% c("kruskal", "wilcox", "t"))
  
  ### Paired test is possible only for wilcox and t
  if(paired){
    if(method == "kruskal"){
      paired <- FALSE
    }
  }
  
  if(is.null(display_statistics)){
    if(method == "t"){
      display_statistics <- c("N", "Mean")
    }else{
      display_statistics <- c("N", "Median")
    }
    
  }else{
    stopifnot(length(display_statistics) >= 1)
    stopifnot(display_statistics %in% c("N", "Median", "Mean", "Min", "Max", "First.Quartile", "Third.Quartile"))
  }
  
  stopifnot(is.data.frame(data))
  stopifnot(nrow(data) > 0)
  
  ### Keep only those variables that are used for the analysis
  data <- data[, c(num_var, cat_var), drop = FALSE]
  
  
  stopifnot(length(cat_var) == 1)
  stopifnot(is.factor(data[, cat_var]))
  stopifnot(nlevels(data[, cat_var]) >= 2)
  
  if(method %in% c("wilcox", "t") && !pairwise){
    stopifnot(nlevels(data[, cat_var]) == 2)
  }
  
  
  stopifnot(length(num_var) == 1)
  stopifnot(is.numeric(data[, num_var]) || is.integer(data[, num_var]))
  
  
  ### Keep non-missing data
  if(!paired){
    data <- data[stats::complete.cases(data[, c(cat_var, num_var)]), , drop = FALSE]
  }
  
  
  ### Drop unused levels except the reference level
  if(drop){
    level_ref <- levels(data[, cat_var])[1]
    data[, cat_var] <- factor(data[, cat_var])
    data[, cat_var] <- factor(data[, cat_var], levels = union(level_ref, levels(data[, cat_var])))
  }
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Calculate summary statistics and do testing
  # --------------------------------------------------------------------------
  
  
  # N = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = length, drop = FALSE)[, 2]
  
  N <- as.numeric(table(data[stats::complete.cases(data[, c(cat_var, num_var)]), cat_var]))
  
  Median = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = median, na.rm = TRUE, drop = FALSE)[, 2]
  Mean = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = mean, na.rm = TRUE, drop = FALSE)[, 2]
  
  Min = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = min, na.rm = TRUE, drop = FALSE)[, 2]
  Max = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = max, na.rm = TRUE, drop = FALSE)[, 2]
  
  First.Quartile = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = quantile, probs = 0.25, na.rm = TRUE, drop = FALSE)[, 2]
  Third.Quartile = stats::aggregate(data[, num_var], list(subgroup = data[, cat_var]), FUN = quantile, probs = 0.75, na.rm = TRUE, drop = FALSE)[, 2]
  
  
  summdf <- data.frame(N, Median, Mean, Min, Max, First.Quartile, Third.Quartile)
  summdf <- summdf[, display_statistics, drop = FALSE]
  rownames(summdf) <- levels(data[, cat_var])
  
  
  if(pairwise){
    
    levels_cat <- levels(data[, cat_var])
    
    test_stats <- lapply(seq_along(levels_cat)[-1], function(i){
      # i = 2
      
      data_sub <- data[data[, cat_var] %in% c(levels_cat[1], levels_cat[i]), , drop = FALSE]
      
      data_sub[, cat_var] <- factor(data_sub[, cat_var], levels = c(levels_cat[1], levels_cat[i]))
      
      tbl <- table(data_sub[, cat_var])
      
      if(sum(tbl > 1) >= 2){
        
        test_res <- NULL
        levels_cat_var <- levels(data_sub[, cat_var])
        
        
        if(method == "wilcox"){
          ## Wilcoxon Rank-Sum test
          
          try(test_res <- stats::wilcox.test(x = data_sub[data_sub[, cat_var] == levels_cat_var[2], num_var], 
            y = data_sub[data_sub[, cat_var] == levels_cat_var[1], num_var], alternative = alternative, paired = paired), silent = TRUE)
          
        }else if(method == "kruskal"){
          ## Kruskal-Wallis H test
          
          try(test_res <- stats::kruskal.test(x = data_sub[, num_var], g = data_sub[, cat_var]), silent = TRUE)
          
        }else if(method == "t"){
          ## t-test
          
          try(test_res <- stats::t.test(x = data_sub[data_sub[, cat_var] == levels_cat_var[2], num_var], 
            y = data_sub[data_sub[, cat_var] == levels_cat_var[1], num_var], alternative = alternative, paired = paired), silent = TRUE)
          
        }
        
        
        if(is.null(test_res)){
          pvalue <- NA
          difference <- NA
        }else{
          pvalue <- test_res$p.value
          if(length(tbl) == 2 && sum(tbl > 1) >= 2){
            if(method == "t"){
              difference <- mean(data_sub[data_sub[, cat_var] == levels_cat_var[2], num_var], na.rm = TRUE) - mean(data_sub[data_sub[, cat_var] == levels_cat_var[1], num_var], na.rm = TRUE)
            }else{
              difference <- median(data_sub[data_sub[, cat_var] == levels_cat_var[2], num_var], na.rm = TRUE) - median(data_sub[data_sub[, cat_var] == levels_cat_var[1], num_var], na.rm = TRUE)
            }
          }else{
            difference <- NA
          }
        }
        
        
      }else{
        pvalue <- NA
        difference <- NA
      }
      
      
      return(data.frame(pvalue = pvalue, difference = difference))
      
    })
    
    test_stats <- rbind.fill(test_stats)
    
    pvalue <- test_stats$pvalue
    difference <- test_stats$difference
    
    
  }else{
      
    
    tbl <- table(data[, cat_var])
    
    if(sum(tbl > 1) >= 2){
      
      test_res <- NULL
      
      if(method == "wilcox"){
        ## Wilcoxon Rank-Sum test
        levels_cat_var <- levels(data[, cat_var])
        
        try(test_res <- stats::wilcox.test(x = data[data[, cat_var] == levels_cat_var[2], num_var], 
          y = data[data[, cat_var] == levels_cat_var[1], num_var], alternative = alternative, paired = paired), silent = TRUE)
        
      }else if(method == "kruskal"){
        ## Kruskal-Wallis H test
        
        try(test_res <- stats::kruskal.test(x = data[, num_var], g = data[, cat_var]), silent = TRUE)
        
      }else if(method == "t"){
        ## t-test
        levels_cat_var <- levels(data[, cat_var])
        
        try(test_res <- stats::t.test(x = data[data[, cat_var] == levels_cat_var[2], num_var], 
          y = data[data[, cat_var] == levels_cat_var[1], num_var], alternative = alternative, paired = paired), silent = TRUE)
        
      }
      
      
      if(is.null(test_res)){
        pvalue <- NA
        difference <- NA
      }else{
        pvalue <- test_res$p.value
        if(length(tbl) == 2 && sum(tbl > 1) >= 2){
          if(method == "t"){
            difference <- Mean[2] - Mean[1]
          }else{
            difference <- Median[2] - Median[1]
          }
        }else{
          difference <- NA
        }
      }
      
      
    }else{
      pvalue <- NA
      difference <- NA
    }
    
    
    
  }
  
  
  
 
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  if(pairwise){
    
    res <- data.frame(covariate = cat_var,
      subgroup = rownames(summdf),
      summdf,
      difference = c(NA, difference),
      pvalue = c(NA, pvalue), 
      stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
    
  }else{
    
    res <- data.frame(covariate = cat_var,
      subgroup = rownames(summdf),
      summdf,
      difference = c(difference, rep(NA, nrow(summdf) - 1)),
      pvalue = c(pvalue, rep(NA, nrow(summdf) - 1)), 
      stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
    
    
  }
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  digits <- rep(2, length(display_statistics))
  digits[display_statistics == "N"] <- 0
  
  if(pairwise){
    
    if(nrow(res) >= 2){
      non_empty <- 2:nrow(res)
    }else{
      non_empty <- NULL
    }
    
    out <- data.frame(Covariate = variable_names[res$covariate], 
      Subgroup = res$subgroup, 
      format_summ_df(summ = summdf, per = "row", digits = digits),
      Difference = format_difference(res$difference, digits = 2, non_empty = non_empty),
      `P-value` = format_pvalues(res$pvalue, non_empty = non_empty), 
      check.names = FALSE, stringsAsFactors = FALSE)
    
    
  }else{
    
    out <- data.frame(Covariate = variable_names[res$covariate], 
      Subgroup = res$subgroup, 
      format_summ_df(summ = summdf, per = "row", digits = digits),
      Difference = format_difference(res$difference, digits = 2, non_empty = 1),
      `P-value` = format_pvalues(res$pvalue, non_empty = 1), 
      check.names = FALSE, stringsAsFactors = FALSE)
    
  }
  
  
  
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
    }else if(method == "t"){
      caption <- paste0("t-test.")
    }
    
    if(paired){
      caption <- paste0("Paired ", caption)
    }
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}



#' @rdname wrapper_kruskal_test_core_col_cat
#' @inheritParams wrapper_kruskal_test_core_col_cat
#' @param strat1_var Name of the first stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_kruskal_test_core_col_cat_strat <- function(data, num_var, cat_var, strat1_var = NULL, strat2_var = NULL, method = "kruskal", alternative = "two.sided", paired = FALSE, variable_names = NULL, caption = NULL, display_statistics = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE, drop = FALSE){
  
  
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
  data <- data[stats::complete.cases(data[, c(strat1_var, strat2_var)]), ]
  
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
      
      
      wrapper_res <- wrapper_kruskal_test_core_col_cat(data = data_strata1, num_var = num_var, cat_var = cat_var, method = method, alternative = alternative, paired = paired, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues, drop = drop)
      
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      
      rownames(res) <- NULL
      rownames(out) <- NULL
      
      wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res))
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if(print_adjpvalues){
    out$`Adj. P-value` <- format_pvalues(stats::p.adjust(res$pvalue, method = "BH"))
  }
  
  ### Set repeating Strata names to empty
  out[out$Covariate == "", variable_names[strat1_var]] <- ""
  out[out$Covariate == "", variable_names[strat2_var]] <- ""
  
  ### Remove dummy columns
  if(strat2_var == "strat2_dummy"){
    res$strat2_dummy <- NULL
    out$`strat2 dummy` <- NULL
  }
  if(strat1_var == "strat1_dummy"){
    res$strat1_dummy <- NULL
    out$`strat1 dummy` <- NULL
  }
  
  ## Update header 
  hdr <- format_header(all_colnames = colnames(out), header_colnames = levels(data[, cat_var]), header_name = variable_names[cat_var])
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
  
}



#' @rdname wrapper_kruskal_test_core_col_num
#' @inheritParams wrapper_kruskal_test_core_col_num
#' @param strat1_var Name of the first stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_kruskal_test_core_col_num_strat <- function(data, num_var, cat_var, strat1_var = NULL, strat2_var = NULL, method = "kruskal", alternative = "two.sided", paired = FALSE, pairwise = FALSE, variable_names = NULL, caption = NULL, display_statistics = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE, drop = FALSE){
  
  if(is.null(display_statistics)){
    if(method == "t"){
      display_statistics <- c("N", "Mean")
    }else{
      display_statistics <- c("N", "Median")
    }
    
  }else{
    stopifnot(length(display_statistics) >= 1)
    stopifnot(display_statistics %in% c("N", "Median", "Mean", "Min", "Max", "First.Quartile", "Third.Quartile"))
  }
  
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
  data <- data[stats::complete.cases(data[, c(strat1_var, strat2_var)]), ]
  
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
      
      
      wrapper_res <- wrapper_kruskal_test_core_col_num(data = data_strata1, num_var = num_var, cat_var = cat_var, method = method, alternative = alternative, paired = paired, pairwise = pairwise, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues, drop = drop)
      
      
      res <- bresults(wrapper_res)
      out <- boutput(wrapper_res)
      
      ## Add info about the strata to the data frames
      
      prefix_df <- data.frame(strata2 = rep(strata2_levels[j], nrow(res)), strata1 = rep(strata1_levels[i], nrow(res)), stringsAsFactors = FALSE)
      
      # To res
      colnames(prefix_df) <- c(strat2_var, strat1_var)
      res <- cbind(prefix_df, res)
      
      # To out
      colnames(prefix_df) <- variable_names[c(strat2_var, strat1_var)]
      out <- cbind(prefix_df, out)
      
      rownames(res) <- NULL
      rownames(out) <- NULL
      
      wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res))
      
      return(wrapper_res)
      
    })
    
    
    ### Merge the results
    
    res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
    out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
    
    rownames(res) <- NULL
    rownames(out) <- NULL
    
    wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]))
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- stats::p.adjust(res$pvalue, method = "BH")
  
  if(print_adjpvalues){
    out$`Adj. P-value` <- format_pvalues(res$adj_pvalue)
  }
  
  
  ### Set repeating Strata names to empty
  out[out$Covariate == "", variable_names[strat1_var]] <- ""
  out[out$Covariate == "", variable_names[strat2_var]] <- ""
  
  ### Remove dummy columns
  if(strat2_var == "strat2_dummy"){
    res$strat2_dummy <- NULL
    out$`strat2 dummy` <- NULL
  }
  if(strat1_var == "strat1_dummy"){
    res$strat1_dummy <- NULL
    out$`strat1 dummy` <- NULL
  }
  
  ## Update header 
  hdr <- format_header(all_colnames = colnames(out), header_colnames = display_statistics, header_name = variable_names[num_var])
  
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
  
}






#' Kruskal-Wallis H test or Wilcoxon Rank-Sum test
#' 
#' Run Kruskal-Wallis H test or Wilcoxon Rank-Sum test for multiple covariates
#' 
#' @inheritParams wrapper_kruskal_test_core_col_cat_strat
#' @param num_vars Vector with names of numerical variables. If it has length >= 1, then 'cat_var' must be of length 1, and stratification subgroups are displayed in columns and statistics in rows.
#' @param cat_vars Vector with names of categorical variables. If it has length >= 1, then 'num_var' must be of length 1, and stratification subgroups are displayed in rows and statistics in columns.
#' @param display_in_column Possible values: "cat", "num".
#' @export
wrapper_kruskal_test <- function(data, num_vars, cat_vars, strat1_var = NULL, strat2_var = NULL, method = "kruskal", alternative = "two.sided", paired = FALSE, pairwise = FALSE, variable_names = NULL, caption = NULL, display_statistics = NULL, display_in_column = "cat", force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE, drop = FALSE){
  
  
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
      
      wrapper_res <- wrapper_kruskal_test_core_col_num_strat(data = data, num_var = num_var, cat_var = cat_var, strat1_var = strat1_var, strat2_var = strat2_var, method = method, alternative = alternative, paired = paired, pairwise = pairwise, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues, drop = drop)
      
      return(wrapper_res)
      
    })
    
    
  }else{
    
    
    wrapper_res <- lapply(1:length(num_vars), function(i){
      # i = 1
      
      num_var <- num_vars[i]
      cat_var <- cat_vars
      
      wrapper_res <- wrapper_kruskal_test_core_col_cat_strat(data = data, num_var = num_var, cat_var = cat_var, strat1_var = strat1_var, strat2_var = strat2_var, method = method, alternative = alternative, paired = paired, variable_names = variable_names, caption = caption, display_statistics = display_statistics, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues, drop = drop)
      
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






































