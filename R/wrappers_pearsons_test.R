




# variable_names = NULL; caption = NULL; margin = 1; force_empty_cols = FALSE; print_pvalues = TRUE




#' Pearson’s Chi-squared test
#' 
#' @param data Data frame.
#' @export
wrapper_core_pearsons_test <- function(data, col_var, row_var, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_pvalues = TRUE){
  
  # --------------------------------------------------------------------------
  # Check about input data and some preprocessing
  # --------------------------------------------------------------------------
  
  stopifnot(is.data.frame(data))
  
  stopifnot(length(col_var) == 1)
  stopifnot(is.factor(data[, col_var]))
  stopifnot(nlevels(data[, col_var]) == 2)
  
  stopifnot(length(row_var) == 1)
  stopifnot(is.factor(data[, row_var]))
  stopifnot(nlevels(data[, row_var]) >= 2)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  
  # --------------------------------------------------------------------------
  # Calculate counts and proportions and do testing
  # --------------------------------------------------------------------------
  
  tbl <- table(data[, row_var], data[, col_var])
  
  margin <- 1
  
  prop <- prop.table(tbl, margin = margin) * 100
  
  tbl_test <- tbl[, rev(seq_len(ncol(tbl)))]
  
  success_level <- colnames(tbl_test)[1]
  
  
  if(sum(margin.table(tbl, margin = 1) >= 1) >= 2 && sum(margin.table(tbl, margin = 2) >= 1) >= 2){
    ## Pearson’s Chi-squared test: get difference of proportions of succeses and CI for difference
    ## Success category is defined by the first column
    
    test_res <- prop.test(tbl_test)
    
    pvalue <- test_res$p.value
    
    if(!is.null(test_res$conf.int)){
      difference <- (test_res$estimate[2] - test_res$estimate[1]) * 100
      CI95_upper <- -test_res$conf.int[1] * 100
      CI95_lower <- -test_res$conf.int[2] * 100
    }else{
      difference <- NA
      CI95_upper <- NA
      CI95_lower <- NA
    }
    

    
  }else{
    
    pvalue <- NA
    difference <- NA
    CI95_lower <- NA
    CI95_upper <- NA
    
  }
  
  
  ### Calculate CIs for proportions of success 
  
  cidf <- lapply(1:nrow(tbl_test), function(k){
    # k = 1
    
    res_binom <- binom.test(as.numeric(tbl_test[k, ]), alternative = "two.sided", conf.level = 0.95)
    
    cidf <- data.frame(CI95_lower = res_binom$conf.int[1] * 100, CI95_upper = res_binom$conf.int[2] * 100)
    
  })
  
  cidf <- rbind.fill(cidf)
  colnames(cidf) <- paste0(c("CI95_lower", "CI95_upper"), "_", success_level)
  
  
  # --------------------------------------------------------------------------
  # Prepare 'res' data frame
  # --------------------------------------------------------------------------
  
  countdf <- as.data.frame.matrix(tbl)
  colnames(countdf) <- paste0("counts_", colnames(countdf))
  
  propdf <- as.data.frame.matrix(prop)
  colnames(propdf) <- paste0("proportions_", colnames(propdf))
  
  
  res <- data.frame(covariate = row_var,
    subgroup = rownames(tbl),
    countdf, propdf, cidf,
    difference = c(difference, rep(NA, nrow(tbl) - 1)),
    CI95_lower = c(CI95_lower, rep(NA, nrow(tbl) - 1)),
    CI95_upper = c(CI95_upper, rep(NA, nrow(tbl) - 1)),
    pvalue = c(pvalue, rep(NA, nrow(tbl) - 1)), 
    stringsAsFactors = FALSE, row.names = NULL, check.names = FALSE)
  
  
  
  
  # --------------------------------------------------------------------------
  # Prepare 'out' data frame
  # --------------------------------------------------------------------------
  
  out <- data.frame(Covariate = variable_names[res$covariate], 
    Subgroup = res$subgroup, 
    
    format_counts_and_props(counts = countdf, props = propdf, digits = 1),
    
    format_CIs.data.frame(cidf, digits = 1, colname = paste0("95% CI ", success_level)),
    
    Difference = format_difference(res$difference, digits = 1),
    
    `95% CI` = format_CIs(res$CI95_lower, res$CI95_upper, digits = 1),
    
    `P-value` = format_pvalues(res$pvalue), 
    
    check.names = FALSE, stringsAsFactors = FALSE)
  
  
  stopifnot(all(sapply(out, class) == "character"))
  
  
  if(!print_pvalues){
    out$`P-value` <- NULL
  }
  
  
  ### If all Difference are empty, do not display that column.
  if(all(out$Difference == "") && !force_empty_cols){
    out$Difference <- NULL
    out$`95% CI` <- NULL
  }
  
  
  ### Set repeating Covariate names to empty
  out$Covariate[indicate_blocks(out, block_vars = "Covariate", return = "empty")] <- ""
  
  
  # --------------------------------------------------------------------------
  # Prepare 'header' data frame
  # --------------------------------------------------------------------------
  
  num_start_cols <- 2
  num_end_cols <- sum(c(paste0("95% CI ", success_level), "Difference", "95% CI", "P-value") %in% colnames(out))
  
  
  header <- c(num_start_cols, nlevels(data[, col_var]), num_end_cols)
  header <- as.integer(header)
  names(header) <- c(" ", variable_names[col_var], " ")
  
  
  # --------------------------------------------------------------------------
  ### Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("Pearson’s Chi-squared test.")
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassTesting(results = res, output = out, caption = caption, header = header)
  
  
  return(bout)
  
  
}





#' @inheritParams wrapper_core_pearsons_test
#' 
#' @param strat1_var Name of the firts stratification variable.
#' @param strat1_var Name of the second stratification variable.
#' @export
wrapper_core_pearsons_test_strat <- function(data, col_var, row_var, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
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
  
  ### Strata cannot include col_var, row_var
  stopifnot(length(intersect(c(strat1_var, strat2_var), c(col_var, row_var))) == 0)
  
  
  ### Keep non-missing data
  all_vars <- c(col_var, row_var, strat1_var, strat2_var)
  data <- data[complete.cases(data[, all_vars]), all_vars, drop = FALSE]
  
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
      
      
      wrapper_res <- wrapper_core_pearsons_test(data = data_strata1, col_var = col_var, row_var = row_var, variable_names = variable_names, caption = caption, force_empty_cols = force_empty_cols, print_pvalues = print_pvalues)
      
      
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
      
      ## Update header by adding 2 corresponding to the two strat variables to the first position
      hdr <- bheader(wrapper_res)
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
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
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
  
  ## Update header by adding 2 corresponding to the two strat variables to the first position
  hdr <- bheader(wrapper_res[[1]])
  hdr[1] <- hdr[1] - hdr_shift
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
  
}







#' @inheritParams wrapper_core_pearsons_test_strat
#' 
#' Pearson’s Chi-squared test
#' 
#' @param row_vars Vector with names of categorical variables.
#' @export
wrapper_pearsons_test <- function(data, col_var, row_vars, strat1_var = NULL, strat2_var = NULL, variable_names = NULL, caption = NULL, force_empty_cols = FALSE, print_pvalues = TRUE, print_adjpvalues = TRUE){
  
  
  # --------------------------------------------------------------------------
  # Checks
  # --------------------------------------------------------------------------
  
  
  stopifnot(length(col_var) == 1)
  stopifnot(length(row_vars) >= 1)
  
  variable_names <- format_variable_names(data = data, variable_names = variable_names)
  
  
  # --------------------------------------------------------------------------
  # Generate the results
  # --------------------------------------------------------------------------
  
  
  wrapper_res <- lapply(1:length(row_vars), function(i){
    # i = 1
    
    row_var <- row_vars[i]
    
    wrapper_res <- wrapper_core_pearsons_test_strat(data, col_var = col_var, row_var = row_var, strat1_var = strat1_var, strat2_var = strat2_var, variable_names = variable_names, caption = caption, force_empty_cols = TRUE, print_pvalues = print_pvalues, print_adjpvalues = print_adjpvalues)
    
    return(wrapper_res)
    
  })
  
  
  ### Merge the results
  
  res <- plyr::rbind.fill(lapply(wrapper_res, bresults))
  out <- plyr::rbind.fill(lapply(wrapper_res, boutput))
  
  
  ## Re-calculate adjusted p-values using the Benjamini & Hochberg method
  res$adj_pvalue <- p.adjust(res$pvalue, method = "BH")
  
  if("Adj. P-value" %in% colnames(out)){
    out$`Adj. P-value` <- format_pvalues(p.adjust(res$pvalue, method = "BH"))
  }
  
  
  hdr <- bheader(wrapper_res[[1]])
  
  ### Replace NAs with "" for "Difference", "95% CI"
  missing_columns <- c("Difference", "95% CI")
  
  for(i in 1:length(missing_columns)){
    if(missing_columns[i] %in% colnames(out)){
      out[is.na(out[, missing_columns[i]]), missing_columns[i]] <- ""
      ### If all are empty, do not display that column.
      if(all(out[, missing_columns[i]] == "") && !force_empty_cols){
        out[, missing_columns[i]] <- NULL
        ### Update header
        hdr[length(hdr)] <- hdr[length(hdr)] - 1
      }
    }
  }
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  wrapper_res <- BclassTesting(results = res, output = out, caption = bcaption(wrapper_res[[1]]), header = hdr)
  
  
  return(wrapper_res)
  
}















