





#' Merge topTable results for multiple contrasts
#' 
#' @param fit Fit should be an object of class `MArrayLM` as produced by `lmFit` and `eBayes` where we can apply `topTable`.
#' @export
wrapper_merge_topTables <- function(fit, contrasts, gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), res_vars = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), sep = "_", pval = 0.05, lfc = c(0, 1)){
  
  
  topTable_out <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    contrast <- contrasts[i]
    
    
    topTable_out <- limma::topTable(fit, coef = contrast, sort.by = "none", number = Inf)
    
    topTable_out <- topTable_out[, c(gene_vars, res_vars), drop = FALSE]
    
    topTable_out$statistic <- -log10(topTable_out$P.Value) * sign(topTable_out$logFC)
    
    # plot(topTable_out$t, topTable_out$statistic)
    
      
    for(p in 1:length(pval)){
      
      for(l in 1:length(lfc)){
        
        topTable_out[, paste0("summary_pval", round(pval[p]*100), "_lfc", lfc[l])] <- as.numeric(limma::decideTests(fit, p.value = pval[p], lfc = lfc[l])[, contrast])
        
      }
      
    }
    
    colnames2change <- ! colnames(topTable_out) %in% gene_vars
    
    colnames(topTable_out)[colnames2change] <- paste(colnames(topTable_out)[colnames2change], contrast, sep = sep)
    
    
    return(topTable_out)
    
  })
  
  
  
  topTable_out <- Reduce(function(...) merge(..., by = gene_vars, all = TRUE, sort = FALSE), topTable_out)
  
  
  return(topTable_out)
  
  
}





#' Create a summary variable indicating significant genes
#' 
#' @param x Data frame of merged topTables
#' @export
wrapper_deside_tests <- function(x, topn = Inf, pval = 0.05, lfc = 0, 
  gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), 
  lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", summary_prefix = NULL, 
  sep = "_",
  order_both_direction = FALSE){
  
  
  stopifnot(length(gene_vars) >= 1)
  stopifnot(all(gene_vars %in% colnames(x)))
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", pval_prefix, sep), "", grep(paste0("^", pval_prefix, sep), colnames(x), value = TRUE))
  
  
  summary_out <- x[, gene_vars, drop = FALSE]
  
  
  for(i in 1:length(contrasts)){
    # i = 1
    
    contrast <- contrasts[i]
    
    pval_var <- paste0(pval_prefix, "_", contrast)
    adjp_var <- paste0(adjp_prefix, "_", contrast)
    lfc_var <- paste0(lfc_prefix, "_", contrast)
    
    
    summary_order <- rep(NA, nrow(x))
    
    
    if(order_both_direction){
      summary_order[!is.na(x[, lfc_var])] <- order(x[!is.na(x[, lfc_var]), pval_var], decreasing = FALSE)
    }else{
      summary_order[x[, lfc_var] >= lfc & !is.na(x[, lfc_var])] <- order(x[x[, lfc_var] >= lfc & !is.na(x[, lfc_var]), pval_var], decreasing = FALSE)
      
      summary_order[x[, lfc_var] <= lfc & !is.na(x[, lfc_var])] <- order(x[x[, lfc_var] <= lfc & !is.na(x[, lfc_var]), pval_var], decreasing = FALSE)
    }
    
    
    
    for(p in 1:length(pval)){
      
      for(l in 1:length(lfc)){
        # p = 1; l = 1
        
        
        if(is.null(summary_prefix)){
          name_out <- paste0("summary_pval", round(pval[p]*100), "_lfc", lfc[l], "_topn", topn, "_", contrast)
        }else{
          name_out <- paste0(summary_prefix, "_", contrast)
        }
        
        
        out <- as.numeric(x[, adjp_var] <= pval[p] & abs(x[, lfc_var]) >= lfc[l] & summary_order <= topn) * sign(x[, lfc_var])
        
        # table(out, useNA = "always")
        
        summary_out[, name_out] <- out
        
        
      }
      
    }
    
  }
  
  
  
  return(summary_out)
  
  
}







#' Extract given statistics for all contrasts available in merged topTable
#' 
#' @param x Data frame of merged topTables
#' @export
wrapper_extract_from_topTable <- function(x, contrasts = NULL, extract_prefix = "logFC", sep = "_"){
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts_identified <- gsub(paste0("^", extract_prefix, sep), "", grep(paste0("^", extract_prefix, sep), colnames(x), value = TRUE))
  
  if(!is.null(contrasts)){
    stopifnot(all(contrasts %in% contrasts_identified))
  }else{
    contrasts <- contrasts_identified
  }
  
  cols <- paste0(extract_prefix, sep, contrasts)
  out <- x[, cols, drop = FALSE]
  colnames(out) <- contrasts
  
  out
  
}




#' Convert results from Cox regression into a topTable
#' 
#' @param x "BclassTesting" object, for example, output of the wrapper_cox_regression_biomarker function. It does not work with log-rank test results.
#' @export
wrapper_bresults_to_topTable <- function(x, contrast_vars, id_cols = "biomarker", statistic_change = "HR", readjust_pvalues = TRUE){
  
  res <- bresults(x)
  
  statistic_change_CI95_lower <- paste0(statistic_change, "_CI95_lower")
  statistic_change_CI95_upper <- paste0(statistic_change, "_CI95_upper")
  
  log_statistic_change <- paste0("log", statistic_change)
  log_statistic_change_CI95_lower <- paste0("log", statistic_change, "_CI95_lower")
  log_statistic_change_CI95_upper <- paste0("log", statistic_change, "_CI95_upper")
  
  
  stopifnot(all(c(statistic_change, statistic_change_CI95_lower, statistic_change_CI95_upper, log_statistic_change, log_statistic_change_CI95_lower, log_statistic_change_CI95_upper, "pvalue", "sign_pvalue") %in% colnames(res)))
  
  if(!readjust_pvalues){
    stopifnot(all(c("adj_pvalue") %in% colnames(res)))
  }
  
  
  if(paste0(statistic_change, "_non_empty") %in% colnames(res)){
    res <- res[res[, paste0(statistic_change, "_non_empty")], ]
  }
  
  
  if(is.null(contrast_vars)){
    res$pooled <- "pooled"
    contrast_vars <- "pooled"
  }
  
  
  ### Format contrasts
  
  for(i in seq_along(contrast_vars)){
    # i = 1
    
    res[, contrast_vars[i]] <- factor(res[, contrast_vars[i]], levels = unique(res[, contrast_vars[i]]))
    
  }
  
  res$contrast <- interaction(res[, contrast_vars, drop = FALSE], sep = "_", lex.order = TRUE)
  
  res <- res[order(res$contrast), ]
  
  table(res$contrast)
  
  
  ### Extract statistics 
  
  
  if(readjust_pvalues){
    
    ### Re-adjust p-values per contrast 
    
    res$adj_pvalue <- res$pvalue
    
    topTable <- pivot_wider(res, id_cols = all_of(id_cols), names_from = all_of("contrast"), values_from = all_of(c(statistic_change, statistic_change_CI95_lower, statistic_change_CI95_upper, log_statistic_change, log_statistic_change_CI95_lower, log_statistic_change_CI95_upper, "pvalue", "sign_pvalue", "adj_pvalue")))
    
    topTable <- mutate_at(topTable, grep("^adj_pvalue", colnames(topTable)), stats::p.adjust, method = "BH")

    
  }else{
    
    
    topTable <- pivot_wider(res, id_cols = all_of(id_cols), names_from = all_of("contrast"), values_from = all_of(c(statistic_change, statistic_change_CI95_lower, statistic_change_CI95_upper, log_statistic_change, log_statistic_change_CI95_lower, log_statistic_change_CI95_upper, "pvalue", "sign_pvalue", "adj_pvalue")))
    
    
  }
  

  ### Update column names to be aligned with what is produced in DGE 
  
  topTable <- data.frame(topTable, check.names = FALSE)
  
  colnames(topTable) <- gsub("^pvalue_", "P.Value_", colnames(topTable))
  colnames(topTable) <- gsub("^adj_pvalue_", "adj.P.Val_", colnames(topTable))
  colnames(topTable) <- gsub("^sign_pvalue_", "sign.P.Val_", colnames(topTable))
  colnames(topTable) <- gsub(paste0("^", statistic_change, "_CI95_"), paste0(statistic_change, ".CI95."), colnames(topTable))
  colnames(topTable) <- gsub(paste0("^", log_statistic_change, "_CI95_"), paste0(log_statistic_change, ".CI95."), colnames(topTable))
  
  
  return(topTable)
  
}






# x <- topTable_DE_Ipat_Dx_Group_endpoints
# contrast <- contrasts[1]
# 
# direction = "up"
# topn = 10
# pval = 0.5
# lfc = 0
# gene_vars = c("Endpoint")
# lfc_prefix = "logFC"
# pval_prefix = "P.Value"
# adjp_prefix = "adj.P.Val"
# stats_prefixes = NULL
# caption = NULL
# sep = "_"




#' Display significantly DE genes
#' 
#' @param x Data frame of merged topTables
#' @param sort_by Possible values: "none", "pval", "lfc".
#' @param direction Possible values: "up", "down", "both".
#' @export
wrapper_dispaly_significant_genes <- function(x, contrast, direction = "up", 
  sort_by = "pval", topn = 20, pval = 0.05, lfc = 0, 
  gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), 
  lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = NULL, sep = "_", 
  caption = NULL){
  
  # -------------------------------------------------------------------------
  # Checks
  # -------------------------------------------------------------------------
  
  stopifnot(length(gene_vars) >= 1)
  stopifnot(all(gene_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  stopifnot(sort_by %in% c("none", "pval", "lfc"))
  
  stopifnot(direction %in% c("up", "down", "both"))
  
  if(direction == "both"){
    direction_print <- "up-, down-"
  }else{
    direction_print <- paste0(direction, "-")  
  }
  
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", pval_prefix, sep), "", grep(paste0("^", pval_prefix, sep), colnames(x), value = TRUE))
  
  stopifnot(contrast %in% contrasts)
  
  
  # -------------------------------------------------------------------------
  # Processing
  # -------------------------------------------------------------------------
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(lfc_prefix, stats_prefixes, unique(c(pval_prefix, adjp_prefix))), sep, contrast)
  
  x <- x[ , c(gene_vars, contrast_vars_display), drop = FALSE]
  
  colnames(x) <- gsub(paste0(sep, contrast, "$"), "", colnames(x))
  
  
  ## Sort
  if(sort_by == "pval"){
    x_sort <- x[order(x[, pval_prefix], decreasing = FALSE), , drop = FALSE]
  }else if(sort_by == "lfc"){
    x_sort <- x[order(x[, lfc_prefix], decreasing = TRUE), , drop = FALSE]
  }else{
    x_sort <- x
  }
  
  
  ## Subset by adj. p-value
  x_sort <- x_sort[x_sort[, adjp_prefix] <= pval & !is.na(x_sort[, adjp_prefix]), , drop = FALSE]
  
  
  ## Subset by direction
  if(direction == "up"){
    x_sort <- x_sort[x_sort[, lfc_prefix] >= lfc, , drop = FALSE]
  }else if(direction == "down"){
    x_sort <- x_sort[x_sort[, lfc_prefix] <= -lfc, , drop = FALSE]
  }else{
    x_sort <- x_sort[abs(x_sort[, lfc_prefix]) >= lfc, , drop = FALSE]
  }
  
  
  if(nrow(x_sort) == 0){
    
    caption <- paste0("There are no ", direction_print, "regulated genes (", adjp_prefix, " <= ", pval, ", |", lfc_prefix, "| >= ", lfc, ") when testing for ", contrast, ".")
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
    return(BclassDE(caption = caption))
    
  }
  
  
  
  res <- x_sort[1:(min(nrow(x_sort), topn)), , drop = FALSE]
  rownames(res) <- NULL
  
  out <- res %>% 
    mutate_at(lfc_prefix, format_difference) %>% 
    mutate_at(unique(c(pval_prefix, adjp_prefix)), format_pvalues2) 
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("List of ", direction_print, "regulated genes (", adjp_prefix, " <= ", pval, ", |", lfc_prefix, "| >= ", lfc, ") when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " genes.")
      
    }
    
  }
  
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}




##############################################################################
# lm
##############################################################################


# require(graphics)
# 
# ## Annette Dobson (1990) "An Introduction to Generalized Linear Models".
# ## Page 9: Plant Weight Data.
# ctl <- c(4.17,5.58,5.18,6.11,4.50,4.61,5.17,4.53,5.33,5.14)
# trt <- c(4.81,4.17,4.41,3.59,5.87,3.83,6.03,4.89,4.32,4.69)
# group <- gl(2, 10, 20, labels = c("Ctl","Trt"))
# weight <- c(ctl, trt)
# lm.D9 <- lm(weight ~ group)
# summary(lm.D9)
# 
# glht_tmp <- multcomp::glht(lm.D9, linfct = t(c(0, 1)))
# summary(glht_tmp)



# data <- pdata
# biomarker_vars <- names(gene_signature_collection[[1]])
# formula <- stats::as.formula("~ 0 + DE_var")


#' Fit a linear model 
#' 
#' @param data Data frame.
#' @param biomarker_vars Variables defining dependent variables.
#' @export
wrapper_lm <- function(data, biomarker_vars, formula, contrast_matrix){
  
  
  stopifnot(length(biomarker_vars) >= 1)
  stopifnot(all(sapply(data[, biomarker_vars], class) == "numeric"))
  
  contrasts <- colnames(contrast_matrix)
  
  
  ## Fit lm for each biomarker separately
  
  out <- lapply(seq_along(biomarker_vars), function(i){
    # i = 1
    
    data$dummy_biomarker_var <- data[, biomarker_vars[i]]
    
    
    formula_tmp <- stats::as.formula(paste0("dummy_biomarker_var", paste0(as.character(formula), collapse = "")))
    
    fit <- stats::lm(formula_tmp, data)
    
    
    ## Fit contrasts one by one
    out_glht <- lapply(seq_along(contrasts), function(k){
      # k = 1
      
      glht_tmp <- multcomp::glht(fit, linfct = t(contrast_matrix[, k, drop = FALSE]))
      
      summary_tmp <- summary(glht_tmp)
      
      
      pval <- summary_tmp$test$pvalues
      
      difference <- as.numeric(summary_tmp$test$coefficients)
      
      out <- data.frame(biomarker = biomarker_vars[i], difference = difference, pval = pval, stringsAsFactors = FALSE)
      colnames(out) <- c("biomarker", paste0("Difference_", contrasts[k]), paste0("P.Value_", contrasts[k]))
      
      out
      
    })
    
    out <- do.call(cbind, out_glht)
    out
    
    
  })
  
  
  out <- plyr::rbind.fill(out)
  
  
  ## Adjust the p-values
  adjp <- data.frame(lapply(out[, paste0("P.Value_", contrasts), drop = FALSE], stats::p.adjust, method = "BH"))
  
  colnames(adjp) <- paste0("adj.P.Val_", contrasts)
  
  out <- cbind(out, adjp)
  
  column_order <- c("biomarker", apply(expand.grid(c("Difference_", "P.Value_", "adj.P.Val_"), contrasts, stringsAsFactors = FALSE), 1, paste0, collapse = ""))
  
  out <- out[, column_order, drop = FALSE]
  
  out
  
  
}





##############################################################################
# lmer
##############################################################################


# data <- pdata
# biomarker_vars <- names(gene_signature_collection[[1]])
# formula <- stats::as.formula("~ 0 + DE_var + (1|PatientID)")


#' Fit a linear mixed model 
#' 
#' @param data Data frame.
#' @param biomarker_vars Variables defining dependent variables.
#' @export
wrapper_lmer <- function(data, biomarker_vars, formula, contrast_matrix){
  
  
  stopifnot(length(biomarker_vars) >= 1)
  stopifnot(all(sapply(data[, biomarker_vars], class) == "numeric"))
  
  contrasts <- colnames(contrast_matrix)
  
  
  ## Fit lm for each biomarker separately
  
  out <- lapply(seq_along(biomarker_vars), function(i){
    # i = 1
    
    data$dummy_biomarker_var <- data[, biomarker_vars[i]]
    
    formula_tmp <- stats::as.formula(paste0("dummy_biomarker_var", paste0(as.character(formula), collapse = "")))
    
    fit <- lme4::lmer(formula_tmp, data)
    
    
    ## Fit contrasts one by one
    out_glht <- lapply(seq_along(contrasts), function(k){
      # k = 1
      
      glht_tmp <- multcomp::glht(fit, linfct = t(contrast_matrix[, k, drop = FALSE]))
      
      summary_tmp <- summary(glht_tmp)
      
      
      pval <- summary_tmp$test$pvalues
      
      difference <- as.numeric(summary_tmp$test$coefficients)
      
      out <- data.frame(biomarker = biomarker_vars[i], difference = difference, pval = pval, stringsAsFactors = FALSE)
      colnames(out) <- c("biomarker", paste0("Difference_", contrasts[k]), paste0("P.Value_", contrasts[k]))
      
      out
      
    })
    
    out <- do.call(cbind, out_glht)
    out
    
    
  })
  
  
  out <- plyr::rbind.fill(out)
  
  
  ## Adjust the p-values
  adjp <- data.frame(lapply(out[, paste0("P.Value_", contrasts), drop = FALSE], stats::p.adjust, method = "BH"))
  
  colnames(adjp) <- paste0("adj.P.Val_", contrasts)
  
  out <- cbind(out, adjp)
  
  column_order <- c("biomarker", apply(expand.grid(c("Difference_", "P.Value_", "adj.P.Val_"), contrasts, stringsAsFactors = FALSE), 1, paste0, collapse = ""))
  
  out <- out[, column_order, drop = FALSE]
  
  out
  
  
}
























































