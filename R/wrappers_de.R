





#' Merge topTable results for multiple contrasts
#' 
#' @param fit Fit should be an object of class MArrayLM as produced by lmFit and eBayes where we can apply topTable.
wrapper_merge_topTables <- function(fit, contrasts, gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), res_vars = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), sep = "_", pval = 0.05, lfc = c(0, 1)){
  
  
  topTable_out <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    contrast <- contrasts[i]
    
    
    topTable_out <- topTable(fit, coef = contrast, sort.by = "none", number = Inf)
    
    topTable_out <- topTable_out[, c(gene_vars, res_vars), drop = FALSE]
    
    
    for(p in 1:length(pval)){
      
      for(l in 1:length(lfc)){
        
        topTable_out[, paste0("summary_pval", round(pval[p]*100), "_lfc", lfc[l])] <- as.numeric(decideTests(fit, p.value = pval[p], lfc = lfc[l])[, contrast])
        
      }
      
    }
    
    summary_vars <- grep("^summary_", colnames(topTable_out), value = TRUE)
    
    colnames2change <- colnames(topTable_out) %in% c(res_vars, summary_vars)
    
    colnames(topTable_out)[colnames2change] <- paste(colnames(topTable_out)[colnames2change], contrast, sep = sep)
    
    
    return(topTable_out)
    
  })
  
  
  
  topTable_out <- Reduce(function(...) merge(..., by = gene_vars, all = TRUE, sort = FALSE), topTable_out)
  
  
  return(topTable_out)
  
  
}





wrapper_extract_from_topTable <- function(x, extract_prefix = "logFC", sep = "_"){
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", extract_prefix, sep), "", grep(paste0("^", extract_prefix, sep), colnames(x), value = TRUE))
  
  cols <- paste0(extract_prefix, sep, contrasts)
  out <- x[, cols, drop = FALSE]
  colnames(out) <- contrasts
  
  out
  
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




wrapper_dispaly_significant_genes <- function(x, contrast, direction = "up", 
  sort_by = "pval", topn = 20, pval = 0.05, lfc = 0, 
  gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
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
  contrasts <- gsub(paste0("^", lfc_prefix, sep), "", grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE))
  
  stopifnot(contrast %in% contrasts)
  
  
  # -------------------------------------------------------------------------
  # Processing
  # -------------------------------------------------------------------------
  
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(lfc_prefix, stats_prefixes, pval_prefix, adjp_prefix), sep, contrast)
  
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
  x_sort <- x_sort[x_sort[, adjp_prefix] <= pval, , drop = FALSE]
  
  
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
    mutate_at(pval_prefix, format_pvalues2) %>% 
    mutate_at(adjp_prefix, format_pvalues2)
  
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("List of ", direction_print, "regulated genes (", adjp_prefix, " <= ", pval, ", |", lfc_prefix, "| >= ", lfc, ") when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " genes.")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}




##############################################################################
# lm
##############################################################################


# data <- pdata
# biomarker_vars <- names(gene_signature_collection[[1]])
# formula <- as.formula("~ 0 + DE_var")



wrapper_lm <- function(data, biomarker_vars, formula, contrast_matrix){
  
  
  stopifnot(length(biomarker_vars) >= 1)
  stopifnot(all(sapply(data[, biomarker_vars], class) == "numeric"))
  
  contrasts <- colnames(contrast_matrix)
  
  
  ## Fit lm for each biomarker separately
  
  out <- lapply(seq_along(biomarker_vars), function(i){
    # i = 1
    
    formula_tmp <- as.formula(paste0(biomarker_vars[i], paste0(as.character(formula), collapse = "")))
    
    fit <- lm(formula_tmp, data)
    
    
    ## Fit contrasts one by one
    out_glht <- lapply(seq_along(contrasts), function(k){
      # k = 1
      
      glht_tmp <- multcomp::glht(fit, linfct = t(contrast_matrix[, k, drop = FALSE]))
      
      summary_tmp <- summary(glht_tmp)
      
      
      pval <- summary_tmp$test$pvalues
      
      lfc <- as.numeric(summary_tmp$test$coefficients)
      
      out <- data.frame(biomarker = biomarker_vars[i], lfc = lfc, pval = pval, stringsAsFactors = FALSE)
      colnames(out) <- c("biomarker", paste0("logFC_", contrasts[k]), paste0("P.Value_", contrasts[k]))
      
      out
      
    })
    
    out <- do.call(cbind, out_glht)
    out
    
    
  })
  
  
  out <- plyr::rbind.fill(out)
  
  
  ## Adjust the p-values
  adjp <- data.frame(lapply(out[, paste0("P.Value_", contrasts), drop = FALSE], p.adjust, method = "BH"))
  
  colnames(adjp) <- paste0("adj.P.Val_", contrasts)
  
  out <- cbind(out, adjp)
  
  column_order <- c("biomarker", apply(expand.grid(c("logFC_", "P.Value_", "adj.P.Val_"), contrasts, stringsAsFactors = FALSE), 1, paste0, collapse = ""))
  
  out <- out[, column_order, drop = FALSE]
  
  out
  
  
}





##############################################################################
# lmer
##############################################################################


# data <- pdata
# biomarker_vars <- names(gene_signature_collection[[1]])
# formula <- as.formula("~ 0 + DE_var + (1|PatientID)")



wrapper_lmer <- function(data, biomarker_vars, formula, contrast_matrix){
  
  
  stopifnot(length(biomarker_vars) >= 1)
  stopifnot(all(sapply(data[, biomarker_vars], class) == "numeric"))
  
  contrasts <- colnames(contrast_matrix)
  
  
  ## Fit lm for each biomarker separately
  
  out <- lapply(seq_along(biomarker_vars), function(i){
    # i = 1
    
    formula_tmp <- as.formula(paste0(biomarker_vars[i], paste0(as.character(formula), collapse = "")))
    
    fit <- lme4::lmer(formula_tmp, data)
    
    
    ## Fit contrasts one by one
    out_glht <- lapply(seq_along(contrasts), function(k){
      # k = 1
      
      glht_tmp <- multcomp::glht(fit, linfct = t(contrast_matrix[, k, drop = FALSE]))
      
      summary_tmp <- summary(glht_tmp)
      
      
      pval <- summary_tmp$test$pvalues
      
      lfc <- as.numeric(summary_tmp$test$coefficients)
      
      out <- data.frame(biomarker = biomarker_vars[i], lfc = lfc, pval = pval, stringsAsFactors = FALSE)
      colnames(out) <- c("biomarker", paste0("logFC_", contrasts[k]), paste0("P.Value_", contrasts[k]))
      
      out
      
    })
    
    out <- do.call(cbind, out_glht)
    out
    
    
  })
  
  
  out <- plyr::rbind.fill(out)
  
  
  ## Adjust the p-values
  adjp <- data.frame(lapply(out[, paste0("P.Value_", contrasts), drop = FALSE], p.adjust, method = "BH"))
  
  colnames(adjp) <- paste0("adj.P.Val_", contrasts)
  
  out <- cbind(out, adjp)
  
  column_order <- c("biomarker", apply(expand.grid(c("logFC_", "P.Value_", "adj.P.Val_"), contrasts, stringsAsFactors = FALSE), 1, paste0, collapse = ""))
  
  out <- out[, column_order, drop = FALSE]
  
  out
  
  
}










##############################################################################
# Beta regression
##############################################################################

# library(betareg)
# library(glmmTMB)



# library(glmmADMB)
# try(fit_tmp <- glmmADMB::glmmadmb(y/total ~ response + day + response:day, family = "beta", data = data_tmp), silent = TRUE)
# try(fit_tmp <- glmmADMB::glmmadmb(prop ~ response + day + response:day + (1|patient_id), family = "beta", data = data_tmp), silent = TRUE)



# biomarker_vars <- names(gene_signature_collection[[1]])
# formula <- as.formula("~ 0 + DE_var + (1|PatientID)")



wrapper_glmm_beta <- function(data, biomarker_vars, formula, contrast_matrix){
  
  
  stopifnot(length(biomarker_vars) >= 1)
  stopifnot(all(sapply(data[, biomarker_vars], class) == "numeric"))
  stopifnot(all(sapply(data[, biomarker_vars], function(x){all(x[!is.na(x)] > 0 & x[!is.na(x)] < 1)})))
  
  
  contrasts <- colnames(contrast_matrix)
  
  
  ## Fit beta mixed model for each biomarker separately
  
  out <- lapply(seq_along(biomarker_vars), function(i){
    # i = 1
    
    data_tmp <- data[complete.cases(data[, biomarker_vars[i]]), , drop = FALSE]
    
    formula_tmp <- as.formula(paste0(biomarker_vars[i], paste0(as.character(formula), collapse = "")))
    
    
    fit <- NULL
    
    ### multcomp::glht does not work for glmmTMB
    # fit <- glmmTMB::glmmTMB(formula = formula_tmp, family = beta_family(link = "logit"), data = data_tmp)
    
    try(fit <- glmmADMB::glmmadmb(formula = formula_tmp, family = "beta", link = "logit", data = data_tmp), silent = TRUE)
    
    
    if(!is.null(fit)){
      
      ## Fit contrasts one by one
      out_glht <- lapply(seq_along(contrasts), function(k){
        # k = 1
        
        glht_tmp <- multcomp::glht(fit, linfct = t(contrast_matrix[, k, drop = FALSE]))
        
        summary_tmp <- summary(glht_tmp)
        
        
        pval <- summary_tmp$test$pvalues
        
        lfc <- as.numeric(summary_tmp$test$coefficients)
        
        out <- data.frame(biomarker = biomarker_vars[i], lfc = lfc, pval = pval, stringsAsFactors = FALSE)
        colnames(out) <- c("biomarker", paste0("logOR_", contrasts[k]), paste0("P.Value_", contrasts[k]))
        
        out
        
      })
      
      out <- do.call(cbind, out_glht)
      
    }else{
      
      out_glht <- lapply(seq_along(contrasts), function(k){
        # k = 1
        
        out <- data.frame(biomarker = biomarker_vars[i], lfc = NA, pval = NA, stringsAsFactors = FALSE)
        colnames(out) <- c("biomarker", paste0("logOR_", contrasts[k]), paste0("P.Value_", contrasts[k]))
        
        out
        
      })
      
      out <- do.call(cbind, out_glht)
      
    }
    
    
    out
    
    
  })
  
  
  out <- plyr::rbind.fill(out)
  
  
  ## Adjust the p-values
  adjp <- data.frame(lapply(out[, paste0("P.Value_", contrasts), drop = FALSE], p.adjust, method = "BH"))
  
  colnames(adjp) <- paste0("adj.P.Val_", contrasts)
  
  out <- cbind(out, adjp)
  
  column_order <- c("biomarker", apply(expand.grid(c("logOR_", "P.Value_", "adj.P.Val_"), contrasts, stringsAsFactors = FALSE), 1, paste0, collapse = ""))
  
  out <- out[, column_order, drop = FALSE]
  
  out
  
  
}





wrapper_betareg <- function(data, biomarker_vars, formula, contrast_matrix){
  
  
  stopifnot(length(biomarker_vars) >= 1)
  stopifnot(all(sapply(data[, biomarker_vars], class) == "numeric"))
  stopifnot(all(sapply(data[, biomarker_vars], function(x){all(x[!is.na(x)] > 0 & x[!is.na(x)] < 1)})))
  
  
  contrasts <- colnames(contrast_matrix)
  
  
  ## Fit beta mixed model for each biomarker separately
  
  out <- lapply(seq_along(biomarker_vars), function(i){
    # i = 1
    
    data_tmp <- data[complete.cases(data[, biomarker_vars[i]]), , drop = FALSE]
    
    formula_tmp <- as.formula(paste0(biomarker_vars[i], paste0(as.character(formula), collapse = "")))
    
    
    fit <- betareg::betareg(formula = formula_tmp, link = "logit", data = data_tmp)
    
    
    ## Fit contrasts one by one
    out_glht <- lapply(seq_along(contrasts), function(k){
      # k = 1
      
      glht_tmp <- multcomp::glht(fit, linfct = t(contrast_matrix[, k, drop = FALSE]))
      
      summary_tmp <- summary(glht_tmp)
      
      
      pval <- summary_tmp$test$pvalues
      
      lfc <- as.numeric(summary_tmp$test$coefficients)
      
      out <- data.frame(biomarker = biomarker_vars[i], lfc = lfc, pval = pval, stringsAsFactors = FALSE)
      colnames(out) <- c("biomarker", paste0("logOR_", contrasts[k]), paste0("P.Value_", contrasts[k]))
      
      out
      
    })
    
    out <- do.call(cbind, out_glht)
    out
    
    
  })
  
  
  out <- plyr::rbind.fill(out)
  
  
  ## Adjust the p-values
  adjp <- data.frame(lapply(out[, paste0("P.Value_", contrasts), drop = FALSE], p.adjust, method = "BH"))
  
  colnames(adjp) <- paste0("adj.P.Val_", contrasts)
  
  out <- cbind(out, adjp)
  
  column_order <- c("biomarker", apply(expand.grid(c("logOR_", "P.Value_", "adj.P.Val_"), contrasts, stringsAsFactors = FALSE), 1, paste0, collapse = ""))
  
  out <- out[, column_order, drop = FALSE]
  
  out
  
  
}














































