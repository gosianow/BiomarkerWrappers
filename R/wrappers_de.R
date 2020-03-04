





#' Merge topTable results for multiple contrasts
#' 
#' @param fit Fit should be an object of class MArrayLM as produced by lmFit and eBayes where we can apply topTable.
wrapper_merge_topTables <- function(fit, contrasts, gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), res_vars = c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val"), sep = "_"){
  
  
  topTable_out <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    contrast <- contrasts[i]
    
    
    topTable_out <- topTable(fit, coef = contrast, sort.by = "none", number = Inf)
    
    topTable_out <- topTable_out[, c(gene_vars, res_vars), drop = FALSE]
    
    topTable_out$summary_LFC0 <- as.numeric(decideTests(fit, p.value = 0.05, lfc = 0)[, contrast])
    topTable_out$summary_LFC1 <- as.numeric(decideTests(fit, p.value = 0.05, lfc = 1)[, contrast])
    
    colnames2change <- colnames(topTable_out) %in% c(res_vars, "summary_LFC0", "summary_LFC1")
    
    colnames(topTable_out)[colnames2change] <- paste(colnames(topTable_out)[colnames2change], contrast, sep = sep)
    
    
    return(topTable_out)
    
  })
  
  
  
  topTable_out <- Reduce(function(...) merge(..., by = gene_vars, all = TRUE, sort = FALSE), topTable_out)
  
  
  return(topTable_out)
  
  
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
# description_prefixes = NULL
# caption = NULL
# sep = "_"



wrapper_dispaly_significant_genes <- function(x, contrast, direction = "up", topn = 20, pval = 0.05, lfc = 0, gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", stats_prefixes = NULL, description_prefixes = NULL, caption = NULL,  sep = "_"){
  
  stopifnot(length(gene_vars) >= 1)
  stopifnot(all(gene_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", lfc_prefix, sep), "", grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE))
  
  stopifnot(contrast %in% contrasts)
  
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(lfc_prefix, stats_prefixes, description_prefixes, pval_prefix, adjp_prefix), sep, contrast)
  
  x <- x[ , c(gene_vars, contrast_vars_display), drop = FALSE]
  
  colnames(x) <- gsub(paste0(sep, contrast, "$"), "", colnames(x))
  
  ## Sort by p-value
  x_sort <- x[order(x[, pval_prefix], decreasing = FALSE), , drop = FALSE]
  ## Subset by adj. p-value
  x_sort <- x_sort[x_sort[, adjp_prefix] < pval, , drop = FALSE]
  ## Subset by LogFC
  if(direction == "up"){
    x_sort <- x_sort[x_sort[, lfc_prefix] > lfc, , drop = FALSE]
  }else{
    x_sort <- x_sort[x_sort[, lfc_prefix] < -lfc, , drop = FALSE]
  }
  
  
  if(nrow(x_sort) == 0){
    
    caption <- paste0("There are no significantly ", direction, "-regulated genes (", adjp_prefix, " < ", pval, ", ", lfc_prefix, " > ", lfc, ") when testing for ", contrast, ".")
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
    return(BclassDE(caption = caption))
    
    
  }else{
    
    res <- x_sort[1:(min(nrow(x_sort), topn)), , drop = FALSE]
    rownames(res) <- NULL
    
    out <- res %>% 
      mutate_at(lfc_prefix, format_difference) %>% 
      mutate_at(pval_prefix, format_pvalues) %>% 
      mutate_at(adjp_prefix, format_pvalues)
    
    if(!is.null(stats_prefixes)){
      out <- mutate_at(out, stats_prefixes, format_difference)
    }
    
  }
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("List of significantly ", direction, "-regulated genes (", adjp_prefix, " < ", pval, ", ", lfc_prefix, " > ", lfc, ") when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " genes.")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}











































