





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
  topn = 20, pval = 0.05, lfc = 0, 
  gene_vars = c("Hgnc_Symbol", "EntrezIDs", "GeneName"), lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = NULL, sep = "_", 
  caption = NULL){
  
  
  
  stopifnot(length(gene_vars) >= 1)
  stopifnot(all(gene_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", lfc_prefix, sep), "", grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE))
  
  stopifnot(contrast %in% contrasts)
  
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(lfc_prefix, stats_prefixes, pval_prefix, adjp_prefix), sep, contrast)
  
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
    
    caption <- paste0("There are no significantly ", direction, "-regulated genes (", adjp_prefix, " < ", pval, ", |", lfc_prefix, "| > ", lfc, ") when testing for ", contrast, ".")
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
    return(BclassDE(caption = caption))
    
  }
  
  
  
  res <- x_sort[1:(min(nrow(x_sort), topn)), , drop = FALSE]
  rownames(res) <- NULL
  
  out <- res %>% 
    mutate_at(lfc_prefix, format_difference) %>% 
    mutate_at(pval_prefix, format_pvalues) %>% 
    mutate_at(adjp_prefix, format_pvalues)
  
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("List of significantly ", direction, "-regulated genes (", adjp_prefix, " < ", pval, ", |", lfc_prefix, "| > ", lfc, ") when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " genes.")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}











wrapper_core_ora <- function(genes, genesets, universe, genesets_extra = NULL, gene_mapping = NULL, 
  method = "hypergeometric", min_GS_size = 10, max_GS_size = 500, display_topn = 20){
  
  # -------------------------------------------------------------------------
  # Checks
  # -------------------------------------------------------------------------
  
  stopifnot(min_GS_size >= 1)
  stopifnot(max_GS_size > min_GS_size)
  
  stopifnot(all(genes %in% universe))
  
  if(!is.null(gene_mapping)){
    stopifnot(all(genes %in% gene_mapping[, 1]))
  }
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  genes <- unique(genes)
  universe <- unique(universe)
  genesets <- lapply(genesets, unique)
  
  ### Keep genes in the gene sets that are in the universe
  
  genesets <- lapply(genesets, intersect, universe)
  
  ### Exclude too small or too large gene sets
  
  size_genesets <- sapply(genesets, length)
  
  genesets <- genesets[size_genesets >= min_GS_size & size_genesets <= max_GS_size]
  
  if(length(genesets) == 0){
    message("There are no common genes between the universe and the genesets.")
    return(NULL)
  }
  
  
  # -------------------------------------------------------------------------
  # Calculate values for the test
  # -------------------------------------------------------------------------
  
  ## Universe
  N <- length(universe)
  ## DE
  n <- length(genes)
  ## Set
  K <- sapply(genesets, length)
  
  genes_intersection <- lapply(genesets, function(x){intersect(genes, x)})
  
  ## DE in Set
  k <- sapply(genes_intersection, length)
  
  
  if(method == "hypergeometric"){
    
    pvalues <- phyper(q = k, m = K, n = N - K, k = n, lower.tail = FALSE)
    
    
  }else{
    
    pvalues <- sapply(1:length(genesets), function(i){
      # i = 1
      
      tbl <- data.frame(not_DE = c(K[i] - k[i], (N - K[i]) - (n - k[i])), DE = c(k[i], n - k[i]), row.names = c("in_Set", "not_in_Set"))
      
      pvalue <- fisher.test(tbl, alternative = "less")$p.value
      
    })
    
  }
  
  # -------------------------------------------------------------------------
  # Prepare output
  # -------------------------------------------------------------------------
  
  
  out <- data.frame(Geneset = names(genesets), stringsAsFactors = FALSE)
  
  if(!is.null(genesets_extra)){
    
    colnames(genesets_extra)[1] <- "Geneset"
    
    out <- left_join(out, genesets_extra, by = "Geneset")
    
  }
  
  
  out$Genes <- sapply(genes_intersection, function(x){
    
    suffix <- ""
    
    if(length(x) > display_topn){
      x <- x[seq_len(display_topn)]
      suffix <- ", ..."
    }
    
    if(!is.null(gene_mapping)){
      
      mm <- match(x, gene_mapping[, 1])
      
      paste0(paste0(gene_mapping[mm, 2], collapse = ", "), suffix)
      
    }else{
      
      paste0(paste0(x, collapse = ", "), suffix)
      
    }
    
  })
  
  
  out$GeneRatio <- paste0(k, " / ", n)
  out$BgRatio <- paste0(K, " / ", N)
  
  out$Observed <- k
  out$Expected <- round(K/N*n, 1)
  
  out$P.Value <- pvalues
  
  out$adj.P.Val <- p.adjust(pvalues, method = "BH")
  
  
  ### Sort by p-value
  
  out <- out[order(out$P.Value, decreasing = FALSE), , drop = FALSE]
  
  out
  
  
}







# x <- topTable
# genesets <- geneset_list_reactome
# 
# universe = NULL; 
# 
# genesets_extra = geneset_extra_reactome 
# gene_mapping = entrez2hgnc
# 
# 
# # universe = NULL; genesets_extra = NULL; gene_mapping = NULL; 
# 
# method = "hypergeometric"; min_GS_size = 10; max_GS_size = 500;
# 
# min_DE_size = 5; topn = Inf; pval = 0.05; lfc = 0;
# gene_var = "GeneID"; lfc_prefix = "logFC"; pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val"; sep = "_";
# display_topn = 10







#' Over-representation analysis (ORA)
#' 
#' @param x 
wrapper_ora <- function(x, genesets, universe = NULL, genesets_extra = NULL, gene_mapping = NULL, 
  method = "hypergeometric", min_GS_size = 10, max_GS_size = 500,
  min_DE_size = 5, topn = Inf, pval = 0.05, lfc = 0,
  gene_var = "EntrezIDs", lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", sep = "_", 
  display_topn = 10){
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  ### Define the universe
  if(is.null(universe)){
    universe <- unique(x[, gene_var])
  }
  
  geneset_vars <- "Geneset"
  
  if(!is.null(genesets_extra)){
    colnames(genesets_extra)[1] <- "Geneset"
    geneset_vars <- colnames(genesets_extra)
  }
  
  
  gene_vars <- gene_var
  
  
  # -------------------------------------------------------------------------
  # Run ORA for each contrast up- and down-regulated genes
  # -------------------------------------------------------------------------
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", lfc_prefix, sep), "", grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE))
  
  directions <- c("up", "down")
  
  
  res_ora <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    res_ora <- lapply(1:length(directions), function(j){
      # i = 1; j = 1
      
      contrast <- contrasts[i]
      direction <- directions[j]
      
      
      genes <- wrapper_dispaly_significant_genes(x, contrast = contrasts, direction = direction, topn = topn, pval = pval, lfc = lfc, gene_vars = gene_var, lfc_prefix = lfc_prefix, pval_prefix = pval_prefix, adjp_prefix = adjp_prefix, sep = sep)
      
      genes <- Bresults(genes)[, gene_var]
      genes <- genes[!is.na(genes)]
      
      if(length(genes) < min_DE_size){
        return(NULL)
      }
      
      ### Run ORA
      
      res_ora <- wrapper_core_ora(genes, genesets = genesets, universe = universe, genesets_extra = genesets_extra, gene_mapping = gene_mapping, 
        method = method, min_GS_size = min_GS_size, max_GS_size = max_GS_size, display_topn = display_topn)
      
      
      ### Add contrast info to the column names
      
      colnames2change <- !colnames(res_ora) %in% geneset_vars
      
      colnames(res_ora)[colnames2change] <- paste(colnames(res_ora)[colnames2change], contrast, direction, sep = sep)
      
      
      return(res_ora)
      
      
    })
    
    res_ora <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_ora)
    
    
  })
  
  
  res_ora <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_ora)
  
  return(res_ora)
  
}








wrapper_dispaly_significant_ora <- function(x, contrast, direction = "up", 
  topn = 20, pval = 0.05,
  geneset_vars = "Geneset", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = c("Genes", "GeneRatio", "BgRatio", "Observed", "Expected"), sep = "_", 
  caption = NULL){
  
  
  stopifnot(length(geneset_vars) >= 1)
  stopifnot(all(geneset_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", pval_prefix, sep), "", grep(paste0("^", pval_prefix, sep), colnames(x), value = TRUE))
  ## We add '$' because we want to match expression at the end of the string
  contrasts <- gsub(paste0(sep, direction, "$"), "", grep(paste0(sep, direction, "$"), contrasts, value = TRUE))
  
  
  if(!contrast %in% contrasts){
    
    caption <- paste0("There was no significant ", direction, "-regulated genes when testing for ", contrast, ", and ORA results are not available.")
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
    return(BclassDE(caption = caption))
    
  }
  
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(stats_prefixes, pval_prefix, adjp_prefix), sep, contrast, sep, direction)
  
  x <- x[ , c(geneset_vars, contrast_vars_display), drop = FALSE]
  
  colnames(x) <- gsub(paste0(sep, contrast, sep, direction, "$"), "", colnames(x))
  
  
  ## Sort by p-value
  x_sort <- x[order(x[, pval_prefix], decreasing = FALSE), , drop = FALSE]
  ## Subset by adj. p-value
  x_sort <- x_sort[x_sort[, adjp_prefix] < pval, , drop = FALSE]
  
  
  if(nrow(x_sort) == 0){
    
    caption <- paste0("There are no significantly over-represented gene sets (", adjp_prefix, " < ", pval, ") by ", direction, "-regulated genes when testing for ", contrast, ".")
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
    return(BclassDE(caption = caption))
    
  }
  
  
  res <- x_sort[1:(min(nrow(x_sort), topn)), , drop = FALSE]
  rownames(res) <- NULL
  
  out <- res %>% 
    mutate_at(pval_prefix, format_pvalues) %>% 
    mutate_at(adjp_prefix, format_pvalues)
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("List of significantly over-represented gene sets (", adjp_prefix, " < ", pval, ") by ", direction, "-regulated genes when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " gene sets.")
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}




















