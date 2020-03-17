



##############################################################################
# ORA
##############################################################################





wrapper_core_ora <- function(genes, genesets, universe, genesets_extra_info = NULL, gene_mapping = NULL, 
  method = "hypergeometric", min_GS_size = 10, max_GS_size = 500, display_topn = 10){
  
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
  
  ## We want to keep the order that is in genes so it has to be the first argument in intersect
  genes_intersection <- lapply(genesets, function(x){intersect(genes, x)})
  
  ## DE in Set
  k <- sapply(genes_intersection, length)
  
  
  if(method == "hypergeometric"){
    
    ## -0.5 like in EGSEA ORA to avoid p-value = 0
    pvalues <- phyper(q = k - 0.5, m = K, n = N - K, k = n, lower.tail = FALSE)
    
    
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
  
  if(!is.null(genesets_extra_info)){
    
    colnames(genesets_extra_info)[1] <- "Geneset"
    
    out <- left_join(out, genesets_extra_info, by = "Geneset")
    
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
# 
# universe = NULL; genesets_extra_info = NULL; gene_mapping = NULL;
# method = "hypergeometric"; min_GS_size = 10; max_GS_size = 500;
# min_DE_size = 5; topn = Inf; pval = 0.05; lfc = 0;
# gene_var = "GeneID"; lfc_prefix = "logFC"; pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val"; sep = "_";
# display_topn = 10






#' Over-representation analysis (ORA)
#' 
#' @param x TopTable.
wrapper_ora <- function(x, genesets, universe = NULL, genesets_extra_info = NULL, gene_mapping = NULL, 
  method = "hypergeometric", min_GS_size = 10, max_GS_size = 500,
  directions = c("up", "down", "both"), min_DE_size = 5, topn = Inf, pval = 0.05, lfc = 0,
  gene_var = "EntrezIDs", lfc_prefix = "logFC", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", sep = "_", 
  display_topn = 10){
  
  
  stopifnot(all(directions %in% c("up", "down", "both")))
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  ### Define the universe
  if(is.null(universe)){
    universe <- unique(x[, gene_var])
  }
  
  geneset_vars <- "Geneset"
  
  if(!is.null(genesets_extra_info)){
    colnames(genesets_extra_info)[1] <- "Geneset"
    geneset_vars <- colnames(genesets_extra_info)
  }
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", lfc_prefix, sep), "", grep(paste0("^", lfc_prefix, sep), colnames(x), value = TRUE))
  
  
  # -------------------------------------------------------------------------
  # Run ORA for each contrast up- and down-regulated genes
  # -------------------------------------------------------------------------
  
  
  res_ora <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    res_ora <- lapply(1:length(directions), function(j){
      # i = 2; j = 2
      
      # print(i)
      # print(j)
      
      contrast <- contrasts[i]
      direction <- directions[j]
      
      
      x_sign <- wrapper_dispaly_significant_genes(x, contrast = contrast, direction = direction, topn = topn, pval = pval, lfc = lfc, gene_vars = gene_var, lfc_prefix = lfc_prefix, pval_prefix = pval_prefix, adjp_prefix = adjp_prefix, sep = sep)
      
      
      if(nrow(bresults(x_sign)) < min_DE_size){
        return(NULL)
      }
      
      genes <- bresults(x_sign)[, gene_var]
      genes <- genes[!is.na(genes)]
      
      
      ### Run ORA
      
      res_ora <- wrapper_core_ora(genes, genesets = genesets, universe = universe, genesets_extra_info = genesets_extra_info, gene_mapping = gene_mapping, 
        method = method, min_GS_size = min_GS_size, max_GS_size = max_GS_size, display_topn = display_topn)
      
      
      ### Add contrast info to the column names
      
      colnames2change <- !colnames(res_ora) %in% geneset_vars
      
      colnames(res_ora)[colnames2change] <- paste(colnames(res_ora)[colnames2change], contrast, direction, sep = sep)
      
      
      return(res_ora)
      
      
    })
    
    

    null_res_ora <- sapply(res_ora, is.null)
    if(all(null_res_ora)){
      return(NULL)
    }
    res_ora <- res_ora[!null_res_ora]
    
    
    res_ora <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_ora)
    
    
  })
  
  
  null_res_ora <- sapply(res_ora, is.null)
  if(all(null_res_ora)){
    return(NULL)
  }
  res_ora <- res_ora[!null_res_ora]
  
  
  res_ora <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_ora)
  
  
}








wrapper_dispaly_significant_ora <- function(x, contrast, direction = "up", 
  topn = 20, pval = 0.05,
  geneset_vars = "Geneset", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = c("Genes", "GeneRatio", "BgRatio", "Observed", "Expected"), sep = "_", 
  caption = NULL){
  
  
  # -------------------------------------------------------------------------
  # Checks
  # -------------------------------------------------------------------------
  
  stopifnot(length(geneset_vars) >= 1)
  stopifnot(all(geneset_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  stopifnot(length(direction) == 1)
  stopifnot(direction %in% c("up", "down", "both"))
  
  if(direction == "both"){
    direction_print <- "up-, down-"
  }else{
    direction_print <- paste0(direction, "-")  
  }
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts_and_directions <- gsub(paste0("^", pval_prefix, sep), "", grep(paste0("^", pval_prefix, sep), colnames(x), value = TRUE))
  ## We add '$' because we want to match expression at the end of the string
  contrasts <- gsub(paste0(sep, direction, "$"), "", grep(paste0(sep, direction, "$"), contrasts_and_directions, value = TRUE))
  
  
  if(!contrast %in% contrasts){
    
    caption <- paste0("There was no or too few ", direction_print, "regulated genes when testing for ", contrast, ", and the ORA analysis was not possible.")
    
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
    
    caption <- paste0("There are no over-represented gene sets (", adjp_prefix, " < ", pval, ") by ", direction_print, "regulated genes when testing for ", contrast, ".")
    
    ## Remove all undescores from the caption because they are problematic when rendering to PDF
    caption <- gsub("_", " ", caption)
    
    return(BclassDE(caption = caption))
    
  }
  
  
  res <- x_sort[1:(min(nrow(x_sort), topn)), , drop = FALSE]
  rownames(res) <- NULL
  
  out <- res %>% 
    mutate_at(pval_prefix, format_pvalues2) %>% 
    mutate_at(adjp_prefix, format_pvalues2)
  
  
  # --------------------------------------------------------------------------
  # Generate caption
  # --------------------------------------------------------------------------
  
  
  if(is.null(caption)){
    
    caption <- paste0("List of over-represented gene sets (", adjp_prefix, " < ", pval, ") by ", direction_print, "regulated genes when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " gene sets.")
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}





























