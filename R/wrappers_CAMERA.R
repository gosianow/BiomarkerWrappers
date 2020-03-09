





##############################################################################
# CAMERA
##############################################################################



# statistic <- topTable$t_LUMB_VS_LUMA
# names(statistic) <- topTable$GeneID
# 
# genesets <- geneset_list_reactome
# genesets_extra_info <- geneset_extra_reactome 
# gene_mapping <- entrez2hgnc
# 
# 
# min_GS_size = 10; max_GS_size = 500; display_topn = 20; statistic_name = "t"




#' Run CAMERA
#' 
#' @param statistic Named vector of t statistics from limma or logFC.
wrapper_core_cameraPR <- function(statistic, genesets, genesets_extra_info = NULL, gene_mapping = NULL, min_GS_size = 10, max_GS_size = 500, display_topn = 10, statistic_name = "t"){
  
  
  # -------------------------------------------------------------------------
  # Checks
  # -------------------------------------------------------------------------
  
  stopifnot(min_GS_size >= 1)
  stopifnot(max_GS_size > min_GS_size)
  
  statistic <- sort(statistic, decreasing = FALSE)
  
  stopifnot(!is.null(names(statistic)))
  
  universe <- names(statistic)
  
  if(!is.null(gene_mapping)){
    stopifnot(all(intersect(unlist(genesets), universe) %in% gene_mapping[, 1]))
  }
  
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
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
  
  ### Generate index for camera
  
  camera_index <- ids2indices(genesets, universe)
  
  # camera_index[["DNA Replication"]]
  
  if(!is.null(gene_mapping)){
    
    mm <- match(universe, gene_mapping[, 1])
    mapped_universe <- gene_mapping[mm, 2]
    
  }else{
    mapped_universe <- universe
  }
  
  
  # -------------------------------------------------------------------------
  # Run cameraPR
  # -------------------------------------------------------------------------
  
  
  camara_out <- limma::cameraPR(statistic = statistic, index = camera_index, inter.gene.cor = 0.05, sort = FALSE)
  
  stopifnot(all(rownames(camara_out) == names(genesets)))
  
  
  # -------------------------------------------------------------------------
  # Prepare output
  # -------------------------------------------------------------------------
  
  
  out <- data.frame(Geneset = names(genesets), stringsAsFactors = FALSE)
  
  if(!is.null(genesets_extra_info)){
    
    colnames(genesets_extra_info)[1] <- "Geneset"
    
    out <- left_join(out, genesets_extra_info, by = "Geneset")
    
  }
  
  out$NGenes <- camara_out$NGenes
  
  out$Genes <- sapply(seq_along(camera_index), function(i){
    # i = 1
    
    x <- camera_index[[i]]
    
    direction <- camara_out$Direction[i]
    
    if(direction == "Up"){
      x <- rev(x)
    }
    
    suffix <- ""
    
    if(length(x) > display_topn){
      x <- x[seq_len(display_topn)]
      suffix <- ", ..."
    }
    
    paste0(paste0(mapped_universe[x], collapse = ", "), suffix)
    
  })
  
  
  out[, paste0("Mean.", statistic_name)] <- round(sapply(camera_index, function(x){
    mean(statistic[x], na.rm = TRUE)
  }), 2)
  
  out$Direction <- camara_out$Direction
  
  out$P.Value <- camara_out$PValue
  
  out$adj.P.Val <- camara_out$FDR
  
  
  ### Sort by p-value
  
  out <- out[order(out$P.Value, decreasing = FALSE), , drop = FALSE]
  
  out
  
  
  
}









# x <- topTable
# 
# genesets <- geneset_list_reactome
# genesets_extra_info <- geneset_extra_reactome 
# gene_mapping <- entrez2hgnc
# 
# 
# min_GS_size = 10; max_GS_size = 500;
# gene_var = "GeneID"; statistic_prefix = "t"; sep = "_"; 
# display_topn = 10



#' Run CAMERA
#' 
#' @param x TopTable
wrapper_cameraPR <- function(x, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
  min_GS_size = 10, max_GS_size = 500,
  gene_var = "EntrezIDs", statistic_prefix = "t", sep = "_", 
  display_topn = 10){
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  geneset_vars <- "Geneset"
  
  if(!is.null(genesets_extra_info)){
    colnames(genesets_extra_info)[1] <- "Geneset"
    geneset_vars <- colnames(genesets_extra_info)
  }
  
  
  # -------------------------------------------------------------------------
  # Run CAMERA for each contrast
  # -------------------------------------------------------------------------
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", statistic_prefix, sep), "", grep(paste0("^", statistic_prefix, sep), colnames(x), value = TRUE))
  
  
  res_camera <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    contrast <- contrasts[i]
    
    statistic <- x[, paste0(statistic_prefix, sep, contrast)]
    names(statistic) <- x[, gene_var]
    
    
    res_camera <- wrapper_core_cameraPR(statistic = statistic, genesets = genesets, genesets_extra_info = genesets_extra_info, gene_mapping = gene_mapping, 
      min_GS_size = min_GS_size, max_GS_size = max_GS_size, display_topn = display_topn, statistic_name = statistic_prefix)
    
    
    ### Add contrast info to the column names
    
    colnames2change <- !colnames(res_camera) %in% geneset_vars
    
    colnames(res_camera)[colnames2change] <- paste(colnames(res_camera)[colnames2change], contrast, sep = sep)
    
    return(res_camera)
    
  })
  
  res_camera <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_camera)
  
  return(res_camera)
  
}






# x <- topTable_camera
# 
# direction = "up"
# topn = 20; pval = 0.05; 
# geneset_vars = "Geneset"; direction_prefix = "Direction"; pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val"; 
# stats_prefixes = c("NGenes", "Genes", "Mean.t"); sep = "_"; 
# caption = NULL


wrapper_dispaly_significant_camera <- function(x, contrast, direction = "up", 
  topn = 20, pval = 0.05, 
  geneset_vars = "Geneset", direction_prefix = "Direction", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = c("NGenes", "Genes", "Mean.t"), sep = "_", 
  caption = NULL){
  
  
  stopifnot(length(geneset_vars) >= 1)
  stopifnot(all(geneset_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  stopifnot(direction %in% c("up", "down", "both"))
  
  if(direction == "both"){
    direction_print <- "up-, down-"
  }else{
    direction_print <- paste0(direction, "-")  
  }
  
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", direction_prefix, sep), "", grep(paste0("^", direction_prefix, sep), colnames(x), value = TRUE))
  
  stopifnot(contrast %in% contrasts)
  
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(direction_prefix, stats_prefixes, pval_prefix, adjp_prefix), sep, contrast)
  
  x <- x[ , c(geneset_vars, contrast_vars_display), drop = FALSE]
  
  colnames(x) <- gsub(paste0(sep, contrast, "$"), "", colnames(x))
  
  ## Sort by p-value
  x_sort <- x[order(x[, pval_prefix], decreasing = FALSE), , drop = FALSE]
  ## Subset by adj. p-value
  x_sort <- x_sort[x_sort[, adjp_prefix] < pval, , drop = FALSE]
  
  
  ## Subset by LogFC
  if(direction == "up"){
    x_sort <- x_sort[x_sort[, direction_prefix] == "Up", , drop = FALSE]
  }else if(direction == "down"){
    x_sort <- x_sort[x_sort[, direction_prefix] == "Down", , drop = FALSE]
  }
  
  
  if(nrow(x_sort) == 0){
    
    caption <- paste0("There are no enriched gene sets (", adjp_prefix, " < ", pval, ") by ", direction_print, "regulated genes when testing for ", contrast, ".")
    
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
    
    caption <- paste0("List of enriched gene sets (", adjp_prefix, " < ", pval, ") by ", direction_print, "regulated genes when testing for ", contrast, ".")
    
    if(nrow(x_sort) >  topn){
      
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " gene sets.")
      
    }
    
  }
  
  ## Remove all undescores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
}



































