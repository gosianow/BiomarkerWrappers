





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
# min_GS_size = 10; max_GS_size = 500; topn_genes = 20; statistic_name = "t"




#' Run CAMERA
#' 
#' @param statistic Named vector of t statistics from limma or logFC.
#' @export
wrapper_cameraPR_core <- function(statistic, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
  name = "", sep = "_",
  min_GS_size = 10, max_GS_size = 500, topn_genes = 20, statistic_name = "t"){
  
  
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
  
  if(name == ""){
    sep = ""
  }
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  genesets <- lapply(genesets, unique)
  
  ### Keep genes in the gene sets that are in the universe
  
  ## We want to keep the order that is in the univarse so it has to be the first argument in intersect
  genesets <- lapply(genesets, function(x){intersect(universe, x)})
  
  ### Exclude too small or too large gene sets
  
  size_genesets <- sapply(genesets, length)
  
  genesets <- genesets[size_genesets >= min_GS_size & size_genesets <= max_GS_size]
  
  if(length(genesets) == 0){
    message("There are no common genes between the universe and the genesets.")
    return(NULL)
  }
  
  ### Generate index for camera
  
  camera_index <- limma::ids2indices(genesets, universe)
  
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
  
  
  camera_out <- limma::cameraPR(statistic = statistic, index = camera_index, inter.gene.cor = 0.05, sort = FALSE)
  
  stopifnot(all(rownames(camera_out) == names(genesets)))
  
  
  # -------------------------------------------------------------------------
  # Prepare output
  # -------------------------------------------------------------------------
  
  
  out <- data.frame(GenesetID = names(genesets), stringsAsFactors = FALSE)
  
  if(!is.null(genesets_extra_info)){
    
    colnames(genesets_extra_info)[1] <- "GenesetID"
    
    out <- dplyr::left_join(out, genesets_extra_info, by = "GenesetID")
    
  }
  
  out[, paste0("Size", sep, name)] <- camera_out$NGenes
  
  ### Names of leading genes
  out[, paste0("Leading.Genes", sep, name)] <- sapply(seq_along(camera_index), function(i){
    # i = 1
    
    x <- camera_index[[i]]
    
    direction <- camera_out$Direction[i]
    
    if(direction == "Up"){
      x <- rev(x)
    }
    
    suffix <- ""
    
    if(length(x) > topn_genes){
      x <- x[seq_len(topn_genes)]
      suffix <- ", ..."
    }
    
    paste0(paste0(mapped_universe[x], collapse = ", "), suffix)
    
  })
  
  
  out[, paste0("Median.", statistic_name, sep, name)] <- round(sapply(camera_index, function(x){
    median(statistic[x], na.rm = TRUE)
  }), 2)
  
  out[, paste0("Direction", sep, name)] <- camera_out$Direction
  
  out[, paste0("P.Value", sep, name)] <- camera_out$PValue
  
  out[, paste0("adj.P.Val", sep, name)] <- camera_out$FDR
  
  
  ### Sort by p-value
  
  out <- out[order(out[, paste0("P.Value", sep, name)], decreasing = FALSE), , drop = FALSE]
  
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
# topn_genes = 10



#' Run CAMERA
#' 
#' @param x TopTable
#' @export
wrapper_cameraPR <- function(x, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
  min_GS_size = 10, max_GS_size = 500,
  gene_var = "EntrezIDs", statistic_prefix = "t", sep = "_", 
  topn_genes = 20){
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  geneset_vars <- "GenesetID"
  
  if(!is.null(genesets_extra_info)){
    colnames(genesets_extra_info)[1] <- "GenesetID"
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
    
    name <- contrast

    res_camera <- wrapper_cameraPR_core(statistic = statistic, genesets = genesets, genesets_extra_info = genesets_extra_info, gene_mapping = gene_mapping, 
      name = name, sep = sep,
      min_GS_size = min_GS_size, max_GS_size = max_GS_size, topn_genes = topn_genes, statistic_name = statistic_prefix)
    

    return(res_camera)
    
  })
  
  res_camera <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_camera)
  
  return(res_camera)
  
}






# x <- topTable_camera
# 
# direction = "up"
# topn = 20; pval = 0.05; 
# geneset_vars = "GenesetID"; direction_prefix = "Direction"; pval_prefix = "P.Value"; adjp_prefix = "adj.P.Val"; 
# stats_prefixes = c("Size", "Leading.Genes", "Mean.t"); sep = "_"; 
# caption = NULL



#' Display significantly enriched gene sets
#' 
#' @param x Data frame with CAMERA results.
#' @export
wrapper_dispaly_significant_camera <- function(x, contrast, direction = "up", 
  sort_by = "pval", topn = 20, pval = 0.05, 
  geneset_vars = "GenesetID", direction_prefix = "Direction", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = c("Size", "Leading.Genes", "Median.t"), sep = "_", 
  caption = NULL){
  
  
  out <- wrapper_dispaly_significant_gsea(x, contrast = contrast, direction = direction, 
    sort_by = sort_by, topn = topn, pval = pval, 
    geneset_vars = geneset_vars, direction_prefix = direction_prefix, 
    pval_prefix = pval_prefix, adjp_prefix = adjp_prefix, 
    stats_prefixes = stats_prefixes, sep = sep, 
    caption = caption)
  
  
  return(out)
  
  
}



































