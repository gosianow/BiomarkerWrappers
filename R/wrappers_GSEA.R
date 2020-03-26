#' @include wrappers_CAMERA.R 
NULL





##############################################################################
# GSEA
##############################################################################



# statistic <- topTable$t_ECR1_VS_parental_NoTrt
# names(statistic) <- topTable$EntrezIDs
# 
# genesets <- geneset_list_reactome
# genesets_extra_info <- geneset_extra_reactome
# gene_mapping <- entrez2hgnc
# 
# 
# min_GS_size = 10; max_GS_size = 500; display_topn = 20; statistic_name = "t"




#' Run GSEA
#' 
#' @param statistic Named vector of t statistics from limma or logFC.
wrapper_core_gsea <- function(statistic, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
  name = "", sep = "_",
  min_GS_size = 10, max_GS_size = 500, display_topn = 10, statistic_name = "t"){
  
  
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
  
  # genesets[[1]]
  # statistic[genesets[[1]]]
  
  
  ### Keep genes in the gene sets that are in the universe
  
  ## We want to keep the order that is in the univarse so it has to be the first argument in intersect
  genesets <- lapply(genesets, function(x){intersect(universe, x)})
  
  # genesets[[1]]
  # statistic[genesets[[1]]]
  
  
  ### Exclude too small or too large gene sets
  
  size_genesets <- sapply(genesets, length)
  
  genesets <- genesets[size_genesets >= min_GS_size & size_genesets <= max_GS_size]
  
  if(length(genesets) == 0){
    message("There are no common genes between the universe and the genesets.")
    return(NULL)
  }
  
  ### Generate index for camera
  
  camera_index <- ids2indices(genesets, universe)
  
  # camera_index[[1]]
  
  if(!is.null(gene_mapping)){
    
    mm <- match(universe, gene_mapping[, 1])
    mapped_universe <- gene_mapping[mm, 2]
    
  }else{
    mapped_universe <- universe
  }
  
  
  # -------------------------------------------------------------------------
  # Run fgseaMultilevel
  # -------------------------------------------------------------------------
  
  ## Actually, fgsea sorts the statistic in the descending order
  fgsea_out <- fgsea::fgseaMultilevel(pathways = genesets, stats = statistic)
  
  ### Make the same order as in genesets
  fgsea_out <- fgsea_out[match(names(genesets), fgsea_out$pathway), , drop = FALSE]
  
  stopifnot(all(fgsea_out$pathway == names(genesets)))
  
  
  # -------------------------------------------------------------------------
  # Prepare output
  # -------------------------------------------------------------------------
  
  
  out <- data.frame(Geneset = names(genesets), stringsAsFactors = FALSE)
  
  if(!is.null(genesets_extra_info)){
    
    colnames(genesets_extra_info)[1] <- "Geneset"
    
    out <- left_join(out, genesets_extra_info, by = "Geneset")
    
  }
  
  out[, paste0("size", sep, name)] <- fgsea_out$size
  
  ### Names of leading genes
  out[, paste0("Genes", sep, name)] <- sapply(seq_along(camera_index), function(i){
    # i = 1
    
    x <- camera_index[[i]]
    
    direction <- ifelse(fgsea_out$ES[i] > 0, "Up", "Down")
    
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
  
  
  out[, paste0("Mean.", statistic_name, sep, name)] <- round(sapply(camera_index, function(x){
    mean(statistic[x], na.rm = TRUE)
  }), 2)
  
  out[, paste0("Direction", sep, name)] <- ifelse(fgsea_out$ES > 0, "Up", "Down")
  out[, paste0("ES", sep, name)] <- round(fgsea_out$ES, 2)
  out[, paste0("NES", sep, name)] <- round(fgsea_out$NES, 2)
  
  out[, paste0("P.Value", sep, name)] <- fgsea_out$pval
  
  out[, paste0("adj.P.Val", sep, name)] <- fgsea_out$padj
  
  
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
# display_topn = 10



#' Run GSEA
#' 
#' @param x TopTable
wrapper_gsea <- function(x, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
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
  # Run GSEA for each contrast
  # -------------------------------------------------------------------------
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", statistic_prefix, sep), "", grep(paste0("^", statistic_prefix, sep), colnames(x), value = TRUE))
  
  
  res_gsea <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    contrast <- contrasts[i]
    
    statistic <- x[, paste0(statistic_prefix, sep, contrast)]
    names(statistic) <- x[, gene_var]
    
    name <- contrast

    res_gsea <- wrapper_core_gsea(statistic = statistic, genesets = genesets, genesets_extra_info = genesets_extra_info, gene_mapping = gene_mapping, 
      name = name, sep = sep,
      min_GS_size = min_GS_size, max_GS_size = max_GS_size, display_topn = display_topn, statistic_name = statistic_prefix)
    
    
    return(res_gsea)
    
  })
  
  res_gsea <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_gsea)
  
  return(res_gsea)
  
}








wrapper_dispaly_significant_gsea <- function(x, contrast, direction = "up", 
  sort_by = "pval", topn = 20, pval = 0.05, 
  geneset_vars = "Geneset", direction_prefix = "Direction", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = c("size", "Genes", "Mean.t"), sep = "_", 
  caption = NULL){
  
  
  out <- wrapper_dispaly_significant_camera(x, contrast = contrast, direction = direction, 
    sort_by = sort_by, topn = topn, pval = pval, 
    geneset_vars = geneset_vars, direction_prefix = direction_prefix, pval_prefix = pval_prefix, adjp_prefix = adjp_prefix, 
    stats_prefixes = stats_prefixes, sep = sep, 
    caption = caption)
  
  return(out)
  
}



































