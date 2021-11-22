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
# min_GS_size = 10; max_GS_size = 500; topn_genes = 20; statistic_name = "t"




#' Run GSEA
#' 
#' @param statistic Named vector of t statistics or logFC from limma.
#' @export
wrapper_gsea_core <- function(statistic, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
  name = "", sep = "_",
  min_GS_size = 10, max_GS_size = 500, topn_genes = 20, statistic_name = "t", scoreType = NULL){
  
  
  # -------------------------------------------------------------------------
  # Checks
  # -------------------------------------------------------------------------
  
  stopifnot(min_GS_size >= 1)
  stopifnot(max_GS_size > min_GS_size)
  
  statistic <- statistic[!is.na(statistic)]
  
  ## fgsea sorts the statistic in the descending order
  statistic <- sort(statistic, decreasing = TRUE)
  
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
  
  ## We want to keep the order that is in the universe so it has to be the first argument in intersect
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
  
  camera_index <- limma::ids2indices(genesets, universe)
  
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
  
  if(is.null(scoreType)){
    scoreType <- "std"
    
    if(all(statistic >= 0)){
      scoreType <- "pos"
    }else if(all(statistic <= 0)){
      scoreType <- "neg"
    }
  }
  
  ## fgsea sorts the statistic in the descending order
  fgsea_out <- fgsea::fgseaMultilevel(pathways = genesets, stats = statistic, eps = 0, scoreType = scoreType)
  
  
  # out <- fgsea_out[order(fgsea_out$pval, decreasing = FALSE), ]
  # plotEnrichment(pathway = genesets[[out$pathway[4]]], stats = statistic)
  
  # table(unlist(lapply(out$leadingEdge, length)))
  
  
  ### Make the same order as in genesets
  fgsea_out <- fgsea_out[match(names(genesets), fgsea_out$pathway), , drop = FALSE]
  
  
  stopifnot(all(fgsea_out$pathway == names(genesets)))
  
  
  # -------------------------------------------------------------------------
  # Prepare output
  # -------------------------------------------------------------------------
  
  
  out <- data.frame(GenesetID = names(genesets), stringsAsFactors = FALSE)
  
  if(!is.null(genesets_extra_info)){
    
    colnames(genesets_extra_info)[1] <- "GenesetID"
    
    out <- dplyr::left_join(out, genesets_extra_info, by = "GenesetID")
    
  }
  
  out[, paste0("Size", sep, name)] <- fgsea_out$size
  
  ### Names of leading genes
  out[, paste0("Leading.Genes", sep, name)] <- sapply(seq_along(camera_index), function(i){
    # i = 1
    
    x <- camera_index[[i]]
    
    direction <- ifelse(fgsea_out$ES[i] > 0, "Up", "Down")
    
    if(direction == "Down"){
      x <- rev(x)
    }
    
    
    ### Keep the leading genes
    
    x <- x[seq_along(fgsea_out$leadingEdge[[i]])]
    
    
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
# statistic_prefix = "t"; sep = "_";
# min_GS_size = 10; max_GS_size = 500;
# topn_genes = 10



#' Run GSEA
#' 
#' @param x TopTable
#' @export
wrapper_gsea <- function(x, genesets, genesets_extra_info = NULL, gene_mapping = NULL, 
  gene_var = "EntrezIDs", statistic_prefix = "t", sep = "_", 
  min_GS_size = 10, max_GS_size = 500,
  topn_genes = 20, scoreType = NULL){
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  geneset_vars <- "GenesetID"
  
  if(!is.null(genesets_extra_info)){
    colnames(genesets_extra_info)[1] <- "GenesetID"
    geneset_vars <- colnames(genesets_extra_info)
  }
  
  
  # -------------------------------------------------------------------------
  # Extract contrasts
  # -------------------------------------------------------------------------
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", statistic_prefix, sep), "", grep(paste0("^", statistic_prefix, sep), colnames(x), value = TRUE))
  
  
  if(length(contrasts) == 0){
    contrasts <- ""
    sep <- ""
  }
  
  
  # -------------------------------------------------------------------------
  # Run GSEA for each contrast
  # -------------------------------------------------------------------------
  
  
  res_gsea <- lapply(1:length(contrasts), function(i){
    # i = 1
    
    contrast <- contrasts[i]
    
    statistic <- x[, paste0(statistic_prefix, sep, contrast)]
    names(statistic) <- x[, gene_var]
    
    name <- contrast
    
    res_gsea <- wrapper_gsea_core(statistic = statistic, genesets = genesets, genesets_extra_info = genesets_extra_info, gene_mapping = gene_mapping, 
      name = name, sep = sep,
      min_GS_size = min_GS_size, max_GS_size = max_GS_size, topn_genes = topn_genes, statistic_name = statistic_prefix, scoreType = scoreType)
    
    
    return(res_gsea)
    
  })
  
  res_gsea <- Reduce(function(...) merge(..., by = geneset_vars, all = TRUE, sort = FALSE), res_gsea)
  
  return(res_gsea)
  
}






#' Display significantly enriched gene sets
#' 
#' @param x Data frame with GSEA results.
#' @export
wrapper_dispaly_significant_gsea <- function(x, contrast, direction = "up", 
  sort_by = "pval", topn = 20, pval = 0.05, 
  geneset_vars = "GenesetID", direction_prefix = "Direction", pval_prefix = "P.Value", adjp_prefix = "adj.P.Val", 
  stats_prefixes = c("Size", "Leading.Genes", "Median.t", "NES"), sep = "_", 
  caption = NULL){
  
  
  # -------------------------------------------------------------------------
  # Checks
  # -------------------------------------------------------------------------
  
  stopifnot(length(geneset_vars) >= 1)
  stopifnot(all(geneset_vars %in% colnames(x)))
  
  stopifnot(topn > 1)
  
  stopifnot(sort_by %in% c("none", "pval"))
  
  stopifnot(length(direction) == 1)
  stopifnot(direction %in% c("up", "down", "both"))
  
  if(direction == "both"){
    direction_print <- "up-, down-"
  }else{
    direction_print <- paste0(direction, "-")  
  }
  
  
  if(contrast == ""){
    sep <- ""
  }
  
  
  # -------------------------------------------------------------------------
  # Preprocessing
  # -------------------------------------------------------------------------
  
  
  ## We add '^' because we want to match expression at the beginning of the string
  contrasts <- gsub(paste0("^", direction_prefix, sep), "", grep(paste0("^", direction_prefix, sep), colnames(x), value = TRUE))
  
  
  stopifnot(contrast %in% contrasts)
  
  
  ## Find columns corresponding to the contrast and subset the data
  ## We add '$' because we want to match expression at the end of the string
  
  contrast_vars_display <- paste0(c(direction_prefix, stats_prefixes, pval_prefix, adjp_prefix), sep, contrast)
  
  x <- x[ , c(geneset_vars, contrast_vars_display), drop = FALSE]
  
  colnames(x) <- gsub(paste0(sep, contrast, "$"), "", colnames(x))
  
  
  ## Sort
  if(sort_by == "pval"){
    x_sort <- x[order(x[, pval_prefix], decreasing = FALSE), , drop = FALSE]
  }else{
    x_sort <- x
  }
  
  ## Subset by adj. p-value
  x_sort <- x_sort[x_sort[, adjp_prefix] <= pval & !is.na(x_sort[, adjp_prefix]), , drop = FALSE]
  
  
  ## Subset by direction
  if(direction == "up"){
    x_sort <- x_sort[x_sort[, direction_prefix] == "Up", , drop = FALSE]
  }else if(direction == "down"){
    x_sort <- x_sort[x_sort[, direction_prefix] == "Down", , drop = FALSE]
  }
  
  
  if(nrow(x_sort) == 0){
    
    caption <- paste0("There are no enriched gene sets (", adjp_prefix, " <= ", pval, ") by ", direction_print, "regulated genes", ifelse(contrast == "", ".", paste0(" when comparing ", contrast, ".")))
    
    ## Remove all underscores from the caption because they are problematic when rendering to PDF
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
    
    caption <- paste0("List of gene sets enriched (", adjp_prefix, " <= ", pval, ") by ", direction_print, "regulated genes", ifelse(contrast == "", ".", paste0(" when comparing ", contrast, ".")))
    
    
    if(nrow(x_sort) >  topn){
      
      caption <- paste0(caption, " Printed ", topn, " out of ", nrow(x_sort), " gene sets.")
      
    }
    
  }
  
  ## Remove all underscores from the caption because they are problematic when rendering to PDF
  caption <- gsub("_", " ", caption)
  
  rownames(res) <- NULL
  rownames(out) <- NULL
  
  bout <- BclassDE(results = res, output = out, caption = caption)
  
  
  return(bout)
  
  
  
  
  
}



































