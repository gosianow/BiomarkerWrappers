% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers_GSEA.R
\name{wrapper_gsea_core}
\alias{wrapper_gsea_core}
\title{Run GSEA}
\usage{
wrapper_gsea_core(
  statistic,
  genesets,
  genesets_extra_info = NULL,
  gene_mapping = NULL,
  name = "",
  sep = "_",
  min_GS_size = 10,
  max_GS_size = 500,
  topn_genes = 20,
  statistic_name = "t",
  scoreType = NULL,
  nPermSimple = 1000
)
}
\arguments{
\item{statistic}{Named vector of t statistics or logFC from limma.}
}
\description{
Run GSEA with fgsea::fgseaMultilevel
}
