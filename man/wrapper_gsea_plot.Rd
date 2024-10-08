% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers_GSEA_plot.R
\name{wrapper_gsea_plot}
\alias{wrapper_gsea_plot}
\title{Plot GSEA statistics ranks}
\usage{
wrapper_gsea_plot(
  x,
  contrast,
  genesets,
  gene_var = "EntrezIDs",
  statistic_prefix = "t",
  sep = "_",
  gsea_results = NULL,
  geneset_var = "GenesetID",
  adjp_var = "adj.P.Val",
  enrichment_score_var = "NES",
  trim_limits = 0.01,
  color_low = "#42399B",
  color_mid = "darkgrey",
  color_high = "#D70131",
  title = "",
  title_size = 10,
  title_width = 100,
  axis_text_y_size = 9,
  axis_text_y_width = 80
)
}
\arguments{
\item{x}{TopTable with statistic, for example, DGE results from limma.}

\item{gsea_results}{TopTable with selected GSEA results obtained by running 'wrapper_dispaly_significant_gsea'.}
}
\description{
Plot GSEA statistics ranks
}
