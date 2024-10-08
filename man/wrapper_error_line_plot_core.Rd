% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrappers_error_line_plot.R
\name{wrapper_error_line_plot_core}
\alias{wrapper_error_line_plot_core}
\alias{wrapper_error_line_plot_core_strat}
\title{Error line plot}
\usage{
wrapper_error_line_plot_core(
  data = NULL,
  x_var,
  y_var,
  color_line_var,
  colors_line = NULL,
  variable_names = NULL,
  title = TRUE,
  subtitle = TRUE,
  xlab = TRUE,
  ylab = TRUE,
  ggdata = NULL
)

wrapper_error_line_plot_core_strat(
  data,
  x_var,
  y_var,
  color_line_var,
  strat1_var = NULL,
  strat2_var = NULL,
  colors_line = NULL,
  variable_names = NULL,
  title = TRUE,
  xlab = TRUE,
  ylab = TRUE,
  strat1_label_both = TRUE,
  strat2_label_both = TRUE,
  legend_colors_line_title = TRUE,
  legend_position = "right",
  xlim = NULL,
  ylim = NULL,
  strat_scales = "fixed",
  strat1_nrow = 1,
  strat1_ncol = NULL,
  strat2_nrow = NULL,
  strat2_ncol = 1,
  less_legends = FALSE
)
}
\arguments{
\item{data}{Data frame.}

\item{strat1_var}{Name of the first stratification variable.}

\item{strat2_var}{Name of the second stratification variable.}
}
\description{
Error line plot
}
