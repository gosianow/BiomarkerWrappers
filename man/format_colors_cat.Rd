% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{format_colors_cat}
\alias{format_colors_cat}
\alias{format_colors}
\alias{format_colors_cat_strata}
\title{Format or create colors for a vector with categorical values}
\usage{
format_colors_cat(
  x,
  colors = NULL,
  palette = NULL,
  rev = FALSE,
  allow_duplicated = TRUE
)

format_colors(
  x,
  colors = NULL,
  palette = NULL,
  rev = FALSE,
  allow_duplicated = TRUE
)

format_colors_cat_strata(x, strata = NULL, palette = NULL, rev = FALSE)
}
\arguments{
\item{x}{Vector of categorical values for which we want to specify colors.}

\item{colors}{Vector of colors longer or equal the number of unique levels of x. Can be named or non-named. If NULL, colors are generated.}

\item{palette}{List of length corresponding the number of levels in strata with colors or names of RColorBrewer palettes. E.g. `list("Blues", "Oranges", "Greens", "Purples", "Reds", "Greys")`}

\item{strata}{Vector of categorical values of stratification groups.}
}
\value{
Named vector of unique colors for all unique values of x.
}
\description{
Make sure that unique named colors are created for each value.
}
\examples{

x <- c("low", "high")

colors <- c("high" = "red", "low" = "grey")
palette <- NULL
allow_duplicated <- FALSE


format_colors_cat(x, colors = colors, palette = palette, allow_duplicated = allow_duplicated)


x <- c("<1", "<2", "<3", "<4")

colors <- NULL
palette <- "RdYlBu"
allow_duplicated <- FALSE


out <- format_colors_cat(x, colors = colors, palette = palette, allow_duplicated = allow_duplicated)

barplot(rep(1, length(out)), col = out)


}
