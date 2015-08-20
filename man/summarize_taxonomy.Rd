% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/taxonomy_based.R
\name{summarize_taxonomy}
\alias{summarize_taxonomy}
\title{Calculate Values for Coarser Taxonomic Groupings}
\usage{
summarize_taxonomy(input, level, relative = TRUE, report_higher_tax = TRUE)
}
\arguments{
\item{input}{The input dataset}

\item{level}{The taxonomy level to summarize by [integer]}

\item{relative}{Return relative abundances or not (default = TRUE)}

\item{report_higher_tax}{Whether or not to return taxonomic strings higher
than the indecated level (default = TRUE)}
}
\description{
Given input as generated from \code{load_taxon_table()},
  calculate relative abundances (or absolute abundances) of taxa at a coarser
  taxonomic threshold. Taxonomy must have been loaded for this function to
  work.
}
\examples{
## Return relative abundances of phyla for each sample
\dontrun{
summarize_taxonomy(input = example_input, level = 2)
}
}
