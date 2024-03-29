% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc_analyses.R
\name{return_top_taxa}
\alias{return_top_taxa}
\title{Return the most abundant taxa in a dataset}
\usage{
return_top_taxa(input, number_taxa)
}
\arguments{
\item{input}{The input dataset as loaded by \code{\link{load_taxa_table}}.}

\item{number_taxa}{The number of top taxa to display.}
}
\description{
Return the n most abundant taxa as calculated by the mean 
sequence counts across all samples.
}
\concept{Misc analyses}
