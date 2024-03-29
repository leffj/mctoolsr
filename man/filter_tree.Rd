% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree-based.R
\name{filter_tree}
\alias{filter_tree}
\title{Filter tips in phylogenetic tree}
\usage{
filter_tree(tree, tip_labels)
}
\arguments{
\item{tree}{The phylogenetic tree from \code{\link{load_tree}}.}

\item{tip_labels}{The tip labels to retain.}
}
\description{
A quick way to keep certain tips in a tree and remove the others
}
\examples{
\dontrun{
filter_tree(loaded_tree, taxa_to_keep)
}
}
\concept{Phylogeny}
