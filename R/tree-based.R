#mctoolsr

## tree-based functions

#' @title Load in a phylogenetic tree
#' @description A wrapper for \code{\link[ape]{read.tree}} in
#'   \code{\link[ape]{ape}}. Reads in a phylogenetic tree in Newic or New
#'   Hampshire format.
#' @param tree_fp The phylogenetic tree file path.
#' @concept Load external data
#' @examples 
#' \dontrun{
#' load_tree("tree file path")
#' }
load_tree = function(tree_fp) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("'ape' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  tree = ape::read.tree(file = tree_fp)
  tree$tip.label = gsub("'", "", tree$tip.label)
  tree
}


#' @title Filter tips in phylogenetic tree
#' @param tree The phylogenetic tree from \code{\link{load_tree}}.
#' @param tip_labels The tip labels to retain.
#' @concept Phylogeny
#' @examples 
#' \dontrun{
#' filter_tree(loaded_tree, taxa_to_keep)
#' }
filter_tree = function(tree, tip_labels) {
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("'ape' package needed for this function to work. Please install it.",
         call. = FALSE)
  }
  tree_f = ape::drop.tip(tree, tree$tip.label[!tree$tip.label %in% tip_labels])
  tree_f
}