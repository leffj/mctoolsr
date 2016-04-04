#mctoolsr

## tree-based functions

#' @title Load in a phylogenetic tree
#' @description A wrapper for \code{read.tree()} in \code{ape}. Reads in a 
#'  phylogenetic tree in Newic or New Hampshire format.
#' @param tree_fp The phylogenetic tree file path.
load_tree = function(tree_fp) {
  tree = ape::read.tree(file = tree_fp)
  tree$tip.label = gsub("'", "", tree$tip.label)
  tree
}


#' @title Filter tips in phylogenetic tree
#' @param tree The phylogenetic tree from \code{load_tree()}.
#' @param tip_labels The tip labels to retain.
filter_tree = function(tree, tip_labels) {
  tree_f = drop.tip(tree, tree$tip.label[!tree$tip.label %in% tip_labels])
  tree_f
}