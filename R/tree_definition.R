#' @title Tree generation
#' @description This function generates a Tree object, with the structure explined later, from a file.
#' First, it reads a file as a data frame with a specific format.
#' Then, it computes the tree structure.
#' @details The object Tree has the following structure:
#'
#' The tree structure is defined as a list composed by a set of features of the tree and a set of nodes.
#' The tree is stored in a file, and a function is provided to read the file and initialize all the data structure. Actually, the features stores in the structure are the following:
#'
#' Number of Nodes of the tree: tree$n_nodes.
#' Number of Primary events: tree$n_primary.
#' Number of Non-primary nodes in the tree, including intermediate nodes and the top event: tree$n_inter.
#' Identifier of the corner node: tree$corner.
#' Identifier of the top event node: tree$top.

#' @param tree_file is a string with the name of the file
#' @return A Tree object containing the tree strcucture
#' @examples Example of use: tree <- create_fault_tree("./tree.conf")
#'
#' @export create_fault_tree

# #############################################################################
# LOAD IN THE FAULT TREE STRUCTURE FROM THE SPECIFIED FILE
# The tree is specified in the file tree.conf and a tree object is created.
# #############################################################################
create_fault_tree <- function (tree_file) {

  data_tree<-read.table(tree_file, header = TRUE, sep = "", quote = "\"'", dec = ".",
                        col.names=c("Id", "Label", "Type", "Primary", "Logic", "Parent"),
                        as.is = FALSE, na.strings = "NA",
                        colClasses = NA, nrows = -1,
                        skip = 0, check.names = TRUE, fill = TRUE,
                        strip.white = FALSE, blank.lines.skip = TRUE,
                        comment.char = "#")

  nodes <- list()
  for (i in 1:nrow(data_tree)){
      nodes[[i]] <- data_tree[i,]
  }

  n_nodes=nrow(data_tree)
  n_primary=nrow(data_tree[data_tree$Type=="leaf" | data_tree$Type=="corner",])
  n_inter=nrow(data_tree[data_tree$Type=="inter" | data_tree$Type=="top",])
  corner=data_tree$Id[data_tree$Type=="corner"]
  top=data_tree$Id[data_tree$Type=="top"]

  children <- compute_children_adjacent_matrix (n_nodes, n_inter, n_primary, nodes)

  tree_definition <- list (n_nodes=nrow(data_tree), n_primary=nrow(data_tree[data_tree$Type=="leaf" | data_tree$Type=="corner",]),
                           n_inter=nrow(data_tree[data_tree$Type=="inter" | data_tree$Type=="top",]), corner=data_tree$Id[data_tree$Type=="corner"],
                          top=data_tree$Id[data_tree$Type=="top"], nodes=nodes, children=children)

  return (tree_definition)
}

# ##################################################################################
# This functions creates the children list using a Matrix of adjacent
# #################################################################################
compute_children_adjacent_matrix <- function (n_nodes, n_inter, n_primary, nodes) {

  child <- matrix(data=0,nrow=n_inter, ncol=(n_nodes-1), byrow=TRUE)
  for (s in 1:(n_nodes-1)) {
    child [nodes[[s]]$Parent - n_primary, s]<-1
  }

  return (child)
}
