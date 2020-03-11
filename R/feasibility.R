feasibility <- function(M){
  # tests if a matrix is feasible i.e. if all species have degree >= 1 and
  no_degree_zero <- !(any(colSums(M) == 0) | any(rowSums(M) == 0)) # do all species interact?
  one_compartment <- sum(M) > sum(dim(M)) # is the network a single connected compartment i.e. there are no isolated subgraphs
  if(no_degree_zero & one_compartment){
    TRUE
  } else {
  FALSE
  }
}
