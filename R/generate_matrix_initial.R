generate_matrix_initial <- function(rows, columns, connectance){
  # generate a feasible matrix of given number of rows, columns, and connectance
  if(is.null(rows)){
    rows <- 20
  }
  if(is.null(columns)){
    columns <- 20
  }
  if(is.null(connectance)){
    connectance <- 0.5
  }
  M <- matrix(0, nrow = rows, ncol = columns) # empty matrix
  nl <- round(prod(dim(M)) * connectance) # number of links
  if(nl <= rows + columns){stop("#Links <= rows + columns. This means the network contains more than one compartment (subsets of the web which are not connected to the main network). Please specify rows, columns and connectance so that #Links > rows + columns.")}
  M[1,] <- 1
  M[,1] <- 1
  M[2,2] <- 1
  nl <- nl - sum(M) # number of links left to allocate
  indices <- which(M == 0, arr.ind = TRUE) # indices of 0 elements of matrix i.e. indices that could have links assigned to them
  cl <- indices[sample(nrow(indices), size = nl, replace = FALSE),,drop = FALSE] # choose indices to allocate as links
  for(i in 1:nrow(cl)){ # fill in randomly chosen links
    cr <- cl[i,1]
    cc <- cl[i,2]
    M[cr,cc] <- 1
  }
  if(!feasibility(M)){stop("Initial matrix is not feasible. This should never happen. Contact the developer.")}
  connectance <- sum(M)/prod(dim(M)) # calculate exact connectance (we rounded earlier)
  return(list(M = M, connectance = connectance))
}
