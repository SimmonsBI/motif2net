neighbour <- function(M, available_mutations, connectance){
  # generates a feasible neighbour solution by either adding a link, removing a link or moving a link
  if(is.null(connectance)){connectance <- sum(M)/prod(dim(M))} # if no specified connectance, adding nodes or removing nodes will maintain the current connectance
  feasible <- FALSE
  Mo <- M # store original matrix
  full <- all(M == 1)
  empty <- all(M == 0)
  if(empty == TRUE){stop("Matrix is empty. This should never happen. Contact the developer.")}
  if(full == TRUE){
    available_mutations <- available_mutations[-which(available_mutations %in% c(1,3))] # is matrix is full, cannot add or move a link
  }
  mutation <- sample_vl(x = available_mutations, size = 1) # pick a random mutation
  while_loop_countdown <- 100 # allow 100 attempts at a given mutation before picking a new one
  while(feasible == FALSE){
    if(while_loop_countdown == 0){
      mutation <- sample_vl(x = available_mutations, size = 1) # pick a new random mutation
      while_loop_countdown <- 100 # reset the countdown
    }
    M <- Mo # reset matrix to original
    if(mutation == 1){ # adds a link
      to_add <- which(M == 0, arr.ind = TRUE) # coordinates of 0s
      wc <- sample.int(nrow(to_add), 1) # pick a row from list of 0 coordinates
      cr <- to_add[wc,"row"] # extract the chosen row
      cc <- to_add[wc,"col"]# extract the chosen column
      M[cr,cc] <- 1 # add a link
    } else if(mutation == 2){ # removes a link
      to_remove <- which(M != 0, arr.ind = TRUE) # coordinates of 1s
      wc <- sample.int(nrow(to_remove), 1) # pick a row from list of 1 coordinates
      cr <- to_remove[wc,"row"] # extract the chosen row
      cc <- to_remove[wc,"col"] # extract the chosen column
      M[cr,cc] <- 0 # remove the link
    } else if(mutation == 3){ # moves one link e.g. i-j -> m-n
      to_remove <- which(M != 0, arr.ind = TRUE) # coordinates of 1s
      to_add <- which(M == 0, arr.ind = TRUE) # coordinates of 0s
      wc <- sample.int(nrow(to_remove), 1) # pick a row from list of 1 coordinates
      cr <- to_remove[wc,"row"] # extract the chosen row
      cc <- to_remove[wc,"col"] # extract the chosen column
      M[cr,cc] <- 0 # remove the link

      wc <- sample.int(nrow(to_add), 1) # pick a row from list of 0 coordinates
      cr <- to_add[wc,"row"] # extract the chosen row
      cc <- to_add[wc,"col"]# extract the chosen column
      M[cr,cc] <- 1 # add a link
    } else if(mutation == 4){ # adds a row node while maintaining connectance
      M <- rbind(M, rep(0, ncol(M)))
      connectance_after_mutation <- sum(M)/prod(dim(M))
      nlinks_to_maintain_connectance <- round((connectance - connectance_after_mutation)*prod(dim(M)))
      if(ncol(M) < nlinks_to_maintain_connectance){stop("ncol(M) < nlinks_to_maintain_connectance, cannot maintain connectance")}
      if(nlinks_to_maintain_connectance > 0){
        M[nrow(M),sample_vl(1:ncol(M),size = nlinks_to_maintain_connectance)] <- 1
      }
    } else if(mutation == 5){ # add a column node while maintaining connectance
      M <- cbind(M, rep(0, nrow(M)))
      connectance_after_mutation <- sum(M)/prod(dim(M))
      nlinks_to_maintain_connectance <- round((connectance - connectance_after_mutation)*prod(dim(M)))
      if(nrow(M) < nlinks_to_maintain_connectance){stop("nrow(M) < nlinks_to_maintain_connectance, cannot maintain connectance")}
      if(nlinks_to_maintain_connectance > 0){
        M[sample_vl(1:nrow(M),size = nlinks_to_maintain_connectance),ncol(M)] <- 1
      }
    } else if(mutation == 6){ # remove a row node while maintaining connectance
      M <- M[-sample_vl(1:nrow(M), size = 1),,drop = FALSE]
      connectance_after_mutation <- sum(M)/prod(dim(M))
      nlinks_to_maintain_connectance <- round((connectance - connectance_after_mutation)*prod(dim(M)))
      if(nlinks_to_maintain_connectance > 0){
        indices <- which(M == 0, arr.ind = TRUE) # indices of 0 elements of matrix i.e. indices that could have links assigned to them
        cl <- indices[sample_vl(nrow(indices), size = nlinks_to_maintain_connectance, replace = FALSE),,drop = FALSE] # choose indices to allocate as links
        for(i in 1:nrow(cl)){ # fill in randomly chosen links
          cr <- cl[i,1]
          cc <- cl[i,2]
          M[cr,cc] <- 1
        }
      } else if (nlinks_to_maintain_connectance < 0){
        indices <- which(M == 1, arr.ind = TRUE) # indices of 0 elements of matrix i.e. indices that could have links assigned to them
        cl <- indices[sample_vl(nrow(indices), size = abs(nlinks_to_maintain_connectance), replace = FALSE),,drop = FALSE] # choose indices to allocate as links
        for(i in 1:nrow(cl)){ # fill in randomly chosen links
          cr <- cl[i,1]
          cc <- cl[i,2]
          M[cr,cc] <- 0
        }
      }
    } else if(mutation == 7){ # remove a column node while maintaining connectance
      M <- M[,-sample_vl(1:ncol(M), size = 1), drop = FALSE]
      connectance_after_mutation <- sum(M)/prod(dim(M))
      nlinks_to_maintain_connectance <- round((connectance - connectance_after_mutation)*prod(dim(M)))
      if(nlinks_to_maintain_connectance > 0){
        indices <- which(M == 0, arr.ind = TRUE) # indices of 0 elements of matrix i.e. indices that could have links assigned to them
        cl <- indices[sample_vl(nrow(indices), size = nlinks_to_maintain_connectance, replace = FALSE),,drop = FALSE] # choose indices to allocate as links
        for(i in 1:nrow(cl)){ # fill in randomly chosen links
          cr <- cl[i,1]
          cc <- cl[i,2]
          M[cr,cc] <- 1
        }
      } else if (nlinks_to_maintain_connectance < 0){
        indices <- which(M == 1, arr.ind = TRUE) # indices of 0 elements of matrix i.e. indices that could have links assigned to them
        cl <- indices[sample_vl(nrow(indices), size = abs(nlinks_to_maintain_connectance), replace = FALSE),,drop = FALSE] # choose indices to allocate as links
        for(i in 1:nrow(cl)){ # fill in randomly chosen links
          cr <- cl[i,1]
          cc <- cl[i,2]
          M[cr,cc] <- 0
        }
      }
    }
    while_loop_countdown <- while_loop_countdown - 1
    if(feasibility(M) == TRUE){feasible <- TRUE}
  }
  M # return the matrix
}
