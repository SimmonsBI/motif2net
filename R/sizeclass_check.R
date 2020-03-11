sizeclass_check <- function(target_motif_distribution){
  if(length(target_motif_distribution) == 44){
    motif_sizes <- c(2,
                     3,3,
                     4,4,4,4,
                     rep(5,10),
                     rep(6, 27)
    )
    all(sapply(split(target_motif_distribution, motif_sizes), function(x) sum(x) == 1))
  } else if(length(target_motif_distribution) == 17){
    motif_sizes <- c(2,
                     3,3,
                     4,4,4,4,
                     rep(5,10)
    )
    all(sapply(split(target_motif_distribution, motif_sizes), function(x) sum(x) == 1))
  }
}
