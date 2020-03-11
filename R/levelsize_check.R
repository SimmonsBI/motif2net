levelsize_check <- function(target_motif_distribution){
  if(length(target_motif_distribution) == 44){
    motif_sizes <- c(1,
                     2,
                     3,
                     4,
                     5,5,
                     7,
                     8,
                     rep(9,4),
                     rep(10,4),
                     11,
                     12,
                     rep(13,6),
                     rep(14,13),
                     rep(15,6),
                     44
    )
    all(sapply(split(target_motif_distribution, motif_sizes), function(x) sum(x) == 1))
  } else if(length(target_motif_distribution) == 17){
    motif_sizes <- c(1,
                     2,
                     3,
                     4,
                     5,5,
                     7,
                     8,
                     rep(9,4),
                     rep(10,4),
                     11
    )
    all(sapply(split(target_motif_distribution, motif_sizes), function(x) sum(x) == 1))
  }
}
