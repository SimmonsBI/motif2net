cost <- function(target_motif_distribution, solution, six_node, distribution_type){
  mcount_out <- mcount(M = solution, six_node = six_node, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)
  dm <- rbind(target_motif_distribution, mcount_out[,distribution_type])
  vegdist(dm, method = "bray")[1]
}
