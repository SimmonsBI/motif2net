cost_maxmotif <- function(target_motif, solution, six_node){
  mcount_out <- mcount(M = solution, six_node = six_node, normalisation = FALSE, mean_weight = FALSE, standard_dev = FALSE)
  -mcount_out[mcount_out$motif == target_motif, "frequency"] # cost = -frequency to ensure frequency is maximised
}
