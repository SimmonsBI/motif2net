sa <- function(mat_initial, target_motif_distribution, distribution_type, six_node, steps = 1000, iterations_at_each_temp = 1, alpha = 0.01, available_mutations, connectance){
  # initialise matrices
  mat_optimum <- mat_initial
  mat_current <- mat_initial
  mat_neighbour <-  mat_initial
  # initialise costs
  cost_initial <- cost(target_motif_distribution = target_motif_distribution, solution = mat_initial, six_node = six_node, distribution_type = distribution_type)
  cost_optimum <- cost_initial
  cost_current <- cost_initial
  cost_neighbour <-  cost_initial
  # initialise temperature
  Temp <- 1
  # begin optimisation
  tick <- 0
  record <- NULL # remove when finished testing
  for(i in 1:steps){
    for(j in 1:iterations_at_each_temp){
      tick <- tick + 1
      # generate neighbour matrix and cost
      mat_neighbour <- neighbour(M = mat_current, available_mutations = available_mutations, connectance = connectance)
      cost_neighbour <- cost(target_motif_distribution = target_motif_distribution, solution = mat_neighbour, six_node = six_node, distribution_type = distribution_type)
      # upate current matrix and cost
      if(cost_neighbour < cost_current){
        mat_current <- mat_neighbour
        cost_current <- cost_neighbour
      } else if(cost_neighbour >= cost_current){
        record <- rbind(record, c(acceptance_probability(cost_current, cost_neighbour, Temp), i)) # remove when finished testing
        if(acceptance_probability(cost_current, cost_neighbour, Temp) > runif(1,0,1)){
          mat_current <- mat_neighbour
          cost_current <- cost_neighbour
        }
      }
      # update optimum matrix and cost
      if(cost_neighbour < cost_optimum){
        mat_optimum <- mat_neighbour
        cost_optimum <- cost_neighbour
      }
      if(tick%%10==0){
        cat("\r", "Progress = ",substr(x = sprintf("%.2f", (i/steps)*100), 1, 4), "% | Optimum = ",substr(x = sprintf("%.5f", cost_optimum), 1, 7), " | Current = ",substr(x = sprintf("%.5f", cost_current), 1, 7), " | Neighbour = ",substr(x = sprintf("%.5f", cost_neighbour), 1, 7)," | Temp = ", substr(x = sprintf("%.5f", Temp), 1, 7)," | size ",substr(x = sprintf("%.1f", sum(dim(mat_neighbour))), 1, 4), sep = "")
      }
    }
    Temp <- Temp*alpha # cool the system
  }
  # assemble output
  optimum_motif_distribution <- bmotif::mcount(M = mat_optimum, six_node = six_node, normalisation = TRUE, mean_weight = FALSE, standard_dev = FALSE)
  motif_counts <- rbind(target_motif_distribution, optimum_motif_distribution[,distribution_type])
  rownames(motif_counts) <- c("target_motif_distribution","optimum_motif_distribution")
  return(list(optimum_cost = cost_optimum, optimum_matrix = mat_optimum, motif_counts = motif_counts, record = record))
}

# optional arguments for how feasible to make it?
# add rows/columns in neighbour solutions?
