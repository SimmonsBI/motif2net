#' @export
motif2net <- function(target_motif_distribution, distribution_type, six_node = FALSE, rows = NULL, columns = NULL, connectance = NULL, steps = 1000, iterations_at_each_temp = 1, alpha = 0.01){
  # argument checks
  if(!length(target_motif_distribution) %in% c(17,44)){stop("target_motif_distribution must be of length 17 (if only considering up to five node motifs) or 44 (if considering up to six node motifs)")}
  if(length(target_motif_distribution) == 17 & six_node == TRUE){stop("target_motif_distribution must be of length 44 if six_node = TRUE because there are 44 motifs up to six nodes. Either change the length of target_motif_distribution or set six_node = FALSE")}
  if(length(target_motif_distribution) == 44 & six_node == FALSE){stop("target_motif_distribution must be of length 17 if six_node = FALSE because there are 17 motifs up to five nodes. Either change the length of target_motif_distribution or set six_node = TRUE")}
  if(!distribution_type %in% c("frequency","normalise_sum","normalise_sizeclass","normalise_levelsize","normalise_nodesets")){stop("distribution type must be one of 'frequency', 'normalise_sum', 'normalise_sizeclass', 'normalise_levelsize' or 'normalise_nodesets'")}
  if(distribution_type == "frequency" & (!all(target_motif_distribution%%1==0))){stop("When distribution_type = 'frequency' all elements of target_motif_distribution must be integers")}
  if(distribution_type == "normalise_sum" & (sum(target_motif_distribution) != 1)){stop("When distribution_type = 'normalise_sum' elements of target_motif_distribution must sum to 1")}
  if(distribution_type == "normalise_sizeclass" & sizeclass_check(target_motif_distribution) == FALSE){stop("When distribution_type = 'normalise_sizeclass' elements of target_motif_distribution must sum to 1 within each motif size class")}
  if(distribution_type == "normalise_levelsize" & levelsize_check(target_motif_distribution) == FALSE){stop("When distribution_type = 'normalise_levelsize' elements of target_motif_distribution must sum to 1 within each motif levelsize class (a group of motifs with a given number of nodes in the top level and the bottom level)")}
  # run analysis
  M <- generate_matrix_initial(rows = rows, columns = columns, connectance = connectance)$M
  available_mutations <- 1:7
  unavailable_mutations <- NULL
  if(!is.null(rows)){unavailable_mutations <- c(unavailable_mutations, 4,6)}
  if(!is.null(columns)){unavailable_mutations <- c(unavailable_mutations, 5,7)}
  if(!is.null(connectance)){unavailable_mutations <- c(unavailable_mutations,1,2)}
  if(!is.null(unavailable_mutations)){available_mutations <- setdiff(available_mutations, unavailable_mutations)}
  sa(mat_initial = M, target_motif_distribution = target_motif_distribution, distribution_type = distribution_type, six_node = six_node, steps = steps, iterations_at_each_temp = iterations_at_each_temp, alpha = alpha, available_mutations = available_mutations, connectance = connectance)
}
