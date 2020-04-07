#' @export
max_motif <- function(target_motif, rows, columns, connectance = NULL, steps = 1000, iterations_at_each_temp = 1, alpha = 0.99){
  # argument checks
  if(!inherits(target_motif, "numeric")){stop("'target_motif' must be a motif ID between 1 and 44 following Figure 3 in Simmons et al (2019), Oikos, 128: 154-170. doi:10.1111/oik.05670")}
  if(!target_motif %in% 1:44){stop("'target_motif' must be a motif ID between 1 and 44 following Figure 3 in Simmons et al (2019), Oikos, 128: 154-170. doi:10.1111/oik.05670")}
  six_node <- ifelse(target_motif >= 18, yes = TRUE, no = FALSE)
  if(!inherits(rows, "numeric")){stop("'rows' must be numeric")}
  if(!inherits(columns, "numeric")){stop("'columns' must be numeric")}
  if(!is.null(connectance)){
    if(!inherits(connectance, "numeric")){stop("'connectance' must be numeric")}
  }
  if(!inherits(steps, "numeric")){stop("'steps' must be numeric")}
  if(!inherits(iterations_at_each_temp, "numeric")){stop("'iterations_at_each_temp' must be numeric")}
  if(!inherits(alpha, "numeric")){stop("'alpha' must be numeric")}
  # run analysis
  M <- generate_matrix_initial(rows = rows, columns = columns, connectance = connectance)$M
  available_mutations <- 1:3
  if(!is.null(connectance)){unavailable_mutations <- c(1,2)}
  if(!is.null(unavailable_mutations)){available_mutations <- setdiff(available_mutations, unavailable_mutations)}
  sa_maxmotif(mat_initial = M, target_motif = target_motif, six_node = six_node, steps = steps, iterations_at_each_temp = iterations_at_each_temp, alpha = alpha, available_mutations = available_mutations, connectance = connectance)
}
