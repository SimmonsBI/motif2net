acceptance_probability <- function(cost_current, cost_neighbour, Temp){
  if(Temp == 0){stop("Temp is 0 which means the acceptance probability cannot be calculated. Choose a slower cooling schedule by increasing 'alpha' or reducing 'steps'")}
  exp(-(cost_neighbour - cost_current)/Temp)
}
