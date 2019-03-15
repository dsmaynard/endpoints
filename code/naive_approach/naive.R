# Naive regression-based approach to endpoint 

library(tidyverse)
# construct operator "not in"
'%ni%' <- Negate('%in%')


prepare_data <- function(filename = "../../data/Kuebbing_plants/natives.csv"){
  # Input
  # filename: file name of a csv file containing the endpoints (rows are observations, columns are species, headers are species names)
  # Output
  # a tibble as the csv file but with an extra column containing the label for the community
  dt <- read_csv(filename) 
  # now label communities
  spnames <- colnames(dt)
  dt <- dt %>% add_column(community = apply(dt, 1, function(x) paste0(spnames[x > 0], collapse = "-")))
  return(dt)
}

fit_endpoints <- function(dt, exclude = NULL){
  # compute the endpoints of all possible endpoints
  # given a matrix of observed endpoints
  # Input
  # dt: tibble produced calling prepare_data()
  # exclude: vector containing the communities to exclude from the fit
  # Output
  # B: fitted matrix
  # fitted: tibble containing, for each community and species, the observed vs predicted values
  #         if exclude == NULL return all communities
  #         otherwise, only the out-of-fit predictions for the excluded communities
  E <- as.matrix(dt %>% filter(community %ni% exclude) %>% select(-community))
  n <- ncol(E)
  B <- matrix(0, n, n)
  for(i in 1:n){
    Ei <- E[E[,i] > 0,, drop = FALSE]
    Bi <- -rowSums(MASS::ginv(Ei))
    B[i,] <- Bi
  }
  # now compute predicted endpoints
  to_predict <- dt
  if (!is.null(exclude)){
    to_predict <- dt %>% filter(community %in% exclude)
  }
  predicted <- to_predict
  predicted <- predicted %>% mutate_at(vars(-community), function(x) (x > 0) * 1) %>% distinct()
  for (i in 1:nrow(predicted)){
    present_species <- predicted[i, 1:n] > 0
    xstar <- -rowSums(solve(B[present_species, present_species, drop = FALSE]))
    predicted[i, c(present_species, FALSE)] <- xstar
    if (any(xstar < 0)){
      predicted[i, c(present_species, FALSE)] <- -1  
    }
  }
  # now organize into observed and predicted
  results <- to_predict %>% gather("species", "observed", -community) %>% filter(observed != 0) %>% 
    left_join(predicted %>% gather("species", "predicted", -community) %>% filter(predicted != 0), 
              by = c("community", "species"))
  return(list(B = B, results = results))
}

predict_out_of_fit <- function(dt){
  # This function predicts each endpoint (with all its replicates) out of fit.
  # The function excludes each endpoint in turn, and use the remaining data to predict it out of fit
  results <- tibble()
  for (i in sort(unique(dt$community))){
    results <- rbind(results, fit_endpoints(dt, exclude = i)$results)
  }
  return(results)
}

stats_and_plot_results <- function(results){
  print("Correlation")
  print(cor(results$observed, results$predicted))
  pl <- ggplot(data = results) + aes(x = observed, y = predicted, colour = species) + geom_point() + geom_abline(slope = 1, intercept = 0)
  return(pl)
}

