# Sum-of-squares-based approach to endpoint 

library(tidyverse)
# construct operator "not in"
'%ni%' <- Negate('%in%')
# penalization constant
PENALIZATION <- 10^90

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

SSQ <- function(Bvec, to_fit_tidy, presence, use_logs){
  # make into matrix
  n <- sqrt(length(Bvec))
  B <- matrix(Bvec, n, n)
  # for each community predict coexistence
  # if no coexistence penalize by setting all density to arbitrarily large value
  for (i in 1:nrow(presence)){
    present_species <- presence[i, 1:n] > 0
    tryCatch(
      {
        xstar <- -rowSums(solve(B[present_species, present_species, drop = FALSE]))
      },
      # ... we've encountered a numerically singular matrix: 
      error=function(error_message) {
        xstar <<- rep(-1, sum(present_species))
      }
    )
    if (any(xstar < 0)){
      xstar[] <- PENALIZATION
    } 
    presence[i, c(present_species, FALSE)] <- as.list(xstar)
  }
  for_ssq <- to_fit_tidy %>% left_join(presence %>% gather("species", "predicted", -community) %>% filter(predicted > 0),
                                       by = c("community", "species"))
  if (use_logs == TRUE) return(sum((log(for_ssq$abundance) - log(for_ssq$predicted))^2))
  return(sum((for_ssq$abundance - for_ssq$predicted)^2))
}

fit_endpoints <- function(dt, # the data to fit
                          exclude = NULL, # vector containing the endpoints to exclude
                          use_logs = FALSE, # use logs for SSQ
                          initial_matrix = NULL, # if not specified, use -I
                          method = "BFGS", # parameters to pass to optim e.g. method = "Nelder-Mead"
                          num_steps = 5000,
                          show_progress = FALSE){
  to_fit <- dt %>% filter(community %ni% exclude) 
  n <- ncol(to_fit) - 1
  to_fit_tidy <- to_fit %>% gather("species", "abundance", -community) %>% filter(abundance > 0)
  # now extract the communities to fit, using just presence/absence
  presence <- to_fit %>% group_by(community) %>% filter(row_number()==1) %>% ungroup()
  presence <- presence %>% mutate_at(vars(-community), function(x) (x > 0) * 1) %>% distinct()
  # initial B:
  # either provided by the user 
  # or taken to be -I
  B <- matrix(0, n, n)
  diag(B) <- -1
  if (!is.null(initial_matrix)) B <- initial_matrix
  
  Bvec <- as.vector(B)
  current_ssq <- SSQ(Bvec, to_fit_tidy, presence, use_logs)
  print("initial B")
  print(B)
  print("initial SSQ")
  print(current_ssq)
  
  tmp <- optim(Bvec, SSQ, to_fit_tidy = to_fit_tidy, presence = presence, use_logs = use_logs,
               method = method,
               control = list(maxit = num_steps, trace = show_progress))
  
  B <- matrix(tmp$par, n, n)
  final_ssq <- tmp$value
  print("final B")
  print(B)
  print("final SSQ")
  print(final_ssq)
  
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
     predicted[i, c(present_species, FALSE)] <- as.list(xstar)
     if (any(xstar < 0)){
       predicted[i, c(present_species, FALSE)] <- as.list(rep(-1, sum(present_species)))  
     }
   }
   # now organize into observed and predicted
   results <- to_predict %>% gather("species", "observed", -community) %>% filter(observed != 0) %>% 
     left_join(predicted %>% gather("species", "predicted", -community) %>% filter(predicted != 0), 
               by = c("community", "species"))
   return(list(B = B, results = results))
}

predict_out_of_fit <- function(dt,
                               use_logs = FALSE, # use logs for SSQ
                               initial_matrix = NULL, # if not specified, use -I
                               method = "BFGS", # parameters to pass to optim e.g. method = "Nelder-Mead"
                               num_steps = 5000,
                               show_progress = FALSE
                               ){
  # This function predicts each endpoint (with all its replicates) out of fit.
  # The function excludes each endpoint in turn, and use the remaining data to predict it out of fit
  results <- tibble()
  for (i in sort(unique(dt$community))){
    print(i)
    results <- rbind(results, fit_endpoints(dt, exclude = i, 
                                            use_logs = use_logs,
                                            initial_matrix = initial_matrix,
                                            method = method,
                                            show_progress = show_progress)$results)
  }
  return(results)
}

stats_and_plot_results <- function(results){
  print("Correlation")
  print(cor(results$observed, results$predicted))
  pl <- ggplot(data = results) + aes(x = observed, y = predicted, colour = species) + geom_point() + geom_abline(slope = 1, intercept = 0)
  return(pl)
}

